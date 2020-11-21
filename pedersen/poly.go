package pedersen

import (
	"crypto/cipher"
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"math"

	"github.com/drand/kyber"
	"github.com/drand/kyber/util/random"
	"github.com/nikkolasg/playsnark/pedersen/commit"
	poly "github.com/nikkolasg/playsnark/polynomial"
	"golang.org/x/crypto/blake2s"
)

// Group has points on it and can create scalar from the scalar fields
type Group = kyber.Group

// Scalar of the field of the curve
type Scalar = kyber.Scalar

// Point in the group (in our case it's elliptic curve so it's a point)
type Point = kyber.Point

type Setup struct {
	commit.Setup
	g      Group
	h      Point
	degree int
}

// NewSetup returns the information necessary for using the pedersen polynomail
// commitment
func NewSetup(g Group, degree int) Setup {
	s := commit.NewSetup(g, degree+1)
	fmt.Println("SO FAR SO GOOD")
	h := g.Point().Pick(oracle(g, degree))
	return Setup{
		Setup:  s,
		h:      h,
		degree: degree,
		g:      g,
	}
}

type Polynomial = []Scalar

type Commitment = Point

// Commit the pedersen commitment of the coefficients
func Commit(s Setup, p poly.Poly, w Scalar) Commitment {
	return commit.Commit(s.Setup, p.Coeffs(), w)
}

// Open strictly follows the algorithm description in Appendix A at
// https://eprint.iacr.org/2020/499.pdf
func Open(s Setup, p poly.Poly, c Commitment, z, w Scalar) Proof {
	// 1. v = p(z)
	v := p.Eval(z)
	// 2. sample random poly p2 such that  p2(z) = 0
	// naive strategy: create a random one, then eval p2(z) = a then substract
	// a from the free coefficient to p2' such  that p2'(z) = p2(z) - a = 0
	p2 := poly.RandomPoly(s.g, s.degree)
	a := p2.Eval(z)
	p2.Set(0, p2.Coeffs()[0].Sub(p2.Coeffs()[0], a))
	assert(p2.Eval(z).Equal(s.g.Scalar().Zero()), "p2(z) != 0")
	// 3. sample commitment randomness w2
	w2 := s.g.Scalar().Pick(random.New())
	// 4. commit to p2 using w2 - hiding commitment
	c2 := commit.Commit(s.Setup, p2.Coeffs(), w2)
	// 5. compute challenge from random oracle
	alpha := s.g.Scalar().Pick(oracleFrom(s.g, c, z, v, c2))
	// 6. then compute p3 = p + alpha * p2
	p2.Scale(alpha)
	p3 := p.Add(p2)
	// p3(z) = p(z) + alpha * p2(z) = v + 0
	assert(p3.Eval(z).Equal(v), "p3 not correct")
	// 7. compute corresponding randomness w3 = w + alpha * w2
	w3 := s.g.Scalar()
	w3 = w3.Mul(alpha, w2)
	w3 = w3.Add(w3, w)
	// 8. compute non hiding commitment to p3 without computing it from scratch
	// Commit(p3) = SUM p3_i * Gi = SUM (p_i + alpha * p2_i) * Gi
	// C + alpha * C2 - w3 * S =
	// w*S + SUM p_i*Gi + alpha * (w2*S + SUM p2_i*Gi) - wS - alpha*w2*S =
	// SUM (p_i + alpha * p2_i) * Gi
	ws := s.g.Point().Mul(w3, s.Setup.S)
	alphac2 := s.g.Point().Mul(alpha, c2)
	c3 := s.g.Point().Add(c, alphac2)
	c3 = c3.Add(c3, ws.Neg(ws))
	assert(c3.Equal(commit.Commit(s.Setup, p3.Coeffs(), nil)), "c3 not correct")
	// --------------
	// Start computing challenges for the recursion and compute new random base
	e0 := s.g.Scalar().Pick(oracleFrom(s.g, c3, z, v))
	h2 := s.g.Point().Mul(e0, s.h)
	fmt.Println("PROOF H2 : ", h2)
	// initialize the vectors for the recursion
	// coeffcients of p3
	cv := p3.Coeffs()
	// powers of z: 1, z, z'2 ...
	zv := make([]Scalar, len(cv))
	acc := s.g.Scalar().One()
	zv[0] = acc.Clone()
	for i := 1; i < len(zv); i++ {
		zv[i] = acc.Mul(acc, z)
	}
	// base points of the pedersen commitment
	gv := s.Setup.Gs
	assert(len(gv) == len(cv), "gs incompatible length")
	ei := e0
	// i: 1 -> log(degree+1) as we are doing a binary dissection of the three
	// precedent vectors -- We assume length is power of two
	max := int(math.Ceil(math.Log2(float64(p.Degree() + 1))))
	// slice of vector that will be filled at each iteration and returned in
	// proof - they're the pedersen commitments of left and right part at each
	// iteration
	var lis []Point
	var ris []Point
	for i := 1; i <= max; i++ {
		// 1.
		// append H2 to the left part of commitment basis for the pedersen
		// commitments for the left part
		l := append(cloneP(leftP(gv)), h2)
		// inner product between the right part of c3 coefficients (polynomial
		// created using the randomness) and left part of powers of z
		leftCz := innerProd(right(cv), left(zv))
		leftSetup := commit.Setup{Gs: l, S: s.Setup.S, G: s.g}
		li := commit.Commit(leftSetup, append(cloneS(right(cv)), leftCz), nil)
		lis = append(lis, li)
		// 2.
		// inner product between the left part of c3 coefficients and right part
		// of power z
		rightCz := innerProd(left(cv), right(zv))
		// commitment basis for the right pedersen commitments
		r := append(cloneP(rightP(gv)), h2)
		rightSetup := commit.Setup{Gs: r, S: s.Setup.S, G: s.g}
		ri := commit.Commit(rightSetup, append(cloneS(left(cv)), rightCz), nil)
		ris = append(ris, ri)
		// 3. get oracle from these commits
		ei = s.g.Scalar().Pick(oracleFrom(s.g, ei, li, ri))
		fmt.Println("PROOF: i", i, " ei = ", ei)
		// 4. construct commitment key for next round:
		// left part of gv + ei * right part of gv
		pgv := cloneP(rightP(gv))
		for i := 0; i < len(pgv); i++ {
			pgv[i] = pgv[i].Mul(ei, pgv[i])
		}
		gv = sumVecP(cloneP(leftP(gv)), pgv)
		// 5. and inputs for next round cv and zv
		// scale right part of cv by negative ei
		pcv := cloneS(right(cv))
		scale(pcv, s.g.Scalar().Neg(ei))
		// left(ci) + ei * right(ci)
		cv = sumVecS(cloneS(left(cv)), pcv)
		// 6. same for z
		pzv := cloneS(right(zv))
		scale(pzv, ei)
		zv = sumVecS(cloneS(left(zv)), pzv)
	}
	assert(len(cv) == 1, "binary dissection is off - len %d", len(cv))
	assert(len(gv) == 1, "binary dissection is off - len %d", len(gv))
	return Proof{
		U:  gv[0],
		C:  cv[0],
		C2: c2,
		W:  w3,
		Ri: ris,
		Li: lis,
	}
}

type Proof struct {
	U Point
	C Scalar
	// Hiding commitment
	C2 Point
	// Commitment randomness
	W  Scalar
	Ri []Point
	Li []Point
}

// Verify takes the setup, c which is the commitment of the polynomial, the
// evaluation point z, the alleged evaluation v, and a proof p of opening.
// It returns true if the proof is correct, false otherwise. The check performs
// the verification routine outlined in the Check section, appendix A.2 of the
// paper.
func Verify(s Setup, C Commitment, z, v Scalar, p Proof) bool {
	// 1. 2.
	d2 := len(s.Setup.Gs) - 1
	// 3. 4. - this is the log(n) check
	h, ok := SuccinctCheck(s.g, s.Setup.S, s.h, d2, C, s.degree, z, v, p)
	if !ok {
		return false
	}
	// 5. this is the linear check as commiting is O(n)
	exp := commit.Commit(s.Setup, h.Coeffs(), nil)
	if !exp.Equal(p.U) {
		return false
	}
	return true
}

func SuccinctCheck(g Group, S, H Point, d2 int, C Point, d int, z, v Scalar, p Proof) (*poly.Poly, bool) {
	// 2.
	if d != d2 {
		fmt.Println("YOLO")
		return nil, false
	}
	// 3. compute challenge alpha as in the proving part
	alpha := g.Scalar().Pick(oracleFrom(g, C, z, v, p.C2))
	// 4. compute non hiding commitment as in the proving part
	ws := g.Point().Mul(p.W, S)
	alphac2 := g.Point().Mul(alpha, p.C2)
	c3 := g.Point().Add(C, alphac2)
	c3 = c3.Add(c3, ws.Neg(ws))
	// 5. Compute first challenge and set new base for H
	max := int(math.Ceil(math.Log2(float64(d + 1))))
	// keep track of all eis generated for later use to construct the verifying
	// polynomial
	var eis = make([]Scalar, 0, max)
	ei := g.Scalar().Pick(oracleFrom(g, c3, z, v))
	H2 := g.Point().Mul(ei, H)
	fmt.Println("CHECK H2: ", H2)
	// 6. Compute first group element C0 = C3 + v * H2
	ci := g.Point().Add(c3, g.Point().Mul(v, H2))
	// 7. Recursion step
	for i := 1; i <= max; i++ {
		// (a) Generate new challenge from previous challenge and L_i and R_i values
		// as in the proof
		li := p.Li[i-1] // we are offsetting by one since there is no R0 or  L0
		ri := p.Ri[i-1]
		ei = g.Scalar().Pick(oracleFrom(g, ei, li, ri))
		fmt.Println("CHECK: i", i, " ei = ", ei)
		eis = append(eis, ei)
		// (b) compute new commitment
		// ei^-1 * Li + C_i-1 + ei * Ri
		nei := g.Scalar().Neg(ei)
		neili := g.Point().Mul(nei, li)
		eiri := g.Point().Mul(ei, ri)
		ci = g.Point().Add(ci, g.Point().Add(neili, eiri))
	}
	// 8. define polynomial h(x) with eis and evaluate it to z. We construct
	// each term as polynomial and multiply them
	// set up an accumulator polynomial init to p(x) = 1 to multiply along the
	// iteration
	// SUM (i:0 -> log(d+1)-1) (1 + ei * X^2^i) =
	// (1 + ei * X) * (...) * (1 + ei * X^2^log(d+1)-1)
	acc := poly.NewZeroPoly(g, 0)
	acc.Set(0, g.Scalar().SetInt64(1))
	for i := 0; i < max; i++ {
		power := 1 << i
		// 1 + e_(log(d+1) - i) * X^2^i
		p := poly.NewZeroPoly(g, power)
		fmt.Printf("constructing power %d -> len(coeffs) = %d\n", power, len(p.Coeffs()))
		p.Set(power, eis[max-1-i]) //  -1 because of array offset
		p.Set(0, g.Scalar().SetInt64(1))
		acc = acc.Mul(p)
	}
	fmt.Println("final poly has degree ", acc.Degree())
	// 9. eval v' = acc(z) * c
	v2 := g.Scalar().Mul(p.C, acc.Eval(z))
	// 10. Compute pedersen commitment over c || v' using special basis
	// i.e. c * U + c * h(z) * H2
	setup := commit.Setup{Gs: []Point{p.U, H2}, G: g, S: S}
	computed := commit.Commit(setup, []Scalar{p.C, v2}, nil)
	expected := ci
	if !computed.Equal(expected) {
		fmt.Println("HERE IT?S GONE WRONG")
		return nil, false
	}
	return &acc, true
}

func left(a []Scalar) []Scalar {
	return a[:half(len(a))]
}

func right(a []Scalar) []Scalar {
	return a[half(len(a)):]
}

func leftP(p []Point) []Point {
	return p[:half(len(p))]
}
func rightP(p []Point) []Point {
	return p[half(len(p)):]
}

func cloneP(a []Point) []Point {
	b := make([]Point, len(a))
	copy(b, a)
	return b
}

func cloneS(s []Scalar) []Scalar {
	b := make([]Scalar, len(s))
	copy(b, s)
	return b
}

func scale(s []Scalar, f Scalar) {
	for i := range s {
		s[i] = s[i].Mul(s[i], f)
	}
}

// half assumes length is a power of two
func half(length int) int {
	return length / 2
}

func innerProd(a, b []Scalar) Scalar {
	assert(len(a) == len(b), "inner prod not possible")
	sum := a[0].Clone().Zero()
	tmp := a[0].Clone()
	for i := 0; i < len(a); i++ {
		sum = sum.Add(sum, tmp.Mul(a[i], b[i]))
	}
	return sum
}

// IN PLACE
// Damn generics in Go ...
func sumVecS(a, b []Scalar) []Scalar {
	assert(len(a) == len(b), "sumVecS not possible")
	for i := range a {
		a[i] = a[i].Add(a[i], b[i])
	}
	return a
}

func sumVecP(a, b []Point) []Point {
	assert(len(a) == len(b), "sumVecP not possible")
	for i := range a {
		a[i] = a[i].Add(a[i], b[i])
	}
	return a
}

func assert(cond bool, msg string, vars ...interface{}) {
	if !cond {
		panic(fmt.Sprintf(msg, vars...))
	}
}

// returns an oracle for deriving a point for pedersen commitment
func oracle(g Group, l int) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-pcommit"))
	h.Write([]byte(g.String()))
	binary.Write(h, binary.LittleEndian, int32(l))
	xof, err := blake2s.NewXOF(0, h.Sum(nil))
	assert(err == nil, "failed oracle")
	return random.New(xof)
}

func oracleFrom(g Group, infos ...kyber.Marshaling) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-polycommit-scalar"))
	h.Write([]byte(g.String()))
	for _, s := range infos {
		s.MarshalTo(h)
	}
	xof, err := blake2s.NewXOF(0, h.Sum(nil))
	assert(err == nil, "failed oracle")
	return random.New(xof)
}
