package pedersen

import (
	"bytes"
	"crypto/cipher"
	"crypto/sha256"
	"encoding/binary"
	"math"

	"github.com/drand/kyber"
	"github.com/drand/kyber/util/random"
	"github.com/nikkolasg/playsnark/pedersen/commit"
	poly "github.com/nikkolasg/playsnark/polynomial"
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
func Commit(s Setup, p poly.Poly, r Scalar) Commitment {
	return commit.Commit(s.Setup, p.Coeffs(), r)
}

// Open strictly follows the algorithm description in Appendix A at
// https://eprint.iacr.org/2020/499.pdf
func Open(s Setup, p poly.Poly, c Commitment, z, w Scalar) Proof {
	v := p.Eval(z)
	// random poly p2 such that  p2(z) = 0
	// naive strategy: create a random one, then eval p2(z) = a then substract
	// a from the free coefficient to p2' such  that p2'(z) = p2(z) - a = 0
	p2 := poly.RandomPoly(s.g, s.degree)
	a := p2.Eval(z)
	p2.Set(0, p2.Coeffs()[0].Sub(p2.Coeffs()[0], a))
	// sample randomness w
	w2 := s.g.Scalar().Pick(random.New())
	// commit to p2 using w - hiding commitment
	c2 := commit.Commit(s.Setup, p2.Coeffs(), w)
	// compute challenge from random oracle
	alpha := s.g.Scalar().Pick(oracleFrom(s.g, c, z, v, c2))
	// then compute p3 = p + alpha * p2
	p2.Scale(alpha)
	p3 := p.Add(p2)
	// compute corresponding randomness w3 = w + alpha * w2
	w3 := s.g.Scalar()
	w3.Mul(alpha, w2)
	w3.Add(w3, w)
	// compute non hiding commitment, i.e. we remove the first hiding of the
	// commitment C
	// c3 = c + alpha * c2 - w * S where w and S are respectively the randomness
	// and the point used for making c in the first place -> we now have a
	// commitment on p3 !
	ws := s.g.Point().Mul(w, s.Setup.S)
	alphac2 := s.g.Point().Mul(alpha, c2)
	c3 := s.g.Point().Add(c, alphac2)
	c3.Add(c3, ws.Neg(ws))
	// Start computing challenges for the recursion and compute new random base
	e0 := s.g.Scalar().Pick(oracleFrom(s.g, c3, z, v))
	h2 := s.g.Point().Mul(e0, s.h)
	// initialize the vectors for the recursion
	// coeffcients of p3
	cv := p3.Coeffs()
	// powers of z: 1, z, z'2 ...
	zv := make([]Scalar, len(cv))
	acc := s.g.Scalar().One()
	for i := 1; i < len(zv); i++ {
		zv[i] = acc.Mul(acc, z)
	}
	// base points of the pedersen commitment
	gv := s.Setup.Gs
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
		// append H2 to the left part of commitment basis for the pedersen
		// commitments for the left part
		l := append(cloneP(gv, 0, half(len(gv))), h2)
		// inner product between the right part of c3 coefficients (polynomial
		// created using the randomness) and left part of powers of z
		leftCz := innerProd(cv[half(len(cv)):], zv[:half(len(zv))])
		leftSetup := commit.Setup{Gs: l, S: s.Setup.S, G: s.g}
		li := commit.Commit(leftSetup, append(cloneS(cv, half(len(cv)), len(cv)), leftCz), w)
		lis = append(lis, li)
		// inner product between the left part of c3 coefficients and right part
		// of power z
		rightCz := innerProd(cv[:half(len(cv))], zv[half(len(zv)):])
		// commitment basis for the right pedersen commitments
		r := append(cloneP(gv, half(len(gv)), len(gv)), h2)
		rightSetup := commit.Setup{Gs: r, S: s.Setup.S, G: s.g}
		ri := commit.Commit(rightSetup, append(cloneS(cv, 0, half(len(cv))), rightCz), w)
		ris = append(ris, ri)
		// get oracle from these commits
		ei = s.g.Scalar().Pick(oracleFrom(s.g, ei, li, ri))
		// construct commitment key for next round:
		// left part of gv + ei * right part of gv
		pgv := cloneP(gv, half(len(gv)), len(gv))
		for i := 0; i < len(pgv); i++ {
			pgv[i] = pgv[i].Mul(ei, pgv[i])
		}
		gv = append(cloneP(gv, 0, half(len(gv))), pgv...)
		// and inputs for next round cv and zv
		// scale right part of cv by negative ei
		pcv := cloneS(cv, half(len(cv)), len(cv))
		scale(pcv, s.g.Scalar().Neg(ei))
		// left(ci) + ei * right(ci)
		cv = append(cloneS(cv, 0, half(len(cv))), pcv...)
		// same for z
		pzv := cloneS(zv, half(len(zv)), len(zv))
		scale(pzv, s.g.Scalar().Neg(ei))
		zv = append(cloneS(zv, 0, half(len(cv))), pzv...)
	}
	return Proof{
		U:  gv,
		C:  cv,
		C2: c2,
		W:  w,
		Ri: ris,
		Li: lis,
	}
}

type Proof struct {
	U []Point
	C []Scalar
	// Hiding commitment
	C2 Point
	// Commitment randomness
	W  Scalar
	Ri []Point
	Li []Point
}

func cloneP(a []Point, low, high int) []Point {
	b := make([]Point, high-low)
	for i := low; i < high; i++ {
		bi := i - low
		b[bi] = a[i]
	}
	return b
}

func cloneS(s []Scalar, low, high int) []Scalar {
	b := make([]Scalar, high-low)
	for i := low; i < high; i++ {
		bi := i - low
		b[bi] = s[i]
	}
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
	if len(a) != len(b) {
		panic("inner prod not possible")
	}
	sum := a[0].Clone().Zero()
	tmp := a[0].Clone()
	for i := 0; i < len(a); i++ {
		sum = sum.Add(sum, tmp.Mul(a[i], b[i]))
	}
	return sum
}

// returns an oracle for deriving a point for pedersen polynomial commitment
func oracle(g Group, degree int) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-polycommit"))
	h.Write([]byte(g.String()))
	binary.Write(h, binary.LittleEndian, degree)
	var b bytes.Buffer
	b.Write(h.Sum(nil))
	return random.New(&b)
}

func oracleFrom(g Group, infos ...kyber.Marshaling) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-polycommit-scalar"))
	h.Write([]byte(g.String()))
	for _, s := range infos {
		s.MarshalTo(h)
	}
	var b bytes.Buffer
	b.Write(h.Sum(nil))
	return random.New(&b)
}
