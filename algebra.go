package playsnark

import (
	"bytes"
	"fmt"
	"sort"

	"github.com/drand/kyber"
)

type Value int

type Vector []Value

type Matrix []Vector

func NewMatrix(rows []Vector) Matrix {
	return Matrix(rows)
}

func (m Matrix) String() string {
	var b bytes.Buffer
	b.WriteString("[\n")
	for _, row := range m {
		b.WriteString("  [")
		for _, v := range row {
			b.WriteString(fmt.Sprintf("%3v", v))
		}
		b.WriteString(" ]\n")
	}
	b.WriteString("]\n")
	return b.String()
}

func (m Matrix) Column(i int) Vector {
	v := make(Vector, 0, len(m))
	for _, row := range m {
		v = append(v, row[i])
	}
	return v
}

func (m Matrix) Transpose() Matrix {
	nbRow := len(m[0])
	m2 := make(Matrix, 0, nbRow)
	for i := 0; i < nbRow; i++ {
		m2 = append(m2, m.Column(i))
	}
	return m2
}

func (m Matrix) Mul(v Vector) Vector {
	var out Vector
	for _, row := range m {
		var acc Value
		for i := range row {
			acc += row[i] * v[i]
		}
		out = append(out, acc)
	}
	return out
}

func (v Vector) Hadamard(v2 Vector) Vector {
	out := make(Vector, 0, len(v))
	for i := range v {
		out = append(out, v[i]*v2[i])
	}
	return out
}

func (v Vector) Sub(v2 Vector) Vector {
	out := make(Vector, len(v))
	for i := range v {
		out = append(out, v[i]-v2[i])
	}
	return out
}

func (v Vector) IsZero() bool {
	for _, e := range v {
		if e != 0 {
			return false
		}
	}
	return true
}

type Poly []Element
type PolyCommit []Commit

func (p Poly) Mul(p2 Poly) Poly {
	l := len(p) + len(p2) - 1
	output := make(Poly, l)
	for i := 0; i < l; i++ {
		output[i] = NewElement()
	}
	for i, v1 := range p {
		for j, v2 := range p2 {
			tmp := NewElement().Mul(v1, v2)
			output[i+j] = output[i+j].Add(output[i+j], tmp)
		}
	}
	return output
}

func (p Poly) Eval(i Element) Element {
	xi := i.Clone()
	v := zero.Clone()
	for j := len(p) - 1; j >= 0; j-- {
		v.Mul(v, xi)
		v.Add(v, p[j])
	}
	return v
}

// Div returns the quotient p / p2 and the remainder using polynomial synthetic
// division
func (p Poly) Div(p2 Poly) (q Poly, r Poly) {
	dividend := p
	divisor := p2
	out := make(Poly, len(dividend))
	for i, c := range dividend {
		out[i] = NewElement().Set(c)
	}
	for i := 0; i < len(dividend)-(len(divisor)-1); i++ {
		out[i].Div(out[i], divisor[0])
		if coef := out[i]; !coef.Equal(zero) {
			var a = NewElement()
			for j := 1; j < len(divisor); j++ {
				out[i+j].Add(out[i+j], a.Mul(a.Neg(divisor[j]), coef))
			}
		}
	}
	separator := len(out) - (len(divisor) - 1)
	return out[:separator], out[separator:]
}

// Long polynomial division
func (p *Poly) Div2(p2 Poly) (q Poly, r Poly) {
	r = p.Clone()
	for len(r) > 0 && len(r) >= len(p2) {
		num := r[len(r)-1].Clone()
		t := num.Div(num, p2[len(p2)-1])
		degreeT := len(r) - len(p2)
		tPoly := newPoly(degreeT)
		tPoly[len(tPoly)-1] = t
		q = q.Add(tPoly)
		// tPoly is n-th coefficient of p / highest of p2 (degree m)
		// (a / b) * x^(n-m)
		// so tPoly * p2 has degree n
		// tPoly * p2 = a * x^n + ... x^n-1 + ...
		// so we can remove the last coefficient of "r" since r always contains
		// the highest coefficient of p not "removed" so far
		r = r.Sub(tPoly.Mul(p2))[:len(r)-1]
	}
	r.Normalize()
	return
}

func (p Poly) Add(p2 Poly) Poly {
	max := len(p)
	if max < len(p2) {
		max = len(p2)
	}

	output := make(Poly, max)
	for i := range p {
		output[i] = NewElement().Set(p[i])
	}
	for i := range p2 {
		if output[i] == nil {
			output[i] = NewElement()
		}
		output[i] = output[i].Add(output[i], p2[i])
	}
	return output
}

func (p Poly) Sub(p2 Poly) Poly {
	max := len(p)
	if max < len(p2) {
		max = len(p2)
	}

	output := make(Poly, max)
	for i := range p {
		output[i] = NewElement().Set(p[i])
	}
	for i := range p2 {
		if output[i] == nil {
			output[i] = NewElement()
		}
		output[i] = output[i].Sub(output[i], p2[i])
	}
	return output
}
func (p Poly) Equal(p2 Poly) bool {
	if len(p) != len(p) {
		return false
	}

	for i := 0; i < len(p); i++ {
		if !p[i].Equal(p2[i]) {
			return false
		}
	}
	return true
}

func (p Poly) Clone() Poly {
	o := make(Poly, len(p))
	for i := 0; i < len(p); i++ {
		o[i] = p[i].Clone()
	}
	return o
}

func newPoly(d int) Poly {
	o := make(Poly, d+1)
	for i := 0; i <= d; i++ {
		o[i] = NewElement()
	}
	return o
}

// Normalize remove all the 0 coefficients from the highest degree downwards
// until it encounters a non zero coefficients (i.e. len(p) will give the degree
// of the coefficient)
func (p Poly) Normalize() Poly {
	maxi := len(p)
	for i := len(p) - 1; i >= 0; i-- {
		if !p[i].Equal(zero) {
			return p[:maxi]
		}
		maxi--
	}
	return p[:maxi]
}

func (p Poly) Degree() int {
	return len(p) - 1
}

type pair struct {
	I int
	V Element
}

// Interpolate takes a list of element [ y_1, ...  y_n] and returns
// a polynomial p such that (note the indices)
// p(1) = y_1, p(2) = y_2, ... p(n) = y_n
// Code largely taken from github.com/drand/kyber
func Interpolate(ys []Element) Poly {
	var pairs []pair
	for i, y := range ys {
		pairs = append(pairs, pair{I: i + 1, V: y})
	}
	x, y := xyScalar(Group, pairs)

	var accPoly = Poly([]Element{zero})
	//den := g.Scalar()
	// Notations follow the Wikipedia article on Lagrange interpolation
	// https://en.wikipedia.org/wiki/Lagrange_polynomial
	for j := range x {
		basis := lagrangeBasis(Group, j, x)
		for i := range basis {
			basis[i] = basis[i].Mul(basis[i], y[j])
		}

		if accPoly == nil {
			accPoly = basis
			continue
		}

		// add all L_j * y_j together
		accPoly = accPoly.Add(basis)
	}
	return accPoly
}

type byIndexScalar []pair

func (s byIndexScalar) Len() int           { return len(s) }
func (s byIndexScalar) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }
func (s byIndexScalar) Less(i, j int) bool { return s[i].I < s[j].I }

// xyScalar returns the list of (x_i, y_i) pairs indexed. The first map returned
// is the list of x_i and the second map is the list of y_i, both indexed in
// their respective map at index i.
func xyScalar(g kyber.Group, shares []pair) (map[int]Element, map[int]Element) {
	// we are sorting first the shares since the shares may be unrelated for
	// some applications. In this case, all participants needs to interpolate on
	// the exact same order shares.
	sorted := make([]pair, 0, len(shares))
	for _, share := range shares {
		sorted = append(sorted, share)
	}
	sort.Sort(byIndexScalar(sorted))

	x := make(map[int]kyber.Scalar)
	y := make(map[int]kyber.Scalar)
	for _, s := range sorted {
		if s.V == nil || s.I < 0 {
			continue
		}
		idx := s.I
		x[idx] = g.Scalar().SetInt64(int64(idx))
		y[idx] = s.V
	}
	return x, y
}

// lagrangeBasis returns a PriPoly containing the Lagrange coefficients for the
// i-th position. xs is a mapping between the indices and the values that the
// interpolation is using, computed with xyScalar().
func lagrangeBasis(g kyber.Group, i int, xs map[int]Element) Poly {
	var basis = Poly([]Element{one.Clone()})
	// compute lagrange basis l_j
	den := g.Scalar().One()
	var acc = g.Scalar().One()
	for m, xm := range xs {
		if i == m {
			continue
		}
		// multiply by x -i
		basis = basis.Mul(Poly([]Element{NewElement().Neg(xm), one.Clone()}))
		den.Sub(xs[i], xm) // den = xi - xm
		den.Inv(den)       // den = 1 / den
		acc.Mul(acc, den)  // acc = acc * den
	}

	// multiply all coefficients by the denominator
	for i := range basis {
		basis[i] = basis[i].Mul(basis[i], acc)
	}
	return basis
}

// blindEvaluation takes a polynomial p(x) and a list of blinded points {s}
// such that the i-th value in blindedPoint is equal to s^i, s being unknown
// from the trusted setup.
// the result is SUM( g^(s^i)^p[i] ) <=> (in addition form) SUM(p[i] * (s^i * g)
// which is equivalent to g^p(s)
// We blindly evaluate for all coefficients of p, blindedPoints can be of higher
// degree it doesn't affect the result, but it must have at least the same
// degree as p
func (p Poly) BlindEval(zero Commit, blindedPoint []Commit) Commit {
	// XXX change that to equality
	if len(p) != len(blindedPoint) {
		panic(fmt.Sprintf("mismatch of length between poly %d and blinded eval points %d", len(p), len(blindedPoint)))
	}
	var acc = zero.Clone()
	var tmp = zero.Clone()
	for i := 0; i < len(p); i++ {
		acc = acc.Add(acc, tmp.Mul(p[i], blindedPoint[i]))
	}
	return acc
}

func (p Poly) Commit(base Commit) PolyCommit {
	var pp = make([]Commit, 0, len(p))
	for _, c := range p {
		bb := base.Clone()
		pp = append(pp, bb.Mul(c, bb))
	}
	return pp
}

// GeneratePowersCommit returns { g^shift * s^i} for i=0...power included
func GeneratePowersCommit(base Commit, e Element, shift Element, power int) []Commit {
	var gi = make([]Commit, 0, power+1)
	gi = append(gi, base.Clone().Mul(shift, nil))
	var si = one.Clone()
	var tmp = NewElement()
	for i := 0; i < power; i++ {
		// s * (tmp) = s * ( s * ( .. ) )
		si = si.Mul(si, e)
		// s^i * shift
		tmp = tmp.Mul(si, shift)
		gi = append(gi, base.Clone().Mul(tmp, nil))
	}
	return gi
}
