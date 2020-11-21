package poly

import (
	"fmt"
	"sort"

	"github.com/drand/kyber"
	"go.dedis.ch/kyber/v3/util/random"
)

// Group has points on it and can create scalar from the scalar fields
type Group = kyber.Group

// Scalar of the field of the curve
type Scalar = kyber.Scalar

// Point in the group (in our case it's elliptic curve so it's a point)
type Point = kyber.Point

type Poly struct {
	c []Scalar
	g Group
}

type PolyCommit []Point

func emptyPoly(g Group, d int) Poly {
	o := make([]Scalar, d+1)
	for i := 0; i <= d; i++ {
		o[i] = g.Scalar()
	}
	return Poly{c: o, g: g}
}

func NewZeroPoly(g Group, degree ...int) Poly {
	if len(degree) > 0 {
		return emptyPoly(g, degree[0])
	}
	return Poly{c: []Scalar{}, g: g}
}

func NewPolyFrom(g Group, coeffs []Scalar) Poly {
	return Poly{c: coeffs, g: g}
}

func (p Poly) Set(pos int, coeff Scalar) {
	p.c[pos] = coeff
}

func (p Poly) Mul(p2 Poly) Poly {
	l := len(p.c) + len(p2.c) - 1
	output := make([]Scalar, l)
	for i := 0; i < l; i++ {
		output[i] = p.g.Scalar().Zero()
	}
	for i, v1 := range p.c {
		for j, v2 := range p2.c {
			tmp := p.g.Scalar().Mul(v1, v2)
			output[i+j] = output[i+j].Add(output[i+j], tmp)
		}
	}
	return Poly{c: output, g: p.g}
}

func (p Poly) Scale(a Scalar) {
	for i := range p.c {
		p.c[i] = p.c[i].Mul(p.c[i], a)
	}
}

func (p Poly) Eval(i Scalar) Scalar {
	xi := i.Clone()
	v := p.g.Scalar().Zero()
	for j := len(p.c) - 1; j >= 0; j-- {
		v.Mul(v, xi)
		v.Add(v, p.c[j])
	}
	return v
}

// Div returns the quotient p / p2 and the remainder using polynomial synthetic
// division
func (p Poly) Div(p2 Poly) (q Poly, r Poly) {
	dividend := p.c
	divisor := p2.c
	out := make([]Scalar, len(dividend))
	for i, c := range dividend {
		out[i] = p.g.Scalar().Set(c)
	}
	for i := 0; i < len(dividend)-(len(divisor)-1); i++ {
		out[i].Div(out[i], divisor[0])
		if coef := out[i]; !coef.Equal(p.g.Scalar().Zero()) {
			var a = p.g.Scalar()
			for j := 1; j < len(divisor); j++ {
				out[i+j].Add(out[i+j], a.Mul(a.Neg(divisor[j]), coef))
			}
		}
	}
	separator := len(out) - (len(divisor) - 1)
	return Poly{c: out[:separator], g: p.g}, Poly{c: out[separator:], g: p.g}
}

// Long polynomial division
func (p *Poly) Div2(p2 Poly) (q Poly, r Poly) {
	r = p.Clone()
	q = NewPolyFrom(p.g, nil)
	for len(r.c) > 0 && len(r.c) >= len(p2.c) {
		num := r.c[len(r.c)-1].Clone()
		t := num.Div(num, p2.c[len(p2.c)-1])
		degreeT := len(r.c) - len(p2.c)
		tPoly := emptyPoly(p.g, degreeT)
		tPoly.c[len(tPoly.c)-1] = t
		q = q.Add(tPoly)
		// tPoly is n-th coefficient of p / highest of p2 (degree m)
		// (a / b) * x^(n-m)
		// so tPoly * p2 has degree n
		// tPoly * p2 = a * x^n + ... x^n-1 + ...
		// so we can remove the last coefficient of "r" since r always contains
		// the highest coefficient of p not "removed" so far
		r = NewPolyFrom(p.g, r.Sub(tPoly.Mul(p2)).c[:len(r.c)-1])
	}
	r.Normalize()
	return
}

func (p Poly) Add(p2 Poly) Poly {
	max := len(p.c)
	if max < len(p2.c) {
		max = len(p2.c)
	}

	output := make([]Scalar, max)
	for i := range p.c {
		output[i] = p.g.Scalar().Set(p.c[i])
	}
	for i := range p2.c {
		if output[i] == nil {
			output[i] = p.g.Scalar()
		}
		output[i] = output[i].Add(output[i], p2.c[i])
	}
	return Poly{c: output, g: p.g}
}

func (p Poly) Sub(p2 Poly) Poly {
	max := len(p.c)
	if max < len(p2.c) {
		max = len(p2.c)
	}

	output := make([]Scalar, max)
	for i := range p.c {
		output[i] = p.g.Scalar().Set(p.c[i])
	}
	for i := range p2.c {
		if output[i] == nil {
			output[i] = p.g.Scalar()
		}
		output[i] = output[i].Sub(output[i], p2.c[i])
	}
	return NewPolyFrom(p.g, output)
}
func (p Poly) Equal(p2 Poly) bool {
	if len(p.c) != len(p.c) {
		return false
	}

	for i := 0; i < len(p.c); i++ {
		if !p.c[i].Equal(p2.c[i]) {
			return false
		}
	}
	return true
}

func (p Poly) Clone() Poly {
	o := make([]Scalar, len(p.c))
	for i := 0; i < len(p.c); i++ {
		o[i] = p.c[i].Clone()
	}
	return Poly{c: o, g: p.g}
}

func (p Poly) Coeffs() []Scalar {
	return p.c
}

// Normalize remove all the 0 coefficients from the highest degree downwards
// until it encounters a non zero coefficients (i.e. len(p) will give the degree
// of the coefficient)
func (p Poly) Normalize() Poly {
	maxi := len(p.c)
	for i := len(p.c) - 1; i >= 0; i-- {
		if !p.c[i].Equal(p.g.Scalar().Zero()) {
			return NewPolyFrom(p.g, p.c[:maxi])
		}
		maxi--
	}
	return NewPolyFrom(p.g, p.c[:maxi])
}

func (p Poly) Degree() int {
	return len(p.c) - 1
}

type pair struct {
	I int
	V Scalar
}

// Interpolate takes a list of element [ y_1, ...  y_n] and returns
// a polynomial p such that (note the indices)
// p(1) = y_1, p(2) = y_2, ... p(n) = y_n
// Code largely taken from github.com/drand/kyber
func Interpolate(g Group, ys []Scalar) Poly {
	var pairs []pair
	for i, y := range ys {
		pairs = append(pairs, pair{I: i + 1, V: y})
	}
	x, y := xyScalar(g, pairs)

	var accPoly Poly
	var setPoly bool
	//den := g.Scalar()
	// Notations follow the Wikipedia article on Lagrange interpolation
	// https://en.wikipedia.org/wiki/Lagrange_polynomial
	for j := range x {
		basis := lagrangeBasis(g, j, x)
		for i := range basis.c {
			basis.c[i] = basis.c[i].Mul(basis.c[i], y[j])
		}

		if !setPoly {
			accPoly = basis
			setPoly = true
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
func xyScalar(g Group, shares []pair) (map[int]Scalar, map[int]Scalar) {
	// we are sorting first the shares since the shares may be unrelated for
	// some applications. In this case, all participants needs to interpolate on
	// the exact same order shares.
	sorted := make([]pair, 0, len(shares))
	for _, share := range shares {
		sorted = append(sorted, share)
	}
	sort.Sort(byIndexScalar(sorted))

	x := make(map[int]Scalar)
	y := make(map[int]Scalar)
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
func lagrangeBasis(g Group, i int, xs map[int]Scalar) Poly {
	var basis = NewPolyFrom(g, []Scalar{g.Scalar().One()})
	// compute lagrange basis l_j
	den := g.Scalar().One()
	var acc = g.Scalar().One()
	for m, xm := range xs {
		if i == m {
			continue
		}
		// multiply by x -i
		basis = basis.Mul(NewPolyFrom(g, []Scalar{g.Scalar().Neg(xm), g.Scalar().One()}))
		den.Sub(xs[i], xm) // den = xi - xm
		den.Inv(den)       // den = 1 / den
		acc.Mul(acc, den)  // acc = acc * den
	}

	// multiply all coefficients by the denominator
	for i := range basis.c {
		basis.c[i] = basis.c[i].Mul(basis.c[i], acc)
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
func (p Poly) BlindEval(zero Point, blindedPoint []Point) Point {
	if len(p.c) != len(blindedPoint) {
		panic(fmt.Sprintf("mismatch of length between poly %d and blinded eval points %d", len(p.c), len(blindedPoint)))
	}
	var acc = zero.Clone().Null()
	var tmp = zero.Clone().Null()
	for i := 0; i < len(p.c); i++ {
		acc = acc.Add(acc, tmp.Mul(p.c[i], blindedPoint[i]))
	}
	return acc
}

func (p Poly) Commit(base Point) PolyCommit {
	var pp = make([]Point, 0, len(p.c))
	for _, c := range p.c {
		bb := base.Clone()
		pp = append(pp, bb.Mul(c, bb))
	}
	return pp
}

func RandomPoly(g Group, d int) Poly {
	var poly = make([]Scalar, 0, d+1)
	for i := 0; i <= d; i++ {
		poly = append(poly, g.Scalar().Pick(random.New()))
	}
	return NewPolyFrom(g, poly)
}
