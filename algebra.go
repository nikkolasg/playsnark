package playsnark

import (
	"github.com/drand/kyber/share"
)

type Value int

type Vector []Value

type Matrix []Vector

func NewMatrix(rows []Vector) Matrix {
	return Matrix(rows)
}

func (m Matrix) Column(i int) Vector {
	v := make(Vector, 0, len(m))
	for _, row := range m {
		v = append(v, row[i])
	}
	return v
}

func (m Matrix) Transpose() Matrix {
	m2 := make(Matrix, 0, len(m[0]))
	for i := range m {
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

func (p Poly) Degree() int {
	for d := len(p) - 1; d >= 0; d-- {
		if !p[d].Equal(zero) {
			return d
		}
	}
	return -1
}

func (p Poly) Mul(p2 Poly) Poly {
	l := p.Degree() + p2.Degree() + 1
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

func (p Poly) Eval(i int) Element {
	ps := share.CoefficientsToPriPoly(Group, p)
	sh := ps.Eval(i)
	return sh.V
}

var zero = NewElement()

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
func (p Poly) Div2(p2 Poly) (q Poly, r Poly) {
	r = p.Clone()
	for len(r) > 0 && r.Degree() >= p2.Degree() {
		num := r[len(r)-1].Clone()
		t := num.Div(num, p2[len(p2)-1])
		degreeT := r.Degree() - p2.Degree()
		tPoly := newPoly(degreeT)
		tPoly[len(tPoly)-1] = t
		q = q.Add(tPoly)
		r = r.Sub(tPoly.Mul(p2))
	}
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
	if p.Degree() != p2.Degree() {
		return false
	}

	for i := 0; i < p.Degree(); i++ {
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
