package poly

import (
	"testing"

	"github.com/drand/kyber/group/edwards25519"
	"github.com/stretchr/testify/require"
)

var group = edwards25519.NewBlakeSHA256Ed25519()

func TestAlgebraEval(t *testing.T) {
	var p Poly = NewPolyFrom(group, []Scalar{
		group.Scalar().SetInt64(1),
		group.Scalar().SetInt64(1),
	})

	res := p.Eval(group.Scalar().SetInt64(1))
	exp := group.Scalar().SetInt64(2)
	require.True(t, res.Equal(exp))
}

func TestAlgebraInterpolate(t *testing.T) {
	var r = RandomPoly(group, 4)
	p := Interpolate(group, r.c)
	for i := 0; i < len(r.c); i++ {
		res := p.Eval(group.Scalar().SetInt64(int64(i + 1)))
		exp := r.c[i]
		require.True(t, res.Equal(exp))
	}
}

// Test the minimal polynomial division theorem in algebra
func TestAlgebraMinimal(t *testing.T) {
	// t(x) = h(x) * z(x)
	// where
	// z(x) is the minimal polynomial
	// <=> z(x) = (x-1) * (x-2)
	// h(x) = 2x
	// therefore
	// t(x) = (x-1)*(x-2)*2x = 2x^3 - 6x^2 + 4x
	var p Poly = NewPolyFrom(group, []Scalar{
		group.Scalar().Zero(),
		group.Scalar().SetInt64(4),
		group.Scalar().Neg(group.Scalar().SetInt64(6)),
		group.Scalar().SetInt64(2),
	})

	var z Poly = NewPolyFrom(group, []Scalar{group.Scalar().One()})
	for i := 1; i <= 2; i++ {
		root := NewPolyFrom(group, []Scalar{
			group.Scalar().Neg(group.Scalar().SetInt64(int64(i))),
			group.Scalar().One(),
		})
		z = z.Mul(root)
	}
	// p / z  should give h(x) without any remainder
	_, rem := p.Div2(z)
	require.True(t, len(rem.Normalize().Coeffs()) == 0)
}

func scalar(i int) Scalar {
	return group.Scalar().SetInt64(int64(i))
}

func TestAlgebraPolyMul(t *testing.T) {
	// p1(x) = 1 + 2x
	var p1 Poly = NewPolyFrom(group, []Scalar{
		scalar(1),
		scalar(2),
	})
	// p2(x) = 3 + x^2
	var p2 Poly = NewPolyFrom(group, []Scalar{
		scalar(3), scalar(0), scalar(1),
	})

	// p3(x) = (1 + 2x) * (3 + x^2)
	// 		 = 3 + 6x + x^2 + 2x^3
	var exp Poly = NewPolyFrom(group, []Scalar{
		scalar(3), scalar(6), scalar(1), scalar(2),
	})

	p3 := p1.Mul(p2)

	require.Equal(t, len(exp.Coeffs()), len(p3.Coeffs()))
	for i := 0; i < len(p3.c); i++ {
		require.Equal(t, p3.c[i], exp.c[i])
	}
}

func TestAlgebraPolyAdd(t *testing.T) {
	type TestVector struct {
		A Poly
		B Poly
	}
	for row, tv := range []TestVector{
		{
			A: RandomPoly(group, 4),
			B: RandomPoly(group, 6),
		},
		{
			A: RandomPoly(group, 6),
			B: RandomPoly(group, 4),
		},
		{
			A: RandomPoly(group, 5),
			B: RandomPoly(group, 5),
		},
	} {
		x := 10
		res := tv.A.Add(tv.B)
		// test for homomorphism
		// A(x) + B(x) = (A + B)(x)
		xf := scalar(x)
		exp := tv.A.Eval(xf)
		exp = exp.Add(exp, tv.B.Eval(xf))
		require.Equal(t, res.Eval(xf), exp, "failed at row %d", row)
	}
}

func TestAlgebraPolyDivManual(t *testing.T) {
	// p2(x) = 1 + 2x
	p2 := NewPolyFrom(group, []Scalar{
		scalar(1), scalar(2),
	})
	// q = 2x
	q := NewPolyFrom(group, []Scalar{
		scalar(0), scalar(2),
	})
	// p1 = p2 * q = (1+2x) * 2x = 2x + 4x^2
	p1 := NewPolyFrom(group, []Scalar{
		scalar(0), scalar(2), scalar(4),
	})

	expQ, expR := p1.Div(p2)
	require.True(t, q.Equal(expQ))
	require.Equal(t, len(expR.Normalize().Coeffs()), 0)
}

func TestAlgebraPolyDiv(t *testing.T) {
	p1 := RandomPoly(group, 5)
	p2 := RandomPoly(group, 4)
	// p1 = q * p2 + r
	q, r := p1.Div2(p2)

	// so p2 * q + r <=> p1
	exp := q.Mul(p2).Add(r)
	require.True(t, p1.Equal(exp))

	h := NewPolyFrom(group, []Scalar{scalar(2), scalar(2)})
	// hp1 = p1 * h
	hp1 := p1.Mul(h)
	// q = hp1 / h = p1
	exp1, rem := hp1.Div2(h)
	require.True(t, exp1.Equal(p1))
	require.True(t, len(rem.Normalize().Coeffs()) == 0)
}
