package playsnark

import (
	"testing"

	"github.com/drand/kyber/util/random"
	"github.com/stretchr/testify/require"
)

func TestPolyMul(t *testing.T) {
	// p1(x) = 1 + 2x
	var p1 Poly = []Element{
		Value(1).ToFieldElement(),
		Value(2).ToFieldElement(),
	}
	// p2(x) = 3 + x^2
	var p2 Poly = []Element{
		Value(3).ToFieldElement(),
		NewElement(),
		Value(1).ToFieldElement(),
	}

	// p3(x) = (1 + 2x) * (3 + x^2)
	// 		 = 3 + 6x + x^2 + 2x^3
	var exp Poly = []Element{
		Value(3).ToFieldElement(),
		Value(6).ToFieldElement(),
		Value(1).ToFieldElement(),
		Value(2).ToFieldElement(),
	}

	p3 := p1.Mul(p2)

	require.Equal(t, exp.Degree(), p3.Degree())
	for i := 0; i < p3.Degree(); i++ {
		require.Equal(t, p3[i], exp[i])
	}
}

func TestPolyAdd(t *testing.T) {
	type TestVector struct {
		A Poly
		B Poly
	}
	for row, tv := range []TestVector{
		{
			A: randomPoly(4),
			B: randomPoly(6),
		},
		{
			A: randomPoly(6),
			B: randomPoly(4),
		},
		{
			A: randomPoly(5),
			B: randomPoly(5),
		},
	} {
		x := 10
		res := tv.A.Add(tv.B)
		// test for homomorphism
		// A(x) + B(x) = (A + B)(x)
		exp := tv.A.Eval(x)
		exp = exp.Add(exp, tv.B.Eval(x))
		require.Equal(t, res.Eval(x), exp, "failed at row %d", row)
	}
}

func TestPolyDivManual(t *testing.T) {
	// p2(x) = 1 + 2x
	p2 := Poly([]Element{
		Value(1).ToFieldElement(),
		Value(2).ToFieldElement(),
	})
	// q = 2x
	q := Poly([]Element{
		NewElement(),
		Value(2).ToFieldElement(),
	})
	// p1 = p2 * q = (1+2x) * 2x = 2x + 4x^2
	p1 := Poly([]Element{
		NewElement(),
		Value(2).ToFieldElement(),
		Value(4).ToFieldElement(),
	})

	expQ, expR := p1.Div(p2)
	require.True(t, q.Equal(expQ))
	require.Equal(t, expR.Degree(), -1)
}

func TestPolyDiv(t *testing.T) {
	p1 := randomPoly(5)
	p2 := randomPoly(4)
	// p1 = q * p2 + r
	q, r := p1.Div(p2)

	// so p2 * q + r <=> p1
	exp := q.Mul(p2).Add(r)
	require.True(t, p1.Equal(exp))
}

func randomPoly(d int) Poly {
	var poly = make(Poly, 0, d+1)
	for i := 0; i <= d; i++ {
		poly = append(poly, NewElement().Pick(random.New()))
	}
	return poly
}
