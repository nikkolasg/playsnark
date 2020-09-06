package playsnark

import (
	"testing"

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
