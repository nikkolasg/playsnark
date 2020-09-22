package playsnark

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestQAPManual(t *testing.T) {
	r1cs := createR1CS()
	qap := ToQAP(r1cs)
	s := createWitness()
	require.Len(t, r1cs.left, 4)
	require.Len(t, r1cs.left.Transpose(), len(s))
	require.Len(t, qap.left, len(s))
	// so for example, all the assignements on the lefts inputs look like this
	// in R1CS:
	// [0, 1, 0, 0, 0, 0]
	// [0, 0, 0, 1, 0, 0]
	// [0, 1, 0, 0, 1, 0]
	// [5, 0, 0, 0, 0, 1]
	// so for example the second column of the transposed matrix represents all
	// the usage of the variable "x" as a left wire in the circuit -> [1,0,1,0]
	// Then we need to interpolate this as a polynomial such that:
	// p(1) = 1, p(2) = 0, p(3) = 1 and p(4) = 0
	// Note we have shifted by one the entries
	xPoly := qap.left[1]
	require.True(t, xPoly.Eval(Value(1).ToFieldElement()).
		Equal(Value(1).ToFieldElement()))
	require.True(t, xPoly.Eval(Value(2).ToFieldElement()).
		Equal(Value(0).ToFieldElement()))
	require.True(t, xPoly.Eval(Value(3).ToFieldElement()).
		Equal(Value(1).ToFieldElement()))
	require.True(t, xPoly.Eval(Value(4).ToFieldElement()).
		Equal(Value(0).ToFieldElement()))

	// Check the the assignement for the first gate
	// We need to evaluate all polynomials for all gates, so x = 1, x=2...,
	// multiply the result by the solution value of the variable and check the
	// usual equation
	for gate := 1; gate <= qap.nbGates; gate++ {
		gateF := Value(gate).ToFieldElement()
		left := zero.Clone()
		right := zero.Clone()
		out := zero.Clone()
		for varIndex, val := range s {
			vfe := val.ToFieldElement()
			prod := NewElement().Mul(vfe, qap.left[varIndex].Eval(gateF))
			left = left.Add(left, prod)

			prod = NewElement().Mul(vfe, qap.right[varIndex].Eval(gateF))
			right = right.Add(right, prod)

			prod = NewElement().Mul(vfe, qap.out[varIndex].Eval(gateF))
			out = out.Add(out, prod)
		}
		// check the equation left * right - out = 0
		res := out.Sub(out, NewElement().Mul(left, right))
		require.True(t, res.Equal(zero), "gate %d", gate)
	}
}

func TestQAPValidity(t *testing.T) {
	s := createWitness()
	r1cs := createR1CS()
	qap := ToQAP(r1cs)
	require.True(t, qap.IsValid(s))
	fmt.Println(qap.nbGates)
	fmt.Println(qap.right)
}

func TestQAPInterpolate(t *testing.T) {
	_ = createWitness()
	r1cs := createR1CS()
	fmt.Println(r1cs.out)
	polys := qapInterpolate(r1cs.out)
	fmt.Println(polys)

}
