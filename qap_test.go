package playsnark

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func TestQAPSnarkManual(t *testing.T) {
	r1cs := createR1CS()
	qap := ToQAP(r1cs)
	// so for example, all the assignements on the lefts inputs look like this
	// in R1CS:
	// [0, 1, 0, 0, 0, 0]
	// [0, 0, 0, 1, 0, 0]
	// [0, 1, 0, 0, 1, 0]
	// [5, 0, 0, 0, 0, 1]
	// so for example the second row of the transposed matrix represents all the
	// usage of the variable "x" as a left wire in the circuit -> [1,0,1,0]
	// Then we need to interpolate this as a polynomial such that:
	// p(0) = 1, p(1) = 0, p(2) = 1 and p(3) = 0
	xPoly := qap.left[1]
	require.True(t, xPoly.Eval(0).Equal(Value(1).ToFieldElement()))
	require.True(t, xPoly.Eval(1).Equal(Value(0).ToFieldElement()))
	require.True(t, xPoly.Eval(2).Equal(Value(1).ToFieldElement()))
	require.True(t, xPoly.Eval(3).Equal(Value(0).ToFieldElement()))
}
