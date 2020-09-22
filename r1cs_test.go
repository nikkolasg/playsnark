package playsnark

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestR1CSEquation(t *testing.T) {
	r1cs := createR1CS()

	// r1cs equation means that given a vector of solution s and matrices Left
	// Right and Out, the following equation is true
	// Left * s^T  x Right * s^T - Out * s^T = 0
	// where x is the hadamard product

	s := createWitness()
	l := r1cs.left.Mul(s)
	r := r1cs.right.Mul(s)
	o := r1cs.out.Mul(s)
	left := l.Hadamard(r).Sub(o)
	require.True(t, left.IsZero())
	require.Len(t, r1cs.left[0], len(s))
	fmt.Println(r1cs.right)
}
