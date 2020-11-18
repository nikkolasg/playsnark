package bproof

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"
)

// IsSatisfied looks if
// 1. All multiplication constraints are satisfied
//		--> the vector a_i * b_i = c_i for all i
// 2. All linear combinations are satisfied
// SUM a_i * w_l_K_i + SUM b_i * w_r_i + SUM c_i * w_o_i = w_k
func (c *ProverCircuit) IsSatisfied(t *testing.T) bool {
	n := c.Circuit.nbVars
	// Check 1
	for i := 0; i < n; i++ {
		r := NewScalar().Mul(c.a[i], c.b[i])
		require.True(t, r.Equal(c.c[i]))
	}

	// Check 2
	c.Flatten()
	q := len(c.constraints)
	// for each linear constraint
	for k := 0; k < q; k++ {
		sum := zero()
		// iterate over all variables and associated coefficient
		wl := c.Circuit.wl[k]
		wr := c.Circuit.wr[k]
		wo := c.Circuit.wo[k]
		wk := c.Circuit.wk[k]
		for i := 0; i < n; i++ {
			// w_l_i * a_i
			sum.Add(sum, NewScalar().Mul(wl[i], c.a[i]))
			// w_r_i * b_i
			sum.Add(sum, NewScalar().Mul(wr[i], c.b[i]))
			// w_o_i * c_i
			sum.Add(sum, NewScalar().Mul(wo[i], c.c[i]))
		}
		require.True(t, sum.Equal(wk))
	}
	return true
}

func TestCircuitProverValid(t *testing.T) {
	prover := GenerateCircuit()
	prover.Flatten()
	fmt.Println(prover)
	prover.IsSatisfied(t)
}
