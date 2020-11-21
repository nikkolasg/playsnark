package pedersen

import (
	"testing"

	"github.com/drand/kyber/group/edwards25519"
	"github.com/drand/kyber/util/random"
	poly "github.com/nikkolasg/playsnark/polynomial"
	"github.com/stretchr/testify/require"
)

var g = edwards25519.NewBlakeSHA256Ed25519()

func TestPedersenPolyCommit(t *testing.T) {
	degree := 3 // d+1 must be power of two
	setup := NewSetup(g, degree)
	// polynomial to commit to
	p := poly.RandomPoly(g, degree)
	// randomness used to commit
	w := g.Scalar().Pick(random.New())
	// commitment to the polynomial
	commitment := Commit(setup, p, w)
	// evaluation point we wish to evaluate and open our polynomial
	z := g.Scalar().Pick(random.New())
	v := p.Eval(z)
	// opening of the commitment for the given evaluation point
	proof := Open(setup, p, commitment, z, w)
	// verification of the correct opening
	require.True(t, Verify(setup, commitment, z, v, proof))
}
