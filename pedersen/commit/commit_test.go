package commit

import (
	"testing"

	"github.com/drand/kyber/group/edwards25519"
	"github.com/drand/kyber/util/random"
	"github.com/stretchr/testify/require"
)

var g = edwards25519.NewBlakeSHA256Ed25519()

func TestCommit(t *testing.T) {
	var l = 10
	setup := NewSetup(g, l)
	msgs := make([]Scalar, l-1)
	for i := range msgs {
		msgs[i] = g.Scalar().Pick(random.New())
	}

	r := g.Scalar().Pick(random.New())
	c1 := Commit(setup, msgs, r)
	c2 := Commit(setup, msgs, r)
	require.True(t, c1.Equal(c2))
	r2 := g.Scalar().Pick(random.New())
	c3 := Commit(setup, msgs, r2)
	require.False(t, c1.Equal(c3))
}
