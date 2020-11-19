package playsnark

import (
	"testing"

	"github.com/drand/kyber/util/random"
	poly "github.com/nikkolasg/playsnark/polynomial"
	"github.com/stretchr/testify/require"
)

func TestAlgebraMatrixTranspose(t *testing.T) {
	var m Matrix = Matrix([]Vector{
		Vector([]Value{Value(1), Value(2), Value(3), Value(4)}),
		Vector([]Value{Value(5), Value(6), Value(7), Value(8)}),
		Vector([]Value{Value(9), Value(10), Value(11), Value(12)}),
	})
	var exp Matrix = Matrix([]Vector{
		Vector([]Value{Value(1), Value(5), Value(9)}),
		Vector([]Value{Value(2), Value(6), Value(10)}),
		Vector([]Value{Value(3), Value(7), Value(11)}),
		Vector([]Value{Value(4), Value(8), Value(12)}),
	})

	tm := m.Transpose()
	require.Len(t, m, 3)
	require.Len(t, tm, 4)
	require.Equal(t, exp, tm)
}

func TestAlgebraBlindEval(t *testing.T) {
	var d = 4
	var p = poly.RandomPoly(Group, 4)
	var x = NewElement().Pick(random.New())
	var px = p.Eval(x)
	var shift = NewElement().Pick(random.New())
	// g^p(x)
	var gpx = Group.Point().Mul(px, nil)
	// g^(p(x) * shift)
	var shiftGpx = Group.Point().Mul(shift, gpx)

	var blindedPoints = GeneratePowersCommit(zeroG1, x, shift, d)
	var res = p.BlindEval(Group.Point().Base(), blindedPoints)
	require.True(t, shiftGpx.Equal(res))
}
