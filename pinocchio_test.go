package playsnark

import (
	"testing"

	"github.com/drand/kyber/util/random"
	"github.com/stretchr/testify/require"
)

func TestPinocchioCombine(t *testing.T) {
	var d = 4
	var p = randomPoly(4)
	var x = NewElement().Pick(random.New())
	var px = p.Eval(x)
	var gpx = NewG1().Mul(px, nil)

	var blindedPoints = generatePowersCommit(zeroG1, x, one, d)
	var res = p.BlindEval(zeroG1, blindedPoints)
	require.True(t, gpx.Equal(res))
}

func TestPinocchioValidProof(t *testing.T) {
	s := createWitness()
	r1cs := createR1CS()
	qap := ToQAP(r1cs)
	diff := qap.nbVars - qap.nbIO
	setup := NewTrustedSetup(qap)
	proof := GenProof(setup, qap, s)
	require.True(t, VerifyProof(setup.VK, qap, proof, s[:diff]))
}
