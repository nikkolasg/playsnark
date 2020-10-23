package playsnark

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func TestGroth16TrustedSetup(t *testing.T) {
	r1cs := createR1CS()
	s := createWitness(r1cs)
	qap := ToQAP(r1cs)
	//diff := qap.nbVars - qap.nbIO
	//tr := NewGroth16TrustedSetup(qap)

	h := qap.Quotient(s)
	require.Equal(t, h.Degree(), qap.z.Degree()-2)
	// z = PROD (x - Xi) for all X
	require.Equal(t, h.Degree(), qap.nbGates-2)
}

func TestGroth16Verify(t *testing.T) {
	r1cs := createR1CS()
	s := createWitness(r1cs)
	qap := ToQAP(r1cs)
	diff := qap.nbVars - qap.nbIO
	tr := NewGroth16TrustedSetup(qap)
	proof := Groth16Prove(tr, qap, s)
	require.True(t, Groth16Verify(tr, qap, proof, s[:diff]))
}

func TestGroth16ProofGen(t *testing.T) {
	r1cs := createR1CS()
	s := createWitness(r1cs)
	qap := ToQAP(r1cs)
	diff := qap.nbVars - qap.nbIO
	tr := NewGroth16TrustedSetup(qap)
	proof := Groth16Prove(tr, qap, s)
	// compute the plain value of A and then put it in the exponent and verify
	// correctness
	var res = NewElement().Zero()
	var x = tr.tw.X
	for i := 0; i < qap.nbVars; i++ {
		uix := qap.left[i].Eval(x)
		res = res.Add(res, uix.Mul(uix, s[i].ToFieldElement()))
	}
	res = res.Add(res, NewElement().Mul(proof.tp.R, tr.tw.Delta))
	res = res.Add(res, tr.tw.Alpha)
	resC := NewG1().Mul(res, nil)
	require.True(t, resC.Equal(proof.A))

	// same for B
	res = NewElement().Zero()
	for i := 0; i < qap.nbVars; i++ {
		vix := qap.right[i].Eval(x)
		res = res.Add(res, vix.Mul(vix, s[i].ToFieldElement()))
	}
	res = res.Add(res, NewElement().Mul(proof.tp.S, tr.tw.Delta))
	res = res.Add(res, tr.tw.Beta)
	resB := res
	resC = NewG2().Mul(res, nil)
	require.True(t, resC.Equal(proof.B))

	// same for C even though a bee more complex
	res = NewElement().Zero()
	for i := diff; i < qap.nbVars; i++ {
		// u_i(x)
		uix := qap.left[i].Eval(x)
		vix := qap.right[i].Eval(x)
		wix := qap.out[i].Eval(x)
		// beta * u_i(x)
		buix := NewElement().Mul(tr.tw.Beta, uix)
		avix := NewElement().Mul(tr.tw.Alpha, vix)
		// sum of beta * u_i(x) + alpha * v_i(x) + w_i(x)
		sum := wix.Add(wix, buix.Add(buix, avix))
		// ai / delta
		ad := NewElement().Div(s[i].ToFieldElement(), tr.tw.Delta)
		tot := sum.Mul(ad, sum)
		res = res.Add(res, tot)
	}

	// h(x) * t(x) part
	h := qap.Quotient(s)
	require.True(t, h.Degree() == len(tr.XiT)-1)
	ht := NewElement().Mul(h.Eval(x), qap.z.Eval(x))
	htd := ht.Div(ht, tr.tw.Delta)
	// res + h(x)*t(x) / delta
	res = res.Add(res, htd)
	resC = NewG1().Mul(res, nil)

	// As
	As := NewG1().Mul(proof.tp.S, proof.A)
	resC = resC.Add(resC, As)
	// Br - we already have DLog of B computed before
	B1 := NewG1().Mul(resB, nil)
	B1r := B1.Mul(proof.tp.R, B1)
	resC = resC.Add(resC, B1r)

	// r*s*delta
	rs := NewElement().Mul(proof.tp.R, proof.tp.S)
	rsd := rs.Mul(rs, tr.tw.Delta)
	rsd = rsd.Neg(rsd)
	rsdC := NewG1().Mul(rsd, nil)
	resC = resC.Add(resC, rsdC)

	require.True(t, resC.Equal(proof.C))
}
