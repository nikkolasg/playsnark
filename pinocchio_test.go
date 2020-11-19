package playsnark

import (
	"fmt"
	"testing"

	"github.com/drand/kyber/util/random"
	poly "github.com/nikkolasg/playsnark/polynomial"
	"github.com/stretchr/testify/require"
)

func TestPinocchioCombine(t *testing.T) {
	var d = 4
	var p = poly.RandomPoly(Group, 4)
	var x = NewElement().Pick(random.New())
	var px = p.Eval(x)
	var gpx = NewG1().Mul(px, nil)

	var blindedPoints = GeneratePowersCommit(zeroG1, x, one, d)
	var res = p.BlindEval(zeroG1, blindedPoints)
	require.True(t, gpx.Equal(res))
}

func TestPinocchioProofValidDivision(t *testing.T) {
	r1cs := createR1CS()
	s := createWitness(r1cs)
	qap := ToQAP(r1cs)
	diff := qap.nbVars - qap.nbIO
	setup := NewPHGR13TrustedSetup(qap)
	proof := PHGR13Prove(setup.EK, qap, s)

	// test GHS
	// compute h(x) then evaluate it blindly at point s
	left, right, out := qap.computeAggregatePoly(s)
	// p(x) = t(x) * h(x)
	px := left.Mul(right).Sub(out)
	// h(x) = p(x) / t(x)
	hx, rem := px.Div2(qap.z)
	if len(rem.Normalize().Coeffs()) > 0 {
		panic("apocalypse")
	}
	hs := hx.Eval(setup.t.s)
	// h(s) * G
	HS := NewG1().Mul(hs, nil)
	require.True(t, proof.hs.Equal(HS))
	// test correct computation of verification key.yts = t(s) * G2
	// use the pairing operation to check the exponents namely
	// e(h(s) * G1, t(s) * y_y * G2) ==  e(h(s) * y_y * t(s)*G1,G2)
	leftP := Pair(HS, setup.VK.yts)
	ryts := NewElement().Mul(qap.z.Eval(setup.t.s), setup.t.ry)
	// h(s) * r_y * t(s)
	r1 := NewG1().Mul(hs, nil)
	r2 := NewG2().Mul(ryts, nil)
	// e(G1,G2)^(h(s) * r_y * t(s))
	rightP := Pair(r1, r2)
	require.True(t, leftP.Equal(rightP))
	// test gvmids
	// compute g_v^(SUM v_k(s) * sol[k]) for k being NON IO
	var vks = NewElement()
	for i, vk := range qap.left[diff:] {
		vks.Add(vks, NewElement().Mul(vk.Eval(setup.t.s), s[diff+i].ToFieldElement()))
	}
	var gvks = setup.t.gv.Clone().Mul(vks, setup.t.gv)
	require.True(t, proof.vss.Equal(gvks))

	// test gwmids
	// compute g_v^(SUM v_k(s) * sol[k]) for k being NON IO
	var wks = NewElement()
	for i, wk := range qap.right[diff:] {
		wks.Add(wks, NewElement().Mul(wk.Eval(setup.t.s), s[diff+i].ToFieldElement()))
	}
	var gwks = setup.t.gw.Clone().Mul(wks, setup.t.gw)
	require.True(t, proof.wss.Equal(gwks))

	// test gymids
	var yks = NewElement()
	for i, yk := range qap.out[diff:] {
		yks.Add(yks, NewElement().Mul(yk.Eval(setup.t.s), s[diff+i].ToFieldElement()))
	}
	var gyks = setup.t.gy.Clone().Mul(yks, setup.t.gy)
	require.True(t, proof.yss.Equal(gyks))

	// test if verifier computes g_v^(SUM v_k(s) * c_k) for all k IO related
	gvkio := computeCommitIOSolution(zeroG1, setup.VK.vs[:diff], s[:diff])
	// compute it manually first by addng all the elements and then committing
	var vkio = NewElement()
	for i, vk := range qap.left[:diff] {
		vkio.Add(vkio, NewElement().Mul(vk.Eval(setup.t.s), s[i].ToFieldElement()))
	}
	var gvkio2 = NewG1().Mul(vkio, setup.t.gv)
	require.True(t, gvkio.Equal(gvkio2))

	// test the same for gw
	gwkio := computeCommitIOSolution(zeroG2, setup.VK.ws[:diff], s[:diff])
	var wkio = NewElement()
	for i, wk := range qap.right[:diff] {
		wkio.Add(wkio, NewElement().Mul(wk.Eval(setup.t.s), s[i].ToFieldElement()))
	}
	var gwkio2 = NewG2().Mul(wkio, setup.t.gw)
	require.True(t, gwkio.Equal(gwkio2))

	// test the same for gy
	gykio := computeCommitIOSolution(zeroG1, setup.VK.ys[:diff], s[:diff])
	var ykio = NewElement()
	for i, yk := range qap.out[:diff] {
		ykio.Add(ykio, NewElement().Mul(yk.Eval(setup.t.s), s[i].ToFieldElement()))
	}
	var gykio2 = NewG1().Mul(ykio, setup.t.gy)
	require.True(t, gykio.Equal(gykio2))

	// test if the addition of the IO and non-IO (proof part) are equal when
	// computing manually v(s) * G_v
	// here we add the io part with the non IO part
	gvs := NewG1().Add(gvkio, proof.vss)
	// here we compute all of the evaluation in the field and then commit
	var vs = NewElement()
	for i, vk := range qap.left {
		vs.Add(vs, NewElement().Mul(vk.Eval(setup.t.s), s[i].ToFieldElement()))
	}
	var gvs2 = NewG1().Mul(vs, setup.t.gv)
	require.True(t, gvs.Equal(gvs2))

	// same for gw
	gws := NewG2().Add(gwkio, proof.wss)
	var ws = NewElement()
	for i, wk := range qap.right {
		ws.Add(ws, NewElement().Mul(wk.Eval(setup.t.s), s[i].ToFieldElement()))
	}
	var gws2 = NewG2().Mul(ws, setup.t.gw)
	require.True(t, gws.Equal(gws2))

	gys := NewG1().Add(gykio, proof.yss)
	// here we compute all of the evaluation in the field and then commit
	var ys = NewElement()
	for i, yk := range qap.out {
		ys.Add(ys, NewElement().Mul(yk.Eval(setup.t.s), s[i].ToFieldElement()))
	}
	var gys2 = NewG1().Mul(ys, setup.t.gy)
	require.True(t, gys.Equal(gys2))

	// try to verify the equation "in the clear" first (and not blindly as the
	// proof is doing).
	// we want to prove that
	// r_v * v(s) * r_w * w(s) =  r_y * h(s) * t(s) + r_y * y(s)
	// r_y * v(s) * w(s) = r_y * (h(s) * t(s) + y(s))
	// v(s) * w(s) - y(s) = h(s) * t(s) = p(s) which is the QAP equation
	// LEFT
	rvvs := NewElement().Mul(setup.t.rv, left.Eval(setup.t.s))
	rwws := NewElement().Mul(setup.t.rw, right.Eval(setup.t.s))
	leftC := NewElement().Mul(rvvs, rwws)
	// RIGHT
	hsts := NewElement().Mul(hx.Eval(setup.t.s), qap.z.Eval(setup.t.s))
	ryhsts := NewElement().Mul(hsts, setup.t.ry)
	rys := NewElement().Mul(out.Eval(setup.t.s), setup.t.ry)
	rightC := NewElement().Add(ryhsts, rys)
	require.True(t, rightC.Equal(leftC))

	// check the same equation but blindly
	// e(G1,G2)^(rv * v(s) * rw * w(s))
	leftS := Pair(gvs, gws)
	leftExp1 := NewG1().Mul(rvvs, nil)
	leftExp2 := NewG2().Mul(rwws, nil)
	leftExp := Pair(leftExp1, leftExp2)
	require.True(t, leftS.Equal(leftExp))
	// put all elements on G1 and looks if it succeeds
	require.True(t, Pair(NewG1().Mul(leftC, nil), NewG2()).Equal(leftS))

	// e(G1,G2)^(h(s) * r_y * t(s)) = e(G1,G2)^(r_ys * p(s))
	right1 := Pair(proof.hs, setup.VK.yts)
	// put all elements on G1 and looks if it succeeds
	rightExpLeft := Pair(NewG1().Mul(ryhsts, nil), NewG2())
	require.True(t, right1.Equal(rightExpLeft))

	// e(G1,G2)^(ry * y(s))
	right2 := Pair(gys, NewG2().Base())
	// put all elements on G1 and looks if it succeeds
	rightExpRight := Pair(NewG1().Mul(rys, nil), NewG2().Base())
	require.True(t, right2.Equal(rightExpRight))

	// e(G1,G2)^(ry * [y(s) + p(s)])
	rightS := right1.Clone().Add(right1, right2)
	// add the two components of the right side
	// e(g1,g2)^[(ry * t(s) * h(s)] * e(g1,g2)^(ry * y(s))
	// <=> e(g1,g2)^[ry*t(s)*h(s) + ry*y(s)]
	// <=> e(g1,g2)^[ry * [t(s) * h(s) + y(s)]]
	// <=> e(g1,g2)^[ry * [p(s) + y(s)]]
	rightExp := rightExpLeft.Add(rightExpLeft, rightExpRight)
	require.True(t, rightS.Equal(rightExp))

	require.True(t, rightExp.Equal(leftExp))
	// r_v * r_w * v(s) * w(s) == r_y * (p(s) +  y(s))
	// r_y * v(s) * w(s)       == r_y * (p(s) + y(s))
	// v(s) * w(s) - y(s)      == p(s) == h(s) * t(s)
	// which is the QAP equation
	require.True(t, leftS.Equal(rightS))
	require.True(t, PHGR13Verify(setup.VK, qap, proof, s[:diff]))
}

func TestPinocchioInvalidProof(t *testing.T) {
	r1cs := createR1CS()
	s := createWitness(r1cs)
	qap := ToQAP(r1cs)
	diff := qap.nbVars - qap.nbIO
	setup := NewPHGR13TrustedSetup(qap)
	proof := PHGR13Prove(setup.EK, qap, s)
	fmt.Println(proof.String())

	// left is e(g^(a_v*v(s) + a_w*w(s) + a_y *y(s)) * beta,g^gamma)
	left := Pair(proof.gz, setup.VK.gamma)
	var vkio = NewElement()
	for i, wk := range qap.left[diff:] {
		vkio.Add(vkio, NewElement().Mul(wk.Eval(setup.t.s), s[i+diff].ToFieldElement()))
	}
	var avvs = NewG1().Mul(NewElement().Mul(vkio, setup.t.rv), nil)

	var wkio = NewElement()
	for i, wk := range qap.right[diff:] {
		wkio.Add(wkio, NewElement().Mul(wk.Eval(setup.t.s), s[i+diff].ToFieldElement()))
	}
	var awws = NewG1().Mul(NewElement().Mul(wkio, setup.t.rw), nil)
	var ykio = NewElement()
	for i, wk := range qap.out[diff:] {
		ykio.Add(ykio, NewElement().Mul(wk.Eval(setup.t.s), s[i+diff].ToFieldElement()))
	}
	var ayys = NewG1().Mul(NewElement().Mul(ykio, setup.t.ry), nil)
	ball := NewG1().Add(avvs, NewG1().Add(awws, ayys))
	ball = ball.Mul(setup.t.beta, ball)
	require.True(t, ball.Equal(proof.gz))

	expLeft := Pair(ball, setup.VK.gamma)
	require.True(t, left.Equal(expLeft))

	// t1 := e(g^(a_v*v(s)) * g^(a_y*y(s)), g^beta*gamma)
	// t2 := e(g^beta*gamma, g^(a_w*w(s))
	// right = t1*t2 = left !
	lt1 := NewG1().Add(proof.vss, proof.yss)
	t1 := Pair(lt1, setup.VK.bgamma2)
	t2 := Pair(setup.VK.bgamma, proof.wss)
	right := zeroGT.Clone().Add(t1, t2)
	require.True(t, right.Equal(left))

	require.True(t, PHGR13Verify(setup.VK, qap, proof, s[:diff]))

	p2 := proof
	p2.yss = NewG1().Pick(random.New())
	require.False(t, PHGR13Verify(setup.VK, qap, p2, s[:diff]))

	p2 = proof
	p2.vss = NewG1().Pick(random.New())
	require.False(t, PHGR13Verify(setup.VK, qap, p2, s[:diff]))

	p2 = proof
	p2.wss = NewG2().Pick(random.New())
	require.False(t, PHGR13Verify(setup.VK, qap, p2, s[:diff]))

	p2 = proof
	p2.hs = NewG1().Pick(random.New())
	require.False(t, PHGR13Verify(setup.VK, qap, p2, s[:diff]))

	p2 = proof
	p2.gz = NewG1().Pick(random.New())
	require.False(t, PHGR13Verify(setup.VK, qap, p2, s[:diff]))

	p2 = proof
	vk2 := setup.VK
	vk2.bgamma2 = NewG2().Pick(random.New())
	require.False(t, PHGR13Verify(vk2, qap, p2, s[:diff]))

	p2 = proof
	vk2 = setup.VK
	vk2.av = NewG2().Pick(random.New())
	require.False(t, PHGR13Verify(vk2, qap, p2, s[:diff]))

	p2 = proof
	vk2 = setup.VK
	vk2.ay = NewG2().Pick(random.New())
	require.False(t, PHGR13Verify(vk2, qap, p2, s[:diff]))

}
