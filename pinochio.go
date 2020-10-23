package playsnark

import (
	"bytes"
	"fmt"

	"github.com/drand/kyber/util/random"
)

// This file implements the logic of the Pinocchio PHGR13 proof system mostly as
// described in the original paper https://eprint.iacr.org/2013/279.pdf with
// some modifications from https://eprint.iacr.org/2013/879.pdf that uses it in
// the asymmetric pairing case, notably with BLS12-381.

// PHGR13Setup contains the information needed for a prover and a verifier.
type PHGR13Setup struct {
	EK PHGR13EvalKey    // required by the prover
	VK PHGR13VerifKey   // required by the verifier
	t  pHGR13ToxicWaste // to be deleted - left for testing
}

// pHGR13ToxicWaste contains the information that should be deleted after a trusted
// setup but that we keep it here for testing reason.
type pHGR13ToxicWaste struct {
	s    Element
	beta Element
	rv   Element
	rw   Element
	ry   Element
	gv   G1
	gw   G2
	gy   G1
}

// PHGR13EvalKey contains the information of the trusted setup needed for the prove to
// create a valid proof
type PHGR13EvalKey struct {
	// These commit evaluation are only for the non IO variables
	// i.e. we take the polynomials that are evaluating only variables that are
	// not inputs to the circuit or outputs
	// g_v ^ (v_k(s))
	vs []G1
	ws []G2
	ys []G1
	// shifted by alpha a_v: we use these to verify that the prover correctly
	// uses the polynomials from the CRS
	vas []G1
	was []G1
	yas []G1

	// g^s^i for i=1...#ofgates (inclusive)
	// We use these to allow the prover to evaluate the polynomials blindly
	gsi []G1

	// g_v^(beta * v_k(s)) only for the non-IO variables
	// We use these to allow the prover to prove he used the same value for the
	// variable for all polynomials (i.e. in the exponent, he used v(s) and w(s)
	// and y(s) and not v(s') with some weird s' for example)
	vbs []G1
	wbs []G2
	ybs []G1
}

// PHGR13VerifKey contains the information from the trusted setup for the
// verifier to verify a valid proof
type PHGR13VerifKey struct {
	// g
	g1 G1
	// g^a_v
	av G2
	// g^a_w
	aw G2
	// g^a_y
	ay G2
	// g^gamma
	gamma G2
	// g^(beta * gamma) in G1
	bgamma G1
	// same but in G2 (optimization from https://eprint.iacr.org/2013/879.pdf
	bgamma2 G2
	// g_y^(t(s)) where t is the minimal polynomial (x-1)(x-2)... See QAP for
	// more details
	yts G2

	// g^v_k(s) for all polynomials (not just the input / output as
	// in the evaluation key)
	vs []G1
	ws []G2
	ys []G1
}

// https://eprint.iacr.org/2013/279.pdf
func NewPHGR13TrustedSetup(qap QAP) PHGR13Setup {
	var ek PHGR13EvalKey
	var vk PHGR13VerifKey
	// s is the private point at which the prover evaluates its QAP polynomials
	// s must thrown away after the trusted setup such that the prover doesn't
	// it, it only evaluates its polynomials blindly to this point
	s := NewElement().Pick(random.New())
	// gsi contains g^(s^i) from i=0 to g^(si^#of gates) included
	ek.gsi = GeneratePowersCommit(zeroG1, s, one.Clone(), qap.z.Degree()-2)
	// alpha for the left right and outputs are for generating the linear
	// combination in the exponent
	av := NewElement().Pick(random.New())
	aw := NewElement().Pick(random.New())
	ay := NewElement().Pick(random.New())
	// random element to form the basis of the commitments of al,al,a
	// that then form the different basis of evaluation of the three polynomials
	// g_v (in G1), g_w (in G2) and g_y (in G1)
	rv := NewElement().Pick(random.New())
	gv := NewG1().Mul(rv, nil)
	rw := NewElement().Pick(random.New())
	gw := NewG2().Mul(rw, nil)
	g1w := NewG1().Mul(rw, nil)
	// ry = rv * rw
	ry := NewElement().Mul(rv, rw)
	gy := NewG1().Mul(ry, nil)
	g2y := NewG2().Mul(ry, nil)
	// Compute the evaluation of the polynomials at the point s and their
	// shifted version. This is to make sure that prover indeed used a
	// the part of the CRS with a polynomial to build up its proof
	diff := qap.nbVars - qap.nbIO
	ek.vs = generateEvalCommit(gv, qap.left[diff:], s, one)
	ek.ws = generateEvalCommit(gw, qap.right[diff:], s, one)
	ek.ys = generateEvalCommit(gy, qap.out[diff:], s, one)
	// compute the same evaluation but shifted by their respective alpha
	ek.vas = generateEvalCommit(gv, qap.left[diff:], s, av)
	ek.was = generateEvalCommit(g1w, qap.right[diff:], s, aw)
	ek.yas = generateEvalCommit(gy, qap.out[diff:], s, ay)

	// Beta and gamma are used to check that same coefficients - same
	// polynomials - were used during the linear combination
	beta := NewElement().Pick(random.New())
	// we now evaluate the commitments of the polynomials shifted by beta
	ek.vbs = generateEvalCommit(gv, qap.left[diff:], s, beta)
	ek.wbs = generateEvalCommit(g1w, qap.right[diff:], s, beta)
	ek.ybs = generateEvalCommit(gy, qap.out[diff:], s, beta)

	gamma := NewElement().Pick(random.New())
	bgamma := NewElement().Mul(gamma, beta)

	vk.g1 = NewG1()
	// g^alpha_v
	vk.av = NewG2().Mul(av, nil)
	// g^alpha_w
	vk.aw = NewG1().Mul(aw, nil)
	// g^alpha_y
	vk.ay = NewG2().Mul(ay, nil)
	// g^gamma
	vk.gamma = NewG2().Mul(gamma, nil)
	// g^(beta*gamma)
	vk.bgamma = NewG1().Mul(bgamma, nil)
	vk.bgamma2 = NewG2().Mul(bgamma, nil)
	// evaluation of the minimal polynomial at the unknonw index s
	ts := qap.z.Eval(s)
	// t(s) * (r_y * G2)
	vk.yts = NewG2().Mul(ts, g2y)
	// g^(v_k(s)) for all k (input/output + intermediate)
	vk.vs = generateEvalCommit(gv, qap.left, s, one)
	vk.ws = generateEvalCommit(gw, qap.right, s, one)
	vk.ys = generateEvalCommit(gy, qap.out, s, one)
	return PHGR13Setup{
		EK: ek,
		VK: vk,
		t: pHGR13ToxicWaste{
			beta: beta,
			s:    s,
			gv:   gv,
			gw:   gw,
			gy:   gy,
			ry:   ry,
			rv:   rv,
			rw:   rw,
		},
	}
}

// PHGR13Proof contains all the information computed by a prover with the
// EvaluationKey and a circuit and its witness
type PHGR13Proof struct {
	// g^(SUM v_k(s)) for non-IO k
	vss G1
	// same this as vss but shifted by alpha_v: g^[alpha_v* (SUM * v_k(s))]
	vass G1
	// g^(SUM w_k(s)) for non-IO k
	wss G2
	// same this as wss but shifted by alpha_w: g^[alpha_w* (SUM * w_k(s))]
	// on G1, since we are in asymmetric pairing case
	// https://eprint.iacr.org/2013/879.pdf
	wass G1

	// g^(SUM y_k(s)) for non-IO k
	yss G1
	// same this as yss but shifted by alpha_y: g^[alpha_y* (SUM * y_k(s))]
	yass G1

	// g^h(s) where h is one of the product of p(x): p(x) = t(x) * h(x)
	hs G1

	// g^a_v * SUM v_k(s) * g^a_w * SUM w_k(s) * g^a_y * SUM y_k(s) for all non
	// IO - used during the linear check
	gz G1
}

// PHGR13Prove takes the evaluation key, the QAP polynomials and the solution
// vector and returns the corresponding proof.
func PHGR13Prove(ek PHGR13EvalKey, qap QAP, solution Vector) PHGR13Proof {
	// compute h(x) then evaluate it blindly at point s
	left, right, out := qap.computeAggregatePoly(solution)
	// p(x) = t(x) * h(x)
	px := left.Mul(right).Sub(out)
	// h(x) = p(x) / t(x)
	hx, rem := px.Div2(qap.z)
	if len(rem.Normalize()) > 0 {
		panic("apocalypse")
	}
	// g^h(s) = SUM(s_i * G)
	ghs := hx.BlindEval(zeroG1, ek.gsi)
	diff := qap.nbVars - qap.nbIO
	// compute g_v^(SUM v_k(s) * sol[k]) for k being NON IO
	// same for y and w
	computeSolCommit := func(zero Commit, evalCommit []Commit) Commit {
		var tmp = zero.Clone()
		var acc = zero.Clone()
		for i, ec := range evalCommit {
			acc = acc.Add(acc, tmp.Mul(solution[diff+i].ToFieldElement(), ec))
		}
		return acc
	}
	// g^(SUM sol[k] * v_k(s))
	gvmids := computeSolCommit(zeroG1, ek.vs)
	gwmids := computeSolCommit(zeroG2, ek.ws)
	gymids := computeSolCommit(zeroG1, ek.ys)
	gvamids := computeSolCommit(zeroG1, ek.vas)
	gwamids := computeSolCommit(zeroG1, ek.was)
	gyamids := computeSolCommit(zeroG1, ek.yas)

	// g^(SUM sol[k] * v_k(s) * beta)
	gvbmids := computeSolCommit(zeroG1, ek.vbs)
	gwbmids := computeSolCommit(zeroG1, ek.wbs)
	gybmids := computeSolCommit(zeroG1, ek.ybs)
	gz := NewG1().Add(gvbmids, NewG1().Add(gwbmids, gybmids))

	return PHGR13Proof{
		hs:   ghs,
		vss:  gvmids,
		wss:  gwmids,
		yss:  gymids,
		vass: gvamids,
		wass: gwamids,
		yass: gyamids,
		gz:   gz,
	}
}

func (p *PHGR13Proof) String() string {
	var b bytes.Buffer
	buff, _ := p.hs.MarshalBinary()
	b.WriteString(fmt.Sprintf("hs: %x\n", buff))
	buff, _ = p.vss.MarshalBinary()
	b.WriteString(fmt.Sprintf("vss: %x\n", buff))
	buff, _ = p.wss.MarshalBinary()
	b.WriteString(fmt.Sprintf("wss: %x\n", buff))
	buff, _ = p.yss.MarshalBinary()
	b.WriteString(fmt.Sprintf("yss: %x\n", buff))
	buff, _ = p.vass.MarshalBinary()
	b.WriteString(fmt.Sprintf("vass: %x\n", buff))
	buff, _ = p.wass.MarshalBinary()
	b.WriteString(fmt.Sprintf("wass: %x\n", buff))
	buff, _ = p.yass.MarshalBinary()
	b.WriteString(fmt.Sprintf("yass: %x\n", buff))
	buff, _ = p.gz.MarshalBinary()
	b.WriteString(fmt.Sprintf("gz: %x\n", buff))
	return b.String()
}

// PHGR13Verify runs the different consistency checks. It needs the trusted
// setup variables (actually only the verification key), the QAP describing
// the circuit, the proof generated by the prover and the inputs and outputs
// expected
func PHGR13Verify(vk PHGR13VerifKey, qap QAP, p PHGR13Proof, io Vector) bool {
	// DIVISION CHECK: we look if the prover correctly evaluated the polynomials
	// such that the QAP equation resolves:
	// r_v * r_w * v(s) * w(s) == r_y * (p(s) +  y(s))
	// r_y * v(s) * w(s)       == r_y * (p(s) + y(s))
	// v(s) * w(s) - y(s)      == p(s) == h(s) * t(s)
	// which is the QAP equation

	// g^v_io(s)^ck where ck are the "valid" coefficients since they're the
	// inputs
	diff := qap.nbVars - qap.nbIO
	{
		gvkio := computeCommitIOSolution(zeroG1, vk.vs[:diff], io)
		gwkio := computeCommitIOSolution(zeroG2, vk.ws[:diff], io)
		gykio := computeCommitIOSolution(zeroG1, vk.ys[:diff], io)

		// first term is reconstructed above from the verification key and the
		// input/output, second term is given by the prover (the intermediate
		// values of the circuit)
		// gv = g^(r_v * SUM_io v_k(s) * c_k) + g^(r_v * SUM_nio v_k(s) * c_k)
		// gv = g^(r_v * SUM_all v_k(s) * c_k)
		// where io represents the indices of the inputs / outputs variables
		// and nio the rest
		gv := NewG1().Add(gvkio, p.vss)
		// g_w_io(s) * g_w_mid(s)
		gw := NewG2().Add(gwkio, p.wss)
		// g_y_io(s) * g_y_mid(s)
		gy := NewG1().Add(gykio, p.yss)

		// e(g1^(r_v * v(s)),g2^(r_w * w(s)) =
		// e(g1,g2)^(r_v * v(s) * r_w * w(s))
		left := Pair(gv, gw)
		// e(g^t(s) , g^h(s)) = e(g1,g2)^(t(s) * h(s))
		//					  = e(g1,g2)^p(s)
		right1 := Pair(p.hs, vk.yts)
		// e(g1^(r_y * y(s)), g2)
		right2 := Pair(gy, zeroG2.Clone().Base())
		right := zeroGT.Clone().Add(right1, right2)
		if !left.Equal(right) {
			return false
		}
	}
	{
		// CRS CHECK: We check that the prover correctly used the
		// g^v(s) from the CRS to form its proof, that they're the polynomials
		// from the CRS
		// To do this we have the "shifted" version with a_v a_w and a_y and thanks
		// to the q power knowledge of exponent assumption, we can tell that the
		// prover is unable to pick a polynomial v' instead of v where the  the
		// following equation pass:
		// e(g^a_v * v(s),g) = e(g^v(s),g^a_v)
		// g^a_v * v(s), g^v(s) term is provided by prover, g^a_v by the CRS
		left := Pair(p.vass, NewG2().Base())
		right := Pair(p.vss, vk.av)
		if !left.Equal(right) {
			return false
		}
		// we do the same for W except here we computed the evaluation of g^w(s) on
		// G2 so we switch the argument order
		left = Pair(p.wass, NewG2().Base())
		right = Pair(vk.aw, p.wss)
		if !left.Equal(right) {
			return false
		}
		// we do the same for Y as in V
		left = Pair(p.yass, NewG2().Base())
		right = Pair(p.yss, vk.ay)
		if !left.Equal(right) {
			return false
		}
	}
	{
		// LINEAR CHECK: we also need to make sure that prover used the same valud
		// into the variable for the three polys v w and y. Wo do this via a linear
		// sum check with randomnized coefficient for each polynomial to avoid
		// malleability.
		// left is e(g^(r_v*v(s) + r_w*w(s) + r_y *y(s)) * beta,g^gamma)
		left := Pair(p.gz, vk.gamma)
		// right is splitted into two GT elements as an optimization of
		// https://eprint.iacr.org/2013/879.pdf because of polynomial w which is in
		// g2
		// t1 := e(g^(r_v*v(s)) * g^(r_y*y(s)), g^beta*gamma)
		// t2 := e(g^beta*gamma, g^(r_w*w(s))
		// right = t1*t2 = left !
		lt1 := NewG1().Add(p.vss, p.yss)
		t1 := Pair(lt1, vk.bgamma2)
		t2 := Pair(vk.bgamma, p.wss)
		right := zeroGT.Clone().Add(t1, t2)
		if !right.Equal(left) {
			fmt.Println("damdam")
			return false
		}
	}
	return true
}

// Evaluate all the polynomials at the given x and return the list of their
// commitments using the base and shifted by "shift". Note by giving shift as
// one, nothing happens !
// { g^(shift * p_i(x)) } for all p_i given
func generateEvalCommit(base Commit, polys []Poly, x, shift Element) []Commit {
	var ret = make([]Commit, 0, len(polys))
	for i := range polys {
		tmp := polys[i].Normalize().Eval(x)
		ret = append(ret, base.Clone().Mul(tmp.Mul(tmp, shift), base))
	}
	return ret
}

func computeCommitIOSolution(base Commit, poly []Commit, io Vector) Commit {
	var acc = base.Clone().Null()
	var tmp = base.Clone()
	for i, gs := range poly {
		// commitment of the evaluation of the polynomial (v,w, and y) to
		// the power of the coefficient of the solution vector (just the
		// public part, input / output)
		// We add all these sums to give back the commitment of the
		// evaluation of the aggregated polynomial as in the QAP case to a
		// random point s
		// (g^v_k(s)) ^ c_k = g^(v_k(s) * c_k)
		// We then compute g^SUM(v_k(s) * c_k) which is equal
		// SUM [v_k(s) * c_k * G] = [SUM v_k(s) * c_k] * G
		tmp := tmp.Mul(io[i].ToFieldElement(), gs)
		acc = acc.Add(acc, tmp)
	}
	return acc
}
