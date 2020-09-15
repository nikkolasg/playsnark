package playsnark

import "github.com/drand/kyber/util/random"

// TrustedSetup contains the information needed for a prover and a verifier.
type TrustedSetup struct {
	EK EvalKey
	VK VerificationKey
}

// EvalKey contains the information of the trusted setup needed for the prove to
// create a valid proof
type EvalKey struct {
	// These commit evaluation are only for the non IO variables
	// i.e. we take the polynomials that are evaluating only variables that are
	// not inputs to the function or outputs
	// g_v ^ (v_k(s))
	vs []G1
	ws []G2
	ys []G1
	// shifted by alpha a_v
	vas []G1
	was []G2
	yas []G1

	// g^s^i
	gsi []G1

	// g_v^(b * v_k(s))
	vbs []G1
	wbs []G2
	ybs []G1
}

// VerificationKey contains the information from the trusted setup for the
// verifier to verify a valid proof
type VerificationKey struct {
	// g
	g1 G1
	// g^a_v
	av G1
	// g^a_w
	aw G2
	// g^a_y
	ay G1
	// g^gamma
	gamma G2
	// g^(beta * gamma) in G1
	bgamma G1
	// g_y^(t(s)) where t is the minimal polynomial
	yts G2

	// g^v_k(s) for all polynomials / variables (not just the input / output as
	// in the evaluation key)
	vs []G1
	ws []G2
	ys []G1
}

// https://eprint.iacr.org/2013/279.pdf
func NewTrustedSetup(qap QAP) TrustedSetup {
	var ek EvalKey
	var vk VerificationKey
	// s is the private point at which the prover evaluates its QAP polynomials
	// s must thrown away after the trusted setup such that the prover doesn't
	// it, it only evaluates its polynomials blindly to this point
	s := NewElement().Pick(random.New())
	// gsi contains g^(s^i) from i=0 to g^(si^#of gates) included
	ek.gsi = generatePowersCommit(zeroG1, s, one.Clone(), qap.nbGates)
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
	// ry = rv * rw
	ry := NewElement().Mul(rv, rw)
	gy := NewG1().Mul(ry, nil)
	// Compute the evaluation of the polynomials at the point s and their
	// shifted version. This is to make sure that prover indeed used a
	// the part of the CRS with a polynomial to build up its proof
	diff := qap.nbVars - qap.nbIO
	ek.vs = generateEvalCommit(gv, qap.left[diff:], s, one)
	ek.ws = generateEvalCommit(gw, qap.right[diff:], s, one)
	ek.ys = generateEvalCommit(gy, qap.right[diff:], s, one)
	// compute the same evaluation but shifted by their respective alpha
	ek.vas = generateEvalCommit(gv, qap.left[diff:], s, av)
	ek.was = generateEvalCommit(gw, qap.right[diff:], s, aw)
	ek.yas = generateEvalCommit(gy, qap.right[diff:], s, ay)

	// Beta and gamma are used to check that same coefficients - same
	// polynomials - were used during the linear combination
	beta := NewElement().Pick(random.New())
	// we now evaluate the commitments of the polynomials shifted by beta
	ek.vbs = generateEvalCommit(gv, qap.left[diff:], s, beta)
	ek.wbs = generateEvalCommit(gv, qap.right[diff:], s, beta)
	ek.ybs = generateEvalCommit(gv, qap.out[diff:], s, beta)

	gamma := NewElement().Pick(random.New())
	bgamma := NewElement().Mul(gamma, beta)

	vk.g1 = NewG1()
	// g^alpha_v
	vk.av = NewG1().Mul(av, nil)
	// g^alpha_w
	vk.aw = NewG2().Mul(aw, nil)
	// g^alpha_y
	vk.ay = NewG1().Mul(ay, nil)
	// g^gamma
	vk.gamma = NewG1().Mul(gamma, nil)
	// g^(beta*gamma)
	vk.bgamma = NewG1().Mul(bgamma, nil)
	// evaluation of the minimal polynomial at the unknonw index s
	ts := qap.z.Eval(s)
	// g_y^t(s)
	vk.yts = NewG1().Mul(ts, nil)
	// g^(v_k(s)) for all k AND g^(v_0(s)) where v(0) = 1 so it's equal to g
	vk.vs = append([]Commit{gv}, generateEvalCommit(gv, qap.left, s, one)...)
	vk.ws = append([]Commit{gw}, generateEvalCommit(gw, qap.right, s, one)...)
	vk.ys = append([]Commit{gy}, generateEvalCommit(gy, qap.out, s, one)...)
	return TrustedSetup{
		EK: ek,
		VK: vk,
	}
}

// Proof contains all the information computed by a prover with the
// EvaluationKey and a circuit and its witness
type Proof struct {
	// g^(SUM v_k(s)) for non-IO k
	vss G1
	// same this as vss but shifted by alpha_v: g^[alpha_v* (SUM * v_k(s))]
	vass G1
	// Same things as vss but shifted by beta
	vbss Commit
	// g^(SUM w_k(s)) for non-IO k
	wss G2
	// same this as wss but shifted by alpha_w: g^[alpha_w* (SUM * w_k(s))]
	wass G2

	// g^(SUM y_k(s)) for non-IO k
	yss G1
	// same this as yss but shifted by alpha_y: g^[alpha_y* (SUM * y_k(s))]
	yass G1

	// g^h(s) where h is one of the product of p(x): p(x) = t(x) * h(x)
	hs G1
}

func GenProof(setup TrustedSetup, qap QAP, solution Vector) Proof {
	// compute h(x) then evaluate it blindly at point s
	left, right, out := qap.computeAggregatePoly(solution)
	// p(x) = t(x) * h(x)
	px := left.Mul(right).Sub(out)
	// h(x) = p(x) / t(x)
	hx, rem := px.Div2(qap.z)
	if len(rem.Normalize()) > 0 {
		panic("apocalypse")
	}
	// g^h(s)
	ghs := hx.BlindEval(zeroG1, setup.EK.gsi)
	diff := qap.nbVars - qap.nbIO
	// compute g^(SUM v_k(s) * sol[k]) for k being NON IO
	computeSolCommit := func(zero Commit, evalCommit []Commit) Commit {
		var tmp = zero.Clone()
		var acc = zero.Clone()
		for i, ec := range evalCommit {
			acc = acc.Add(acc, tmp.Mul(solution[diff+i].ToFieldElement(), ec))
		}
		return acc
	}
	// g^(SUM sol[k] * v_k(s))
	gvmids := computeSolCommit(zeroG1, setup.EK.vs)
	gwmids := computeSolCommit(zeroG2, setup.EK.ws)
	gymids := computeSolCommit(zeroG1, setup.EK.ys)
	return Proof{
		hs:  ghs,
		vss: gvmids,
		wss: gwmids,
		yss: gymids,
	}
}

// VerifyProof runs the different consistency checks. It needs the trusted
// setup variables (actually only the verification key), the QAP describing
// the circuit, the proof generated by the prover and the inputs and outputs
// expected
func VerifyProof(vk VerificationKey, qap QAP, p Proof, io Vector) bool {
	// g^v_io(s)^ck where ck are the "valid" coefficients since they're the
	// inputs
	diff := qap.nbVars - qap.nbIO
	evalCommit := func(base Commit, poly []Commit) Commit {
		var acc = base.Clone().Null()
		for i, gs := range vk.vs[:diff] {
			// commitment of the evaluation of the polynomial (v,w, and y) to
			// the power of the coefficient of the solution vector (just the
			// public part, input / output)
			// We add all these sums to give back the commitment of the
			// evaluation of the aggregated polynomial as in the QAP case to a
			// random point s
			// (g^v_k(s)) ^ c_k
			tmp := base.Clone().Mul(io[i].ToFieldElement(), gs)
			acc = acc.Add(acc, tmp)
		}
		return acc
	}
	gvkio := evalCommit(zeroG1, vk.vs[:diff])
	gwkio := evalCommit(zeroG2, vk.ws[:diff])
	gykio := evalCommit(zeroG2, vk.ys[:diff])

	// g_v_io(s) * g_v_mid(s)
	// first term is reconstructed above from the verification key and the
	// input/output, second term is given by the prover (the intermediate
	// values of the circuit)
	gv := NewG1().Add(gvkio, p.vss)
	// g_w_io(s) * g_w_mid(s)
	gw := NewG2().Add(gwkio, p.wss)
	// g_y_io(s) * g_y_mid(s)
	gy := NewG1().Add(gykio, p.yss)

	// NOTE: we put all elements on G1 so during the pairing, the G2 element is
	// always G_2
	// e(g1^(r_v * v(s)),g2^(r_w * w(s)) =
	// e(g1,g2)^(r_v * v(s) * r_w * w(s))
	left := Pair(gv, gw)
	// e(g^t(s) , g^h(s)) = e(g1,g2)^(t(s) * h(s))
	//					  = e(g1,g2)^p(s)
	right1 := Pair(vk.yts, p.hs)
	// e(g1^(r_y * y(s)), g2)
	right2 := Pair(gy, zeroG2.Clone().Base())
	right := zeroGT.Clone().Add(right1, right2)
	return left.Equal(right)
}

// generatePowersCommit returns { shift * s^i} for i=0...power included
func generatePowersCommit(base Commit, e Element, shift Element, power int) []Commit {
	var gi = make([]Commit, 0, power+1)
	gi = append(gi, base.Clone().Base())
	var tmp = one.Clone()
	for i := 0; i < power; i++ {
		tmp = tmp.Mul(tmp, e)
		tmp = tmp.Mul(tmp, shift)
		gi = append(gi, base.Clone().Mul(tmp, nil))
	}
	return gi
}

// Evaluate all the polynomials at the given x and return the list of their
// commitments using the base and shifted by "shift". Note by giving shift as
// one, nothing happens !
// { g^(shift * p_i(x)) } for all p_i given
//
func generateEvalCommit(base Commit, polys []Poly, x, shift Element) []Commit {
	var ret = make([]Commit, 0, len(polys))
	for i := range polys {
		tmp := polys[i].Eval(x)
		ret = append(ret, base.Clone().Mul(tmp.Mul(tmp, shift), base))
	}
	return ret
}
