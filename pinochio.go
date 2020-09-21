package playsnark

import (
	"github.com/drand/kyber/util/random"
)

// TrustedSetup contains the information needed for a prover and a verifier.
type TrustedSetup struct {
	EK EvalKey
	VK VerificationKey
	t  ToxicWaste
}

// ToxicWaste contains the information that should be deleted after a trusted
// setup but that we keep it here for testing reason.
type ToxicWaste struct {
	s  Element
	ry Element
	rv Element
	rw Element
	gv G1
	gw G2
	gy G1
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
	was []G1
	yas []G1

	// g^s^i for i=1...#ofgates (inclusive)
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
	av G2
	// g^a_w
	aw G2
	// g^a_y
	ay G2
	// g^gamma
	gamma G2
	// g^(beta * gamma) in G1
	bgamma G1
	// g_y^(t(s)) where t is the minimal polynomial
	yts G2

	// g^v_k(s) for all polynomials (not just the input / output as
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
	ek.ys = generateEvalCommit(gy, qap.right[diff:], s, one)
	// compute the same evaluation but shifted by their respective alpha
	ek.vas = generateEvalCommit(gv, qap.left[diff:], s, av)
	ek.was = generateEvalCommit(g1w, qap.right[diff:], s, aw)
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
	vk.av = NewG2().Mul(av, nil)
	// g^alpha_w
	vk.aw = NewG1().Mul(aw, nil)
	// g^alpha_y
	vk.ay = NewG2().Mul(ay, nil)
	// g^gamma
	vk.gamma = NewG1().Mul(gamma, nil)
	// g^(beta*gamma)
	vk.bgamma = NewG1().Mul(bgamma, nil)
	// evaluation of the minimal polynomial at the unknonw index s
	ts := qap.z.Eval(s)
	// t(s) * (r_y * G2)
	vk.yts = NewG2().Mul(ts, g2y)
	// g^(v_k(s)) for all k (input/output + intermediate)
	vk.vs = generateEvalCommit(gv, qap.left, s, one)
	vk.ws = generateEvalCommit(gw, qap.right, s, one)
	vk.ys = generateEvalCommit(gy, qap.out, s, one)
	return TrustedSetup{
		EK: ek,
		VK: vk,
		t: ToxicWaste{
			s:  s,
			gv: gv,
			gw: gw,
			gy: gy,
			ry: ry,
			rv: rv,
			rw: rw,
		},
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
	wss G1
	// same this as wss but shifted by alpha_w: g^[alpha_w* (SUM * w_k(s))]
	wass G1

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
	// g^h(s) = SUM(s_i * G
	ghs := hx.BlindEval(zeroG1, setup.EK.gsi)
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
	gvmids := computeSolCommit(zeroG1, setup.EK.vs)
	gwmids := computeSolCommit(zeroG2, setup.EK.ws)
	gymids := computeSolCommit(zeroG1, setup.EK.ys)
	gvamids := computeSolCommit(zeroG1, setup.EK.vas)
	gwamids := computeSolCommit(zeroG1, setup.EK.was)
	gyamids := computeSolCommit(zeroG1, setup.EK.yas)

	return Proof{
		hs:   ghs,
		vss:  gvmids,
		wss:  gwmids,
		yss:  gymids,
		vass: gvamids,
		wass: gwamids,
		yass: gyamids,
	}
}

// VerifyProof runs the different consistency checks. It needs the trusted
// setup variables (actually only the verification key), the QAP describing
// the circuit, the proof generated by the prover and the inputs and outputs
// expected
func VerifyProof(vk VerificationKey, qap QAP, p Proof, io Vector) bool {
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
		// gv = (SUM_io v_k(s) * c_k) * Gv + (SUM_nio v_k(s) * c_k) * Gv
		// gv = (SUM_all v_k(s) * c_k) * Gv
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
	// LINEAR COMBINATION CHECK: We check that the prover correctly used the
	// g^v(s) etc to form its proof and not any other polynomial by forcing him
	// to do a random linear combination with coefficients "glued" with the
	// polynomials (v, w, y) in the CRS.
	// --> we check that : e(g^a_v * v(s),g) = e(g^v(s),g^a_v)
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
	return true
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
func generateEvalCommit(base Commit, polys []Poly, x, shift Element) []Commit {
	var ret = make([]Commit, 0, len(polys))
	for i := range polys {
		tmp := polys[i].Eval(x)
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
