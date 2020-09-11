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
	vs []Commit
	ws []Commit
	ys []Commit
	// shifted by alpha a_v
	vas []Commit
	was []Commit
	yas []Commit

	// g^s^i
	gsi []Commit

	// g_v^(b * v_k(s))
	vbs []Commit
	wbs []Commit
	ybs []Commit
}

// VerificationKey contains the information from the trusted setup for the
// verifier to verify a valid proof
type VerificationKey struct {
	// g
	g1 Commit
	// g^a_v
	av Commit
	// g^a_w
	aw Commit
	// g^a_y
	ay Commit
	// g^gamma
	gamma Commit
	// g^(beta * gamma)
	bgamma Commit
	// g_y^(t(s)) where t is the minimal polynomial
	yts Commit

	// g^v_k(s) for all polynomials / variables (not just the input / output as
	// in the evaluation key)
	vs []Commit
	ws []Commit
	ys []Commit
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
	ek.gsi = generatePowersCommit(s, one.Clone(), qap.nbGates)
	// alpha for the left right and outputs are for generating the linear
	// combination in the exponent
	av := NewElement().Pick(random.New())
	aw := NewElement().Pick(random.New())
	ay := NewElement().Pick(random.New())
	// random element to form the basis of the commitments of al,al,a
	rv := NewElement().Pick(random.New())
	gv := NewCommit().Mul(rv, nil)
	rw := NewElement().Pick(random.New())
	gw := NewCommit().Mul(rw, nil)
	// ry = rv * rw
	ry := NewElement().Mul(rv, rw)
	gy := NewCommit().Mul(ry, nil)
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
	ek.vbs = generateEvalCommit(gv, qap.left[diff:], s, beta)
	ek.wbs = generateEvalCommit(gv, qap.right[diff:], s, beta)
	ek.ybs = generateEvalCommit(gv, qap.out[diff:], s, beta)

	gamma := NewElement().Pick(random.New())
	bgamma := NewElement().Mul(gamma, beta)

	vk.g1 = NewCommit()
	// g^alpha_v
	vk.av = NewCommit().Mul(av, nil)
	// g^alpha_w
	vk.aw = NewCommit().Mul(aw, nil)
	// g^alpha_y
	vk.ay = NewCommit().Mul(ay, nil)
	// g^gamma
	vk.gamma = NewCommit().Mul(gamma, nil)
	// g^(beta*gamma)
	vk.bgamma = NewCommit().Mul(bgamma, nil)
	// evaluation of the minimal polynomial at the unknonw index s
	ts := qap.z.Eval(s)
	// g_y^t(s)
	vk.yts = NewCommit().Mul(ts, nil)
	// g^(v_k(s)) for all k
	vk.vs = generateEvalCommit(gv, qap.left, s, one)
	vk.ws = generateEvalCommit(gw, qap.right, s, one)
	vk.ys = generateEvalCommit(gy, qap.out, s, one)
	return TrustedSetup{
		EK: ek,
		VK: vk,
	}
}

func generatePowersCommit(e Element, shift Element, power int) []Commit {
	var gi = make([]Commit, 0, power+1)
	gi = append(gi, NewCommit())
	var tmp = one.Clone()
	for i := 0; i <= power; i++ {
		tmp = tmp.Mul(tmp, e)
		tmp = tmp.Mul(tmp, shift)
		gi = append(gi, NewCommit().Mul(tmp, nil))
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
