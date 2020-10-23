package playsnark

import (
	"github.com/drand/kyber/util/random"
)

// Implements Groth16 paper https://eprint.iacr.org/2016/260.pdf
// In the paper, m represents the number of variables and n the number of
// constraints / equation. Note this implementation does not perform all
// optimization listed in the paper such as separating the verification from
// evaluation key and the pre-pairing result from the trusted setup.

// groth16ToxicWaste contains the results that must be delete after a trusted
// setup. It is kept here for testing and learning purpose.
type groth16ToxicWaste struct {
	// Elements that are going to be committed on G1 for the "public" part of the
	// trusted setup
	Alpha Element
	Beta  Element
	Delta Element
	X     Element
	IoLP  []Element
	NioLP []Element

	Gamma Element
}

// Groth16Setup contains all the information created during a trusted setup
// required by the prover and the verifier.
type Groth16Setup struct {
	tw groth16ToxicWaste
	// Alpha and beta are required to make sure the computation of the proof
	// elements A B and C are consistent with each other w.r.t. the intermediate
	// variables used, i.e. they used the same a_i inside their computation.
	// three computations
	Alpha G1
	Beta  G1
	// Delta and gamma (in G2 later) forces independence of computations for A
	// and B such that results can only be balanced by C and nothing else.
	Delta G1
	// {x^i} for i:0->nbGates-1
	Xi []G1
	// (beta*u_i(x) + alpha*v_i(x) + w_i(x)) / gamma for io related variable
	// on G1
	IoLP []G1
	// (beta*u_i(x) + alpha*v_i(x) + w_i(x)) / delta for non-io / intermediate
	// related variable on G1
	NioLP []G1

	// XiT are { x^i * t(x) / delta } for i:0 -> nbGates-2 where t(x) is the
	// minimal polynomial of the QAP equation - used by the prover to evaluate
	// h(x) * t(x) / delta blindly
	XiT []G1

	// Same element but in G2
	Beta2  G2
	Delta2 G2
	Gamma  G2
	// {x^i} for i:0->nbGates-1
	Xi2 []G1
}

// NewGroth16TrustedSetup returns a setup for the given circuit
func NewGroth16TrustedSetup(qap QAP) Groth16Setup {
	var tw groth16ToxicWaste
	var tr Groth16Setup
	tw.Alpha = NewElement().Pick(random.New())
	tr.Alpha = NewG1().Mul(tw.Alpha, nil)

	tw.Beta = NewElement().Pick(random.New())
	tr.Beta = NewG1().Mul(tw.Beta, nil)
	tr.Beta2 = NewG2().Mul(tw.Beta, nil)

	tw.Delta = NewElement().Pick(random.New())
	tr.Delta = NewG1().Mul(tw.Delta, nil)
	tr.Delta2 = NewG2().Mul(tw.Delta, nil)

	tw.X = NewElement().Pick(random.New())
	tr.Xi = GeneratePowersCommit(zeroG1, tw.X, one.Clone(), qap.nbGates-1)
	tr.Xi2 = GeneratePowersCommit(zeroG2, tw.X, one.Clone(), qap.nbGates-1)

	tw.Gamma = NewElement().Pick(random.New())
	tr.Gamma = NewG2().Mul(tw.Gamma, nil)
	// diff marks the separation between IO poly variables and intermediates
	// ones
	diff := qap.nbVars - qap.nbIO
	// (beta*u_i(x) + alpha*v_i(x) + w_i(x)) / gamma for io related variable
	// poly
	tw.IoLP, tr.IoLP = fullLinearPoly(qap, 0, diff, tw.X, tw.Alpha, tw.Beta, tw.Gamma)
	// same for intermediate variables, "non-io", and divided by delta
	tw.NioLP, tr.NioLP = fullLinearPoly(qap, diff, qap.nbVars, tw.X, tw.Alpha, tw.Beta, tw.Delta)

	// XiT are { x^i * t(x) / delta } for i:0 -> nbGates-2 where t(x) is the
	tx := qap.z.Eval(tw.X)
	txd := NewElement().Div(tx, tw.Delta)
	power := qap.nbGates - 2
	tr.XiT = GeneratePowersCommit(zeroG1, tw.X, txd, power)

	tr.tw = tw
	return tr
}

// Groth16Proof contains the three elements required by the verifier as well as
// the private values used by the prover (which must not be given to verifier as
// this would break zero knowledge)
type Groth16Proof struct {
	tp groth16ToxicProof
	A  G1
	B  G2
	C  G1
}

// This struct contains the two private fields randomly sampled by the prover
// during computation of the proof, required to provide zero knowledge.
type groth16ToxicProof struct {
	R Element
	S Element
}

// Groth16Prove proofs it knows a solution sol for the given circuit and returns
// the proof
func Groth16Prove(tr Groth16Setup, q QAP, sol Vector) Groth16Proof {
	// The proof code is structured in three pieces, for generating the three
	// elements of the proofs A B and C.
	//
	// the following is just an helper function to compute sum of blindly
	// evaluated polynomials
	// basis^SUM(a_i * poly(x))
	// For SUM(a_i * u_i(x)) since we dont know x, we need to use G1^x^i
	// directly
	// g^u_i(x) = SUM_j( g^(x^i)^coeff(u_i,j)) = g^(u0 *x^0 + u1*x^1 + ...)
	// We then multiply by the solution vector to have g^(s_i * u_i(x))
	// and iteratively sum the results for all variables
	sumBlind := func(basis Commit, polys []Poly, xi []Commit) Commit {
		var sum = basis.Clone().Null()
		for i := 0; i < q.nbVars; i++ {
			uix := polys[i].BlindEval(basis.Clone().Null(), xi)
			sum = sum.Add(sum, uix.Mul(sol[i].ToFieldElement(), uix))
		}
		return sum
	}
	// ----------------------------------------------------
	// Compute A = G1^(alpha + SUM(a_i * u_i(x)) + r*delta)
	// we compute each part directly in the exponent thx to the trusted setup
	//
	var A = sumBlind(zeroG1, q.left, tr.Xi)
	//  Pick r and then compute g^(r * delta)
	r := NewElement().Pick(random.New())
	rd := NewG1().Mul(r, tr.Delta)
	// A = G1^(alpha + SUM(a_i * u_i(x)) + r*delta)
	A = A.Add(A, rd)
	A = A.Add(tr.Alpha, A)

	// ----------------------------------------------
	// We do something similar for B expcet in it's G2
	// B = G2^(beta + SUM(a_i * v_i(x)) + s*delta
	var B = sumBlind(zeroG2, q.right, tr.Xi2)
	s := NewElement().Pick(random.New())
	sd := NewG2().Mul(s, tr.Delta2)
	B = B.Add(B, sd)
	B = B.Add(tr.Beta2, B)

	// ------------------------------------------------------------------
	// C is a bit more complex and we need to use different parts of the trusted
	// setup to construct it
	// C = G1^((1/delta * SUM(a_i * NioLP) + h(x)t(x)) + A*s + B*r - r*s*delta
	// where NioLP = beta*u_i(x) + alpha * v_i(x) + w_i(x) where we only sum
	// over the intermediates / non IO variables.
	//
	C := NewG1().Null()
	// for the part with NioLP we use the NioLP part of the trusted setup and
	// multiply every entry by the piecewise solution element
	nio := NewG1().Null()
	// we only take variables which are _not_ io
	diff := q.nbVars - q.nbIO
	for i := range tr.NioLP {
		nio = nio.Add(nio, NewG1().Mul(sol[i+diff].ToFieldElement(), tr.NioLP[i]))
	}
	C = C.Add(C, nio)
	// we can compute h(x)t(x)/delta from the XiT part of the trusted setup
	// We can construct h(x) thanks to x^i and since we want to multiply by t(x)
	// and divide by delta, then we directly use x^i * t(x) / delta which is XiT
	// we first compute polynomial h so we get the coefficients
	h := q.Quotient(sol)
	htd := h.BlindEval(zeroG1, tr.XiT)
	C = C.Add(C, htd)

	//  As is simple multiplication
	As := NewG1().Mul(s, A)
	C = C.Add(C, As)
	// Br forces us to recompute B in G1 group though
	B1 := sumBlind(zeroG1, q.right, tr.Xi)
	sd1 := NewG1().Mul(s, tr.Delta)
	B1 = B1.Add(B1, sd1)
	B1 = B1.Add(B1, tr.Beta)
	Br := NewG1().Mul(r, B1)
	C = C.Add(C, Br)
	//  - r*s*delta
	rsd := NewG1().Mul(NewElement().Mul(r, s), tr.Delta)
	C = C.Add(C, rsd.Neg(rsd))

	return Groth16Proof{
		tp: groth16ToxicProof{
			R: r,
			S: s,
		},
		A: A,
		B: B,
		C: C,
	}
}

// Groth16Verify returns true if the proof is valid
func Groth16Verify(tr Groth16Setup, q QAP, p Groth16Proof, io Vector) bool {
	// Proof verification consists in 4 pairings (without optimizations) and one
	// equation check:
	// left side :  e(A * B)
	left := Pair(p.A, p.B)
	// right side: a * b * c
	//  	a. e(alpha, beta)
	//		b. e(SUM IoLP, gamma)
	//		c. e(C1,  delta)

	a := Pair(tr.Alpha, tr.Beta2)
	b1 := NewG1().Null()
	for i, iolp := range tr.IoLP {
		b1 = b1.Add(b1, NewG1().Mul(io[i].ToFieldElement(), iolp))
	}
	b := Pair(b1, tr.Gamma)
	c := Pair(p.C, tr.Delta2)
	right := a.Add(a, b.Add(b, c))
	return left.Equal(right)
}

// in g1, (beta*u_i(x) + alpha*v_i(x) + w_i(x)) for the i-th poly variable.
// I call this relation "linearPoly". fullLinearPoly iterates over multiples
// variables and returns the list and its commitment
func linearPolyForVar(qap QAP, i int, x, alpha, beta Element) Element {
	// u_i(x)
	ui := qap.left[i].Eval(x)
	// beta * u_i(x)
	bui := NewElement().Mul(ui, beta)
	// v_i(x)
	vi := qap.right[i].Eval(x)
	// alpha * v_i(x)
	avi := NewElement().Mul(vi, alpha)
	wi := qap.out[i].Eval(x)
	return wi.Add(wi, NewElement().Add(bui, avi))
}

// call linearPolyForVar for variable between min and max, and use the "div" as
// divider. Div should be gamma for the IO related variables and delta for the
// non IO related variables
func fullLinearPoly(qap QAP, min, max int, x, alpha, beta, div Element) ([]Element, []G1) {
	var length = max - min
	var lps = make([]Element, 0, length)
	var commitLps = make([]G1, 0, length)
	for i := min; i < max; i++ {
		var lp = NewElement().Div(linearPolyForVar(qap, i, x, alpha, beta), div)
		lps = append(lps, lp)
		commitLps = append(commitLps, NewG1().Mul(lp, nil))
	}
	return lps, commitLps
}
