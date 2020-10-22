package playsnark

import "github.com/drand/kyber/util/random"

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

type Groth16TrustedSetup struct {
	tw    groth16ToxicWaste
	Alpha G1
	Beta  G1
	Delta G1
	// {x^i} for i:0->nbGates-1
	Xi []G1
	// (beta*u_i(x) + alpha*v_i(x) + w_i(x)) / gamma for io related variable
	// on G1
	IoLP []G1
	// (beta*u_i(x) + alpha*v_i(x) + w_i(x)) / gamma for non-io / intermediate
	// related variable on G1
	NioLP []G1

	// XiT are { x^i * t(x) / delta } for i:0 -> nbGates-2 where t(x) is the
	// minimal polynomial of the QAP equation
	XiT []G1

	// Same element but in G2
	Beta2  G2
	Delta2 G2
	Gamma  G2
	// {x^i} for i:0->nbGates-1
	Xi2 []G1
}

func NewGroth16TrustedSetup(qap QAP) Groth16TrustedSetup {
	var tw groth16ToxicWaste
	var tr Groth16TrustedSetup
	tw.Alpha = NewElement().Pick(random.New())
	tr.Alpha = NewG1().Mul(tw.Alpha, nil)

	tw.Beta = NewElement().Pick(random.New())
	tr.Beta = NewG1().Mul(tw.Beta, nil)
	tr.Beta2 = NewG2().Mul(tw.Beta, nil)

	tw.Delta = NewElement().Pick(random.New())
	tr.Delta = NewG1().Mul(tw.Delta, nil)
	tr.Delta2 = NewG2().Mul(tw.Delta, nil)

	tw.X = NewElement().Pick(random.New())
	tr.Xi = generatePowersCommit(zeroG1, tw.X, one.Clone(), qap.nbGates)
	tr.Xi2 = generatePowersCommit(zeroG2, tw.X, one.Clone(), qap.nbGates)

	tw.Gamma = NewElement().Pick(random.New())
	tr.Gamma = NewG2().Mul(tw.Gamma, nil)
	// diff marks the separation between IO poly variables and intermediates
	// ones
	diff := qap.nbVars - qap.nbIO
	// (beta*u_i(x) + alpha*v_i(x) + w_i(x)) / gamma for io related variable
	// poly
	tw.IoLP, tr.IoLP = fullLinearPoly(qap, 0, diff, tw.X, tw.Alpha, tw.Beta, tw.Gamma)
	// same for intermediate variables, "non-io", and divided by delta
	tw.NioLP, tr.NioLP = fullLinearPoly(qap, diff, len(qap.right), tw.X, tw.Alpha, tw.Beta, tw.Delta)

	// XiT are { x^i * t(x) / delta } for i:0 -> nbGates-2 where t(x) is the
	tx := qap.z.Eval(tw.X)
	txd := NewElement().Div(tx, tw.Delta)
	power := len(qap.z) - 2
	tr.XiT = generatePowersCommit(zeroG1, tw.X, txd, power)

	tr.tw = tw
	return tr
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
