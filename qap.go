package playsnark

import (
	"github.com/drand/kyber/share"
)

// QAP is the polynomial version of a R1CS circuit. It contains the left right
// and outputs polynomials. There are exactly n polynomials for each, where n is
// the number of variables.
type QAP struct {
	// in pinocchio paper, this is the variable m
	// the paper decomposes this into input/output variables
	// with N = n + n' where n is input and n' is output
	// then we have all variables = {1...N, N+1 ... m} - it contains the
	// intermediate variables
	nbVars int
	// nbIO is N in the previous comment
	nbIO int
	// nbGates is the variable "d" in the pinochio paper
	nbGates int
	// left right and out
	left  []Poly
	right []Poly
	out   []Poly
	// z is the minimal polynomial (x-1)(x-2)(x-3)...
	z Poly
}

// ToQAP takes a R1CS circuit description and turns it into its polynomial QAP
// form.
// The first thing it does is transposing each of the matrix of the R1CS so each
// row represents the pairs of point that we want to interpolate.
// It thens transforms these rows into their finite field version and
// interpolate the polynomial.
func ToQAP(circuit R1CS) QAP {
	left := qapInterpolate(circuit.left)
	right := qapInterpolate(circuit.right)
	out := qapInterpolate(circuit.out)
	nbVar := len(circuit.vars)
	nbGates := len(circuit.left)
	var z Poly
	for i := 1; i <= nbGates; i++ {
		// you multiply by (x - i) with i being as high as the number of gates,
		// because the polynomials left,right and out vanishes on these inputs
		// -i + x
		xi := Poly([]Element{
			NewElement().Neg(Value(i).ToFieldElement()),
			one.Clone(),
		})
		if z == nil {
			z = xi
		} else {
			z = z.Mul(xi)
		}
	}
	return QAP{
		nbVars:  nbVar,
		nbGates: nbGates,
		nbIO:    circuit.nbIO(),
		left:    left,
		right:   right,
		out:     out,
		z:       z,
	}
}

func qapInterpolate(m Matrix) []Poly {
	// once transposed, a row of this matrix represents all the usage of this
	// variable at each step of the circuit
	// so for example, all the assignements on the lefts inputs look like this
	// in R1CS:
	//  c  x out v  u  w
	// [0, 1, 0, 0, 0, 0] g1
	// [0, 0, 0, 1, 0, 0] g2
	// [0, 1, 0, 0, 1, 0] g3
	// [5, 0, 0, 0, 0, 1] g4
	// so for example the second row of the transposed matrix represents all the
	// usage of the variable "x" as a left wire in the circuit -> [1,0,1,0]
	// Then we need to interpolate this as a polynomial such that:
	// p(1) = 1, p(2) = 0, p(3) = 1 and p(4) = 0
	// We return one polynomial for each variable
	t := m.Transpose()
	out := make([]Poly, 0, len(t))
	for _, variable := range t {
		var ys = make([]Element, 0, len(t))
		for _, v := range variable {
			ys = append(ys, v.ToFieldElement())
		}
		poly := Interpolate(ys)
		out = append(out, poly)
	}
	return out
}

func toPoly(s *share.PriPoly) Poly {
	pp := s.Coefficients()
	return Poly(pp)
}

// Verify takes the vector of solution and makes the dot product with it such
// that
// t = (Left . s) * (Right . s) - (Out . s)
// z = minimal polynomial = (x - 1) * (x - 2) .. (x - #ofgates)
// h = t/z with no remainder ! if there is no remainder, that means
// the polynomial t vanishes on all the points corresponding to the gate since
// z is a factor, hence the solution is correct
func (q *QAP) IsValid(sol Vector) bool {
	q.sanityCheck(sol)
	// We need to multiply each entry of the solution with the corresponding
	// polynomial.
	// The original python code is short and self explanatory:
	// for rval, a in zip(r, new_A):
	//  Apoly = add_polys(Apoly, multiply_polys([rval], a))
	//
	// Let me try to debunk that, by taking the left polynomials as example
	//  * There are 6 left polynomials, one for each variable, that set or unset
	//  the variable depending on the gate it is evaluated for.
	//  * The solution vector is of length 6, giving the value to set for each
	//  variable
	//  * We want to make a polynomial that evaluates all these polynomials at
	//  once that when "multiplied" (dot product) by the solution vector gives
	//  assign the right value to the right variable at the right gate, as in
	//  R1CS.
	//  * So we take the polynomial that corresponds to the first variable,
	//  multiply by the value from the solution vector, and add it to a
	//  acumulator polynomial.
	//  * We do the same for the second variable and we add it to the
	//  accumulator etc
	//  * That means when we evaluate solution at the first gate, so x=1, we
	//  will look diretly at all assignements of all variables in _one_ polynomial
	//  which is already "multiplied" by the solution vector
	//   ==> Ai(x) * si + Aj(x) * sj ...
	//   When evaluated to 1, first term will give the value of the const at the
	//   first gate, second term will give the value of the input "x" at the
	//   first gate etc leading to the same values as in R1CS.
	// Left one
	left, right, out := q.computeAggregatePoly(sol)
	// Now we are using these precedent polynomials in the satisfying equation
	// t(x) = left(x) * right(x) - out(x) = h(x)z(x)
	t := left.Mul(right).Sub(out)

	// now we try to verify if the equation is really satisfied by dividing eq /
	// z(x) : we should find no remainder
	_, rem := t.Div2(q.z)
	if len(rem.Normalize()) > 0 {
		return false
	}
	return true
}

func (q QAP) Quotient(sol Vector) Poly {
	// compute h(x) then evaluate it blindly at point s
	left, right, out := q.computeAggregatePoly(sol)
	// p(x) = t(x) * h(x)
	px := left.Mul(right).Sub(out)
	// h(x) = p(x) / t(x)
	hx, rem := px.Div2(q.z)
	if len(rem.Normalize()) > 0 {
		panic("apocalypse")
	}
	return hx
}

func (q QAP) computeAggregatePoly(sol Vector) (left Poly, right Poly, out Poly) {
	left = Poly([]Element{})
	right = Poly([]Element{})
	out = Poly([]Element{})
	for varIndex, val := range sol {
		polyVal := Poly([]Element{val.ToFieldElement()})
		left = left.Add(q.left[varIndex].Mul(polyVal))
		right = right.Add(q.right[varIndex].Mul(polyVal))
		out = out.Add(q.out[varIndex].Mul(polyVal))
	}
	return
}

func (q QAP) sanityCheck(sol Vector) {
	if len(sol) != len(q.left) {
		panic("different number of solution variables than left polynomials")
	}

	if len(sol) != len(q.right) {
		panic("different numberof solution variables than right polynomials")
	}

	if len(sol) != len(q.out) {
		panic("different numbers of solutions variables than out polynomials")
	}
}
