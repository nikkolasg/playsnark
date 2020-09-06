package playsnark

import "github.com/drand/kyber/share"

// QAP is the polynomial version of a R1CS circuit. It contains the left right
// and outputs polynomials. There are exactly n polynomials for each, where n is
// the number of variables.
type QAP struct {
	left  []Poly
	right []Poly
	out   []Poly
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
	return QAP{
		left:  left,
		right: right,
		out:   out,
	}
}

type pair = share.PriShare

func qapInterpolate(m Matrix) []Poly {
	// once transposed, a row of this matrix represents all the usage of this
	// variable at each step of the circuit
	// so for example, all the assignements on the lefts inputs look like this
	// in R1CS:
	// [0, 1, 0, 0, 0, 0]
	// [0, 0, 0, 1, 0, 0]
	// [0, 1, 0, 0, 1, 0]
	// [5, 0, 0, 0, 0, 1]
	// so for example the second row of the transposed matrix represents all the
	// usage of the variable "x" as a left wire in the circuit -> [1,0,1,0]
	// Then we need to interpolate this as a polynomial such that:
	// p(0) = 1, p(1) = 0, p(2) = 1 and p(3) = 0
	// We return one polynomial for each variable
	t := m.Transpose()
	out := make([]Poly, 0, len(t))
	for _, variable := range t {
		pairs := make([]*pair, 0, len(variable))
		for i, v := range variable {
			pairs = append(pairs, &pair{I: i, V: v.ToFieldElement()})
		}
		// lagrange interpolation - uses an existing implementation
		sharePoly, err := share.RecoverPriPoly(Group, pairs, len(pairs), len(pairs))
		if err != nil {
			panic(err)
		}
		out = append(out, toPoly(sharePoly))
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
func (q QAP) Verify(sol Vector) bool {
	// create Z(x) = 1 so we can multiply it easily afterwards
	z := Poly([]Element{NewElement().One()})
	for i := 1; i <= len(q.left); i++ {
		// you multiply by (x - i) with i being as high as the number of gates,
		// which is exactly the number of coefficients of A
		minusi := NewElement().SetInt64(int64(i))
		minusi = minusi.Neg(minusi)
		xi := Poly([]Element{minusi, NewElement().One()})
		z = z.Mul(xi)
	}

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
	//   first gate etc
	// Left one
	left := Poly([]Element{NewElement().One()})
	for i, sval := range sol {
		left = left.Add(q.left[i].Mul(Poly([]Element{Value(sval).ToFieldElement()})))
	}

	right := Poly([]Element{NewElement().One()})
	for i, sval := range sol {
		right = right.Add(q.right[i].Mul(Poly([]Element{Value(sval).ToFieldElement()})))
	}

	out := Poly([]Element{NewElement().One()})
	for i, sval := range sol {
		out = out.Add(q.out[i].Mul(Poly([]Element{Value(sval).ToFieldElement()})))
	}

	// Now we are using these precedent polynomials in the satisfying equation
	// left(x) * right(x) - out(x) = h(x)z(x)
	eq := left.Mul(right).Sub(out)

	// now we try to verify if the equation is really satisfied by dividing eq /
	// z(x) : we should find no remainder
	_, rem := eq.Div(z)
	if rem.Degree() >= 0 {
		return false
	}
	return true
}
