package playsnark

import "github.com/drand/kyber/share"

// QAP is the polynomial version of a R1CS circuit
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
