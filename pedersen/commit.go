package pedersen

import (
	"github.com/drand/kyber"
	"github.com/drand/kyber/util/random"
)

// Group has points on it and can create scalar from the scalar fields
type Group = kyber.Group

// Scalar of the field of the curve
type Scalar = kyber.Scalar

// Point in the group (in our case it's elliptic curve so it's a point)
type Point = kyber.Point

// Setup contains the setup information for pedersen commitments
type Setup struct {
	gs []Point
	h  Point
	g  Group
}

// NewSetup returns the setup necessary to compute a commitment
func NewSetup(g Group, l int) Setup {
	gs := make([]Point, l)
	for i := 0; i < l; i++ {
		gs[i] = g.Point().Pick(random.New())
	}
	h := g.Point().Pick(random.New())
	return Setup{
		gs: gs,
		h:  h,
		g:  g,
	}
}

// Commit compute the commitment of the points given the setup and the
// randomness to use.
func Commit(s Setup, msgs []Scalar, r Scalar) Point {
	var acc = s.g.Point().Mul(r, s.h)
	for i, m := range msgs {
		gxi := s.g.Point().Mul(m, s.gs[i])
		acc.Add(acc, gxi)
	}
	return acc
}
