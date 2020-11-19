package commit

import (
	"bytes"
	"crypto/cipher"
	"crypto/sha256"
	"encoding/binary"

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
	s  Point
	g  Group
}

// NewSetup returns the setup necessary to compute a commitment
func NewSetup(g Group, l int) Setup {
	gs := make([]Point, l)
	for i := 0; i < l; i++ {
		gs[i] = g.Point().Pick(oracle(g, l, i))
	}
	s := g.Point().Pick(oracle(g, l, l+1))
	return Setup{
		gs: gs,
		s:  s,
		g:  g,
	}
}

// returns an oracle for deriving a point for pedersen commitment
func oracle(g Group, l, pos int) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-commit"))
	h.Write([]byte(g.String()))
	binary.Write(h, binary.LittleEndian, l)
	binary.Write(h, binary.LittleEndian, pos)
	var b bytes.Buffer
	b.Write(h.Sum(nil))
	return random.New(&b)
}

// Commit compute the commitment of the points given the setup and the
// randomness to use.
func Commit(s Setup, msgs []Scalar, r Scalar) Point {
	var acc = s.g.Point().Mul(r, s.s)
	for i, m := range msgs {
		gxi := s.g.Point().Mul(m, s.gs[i])
		acc.Add(acc, gxi)
	}
	return acc
}
