package commit

import (
	"crypto/cipher"
	"crypto/sha256"
	"encoding/binary"

	"github.com/drand/kyber"
	"github.com/drand/kyber/util/random"
	"golang.org/x/crypto/blake2s"
)

// Group has points on it and can create scalar from the scalar fields
type Group = kyber.Group

// Scalar of the field of the curve
type Scalar = kyber.Scalar

// Point in the group (in our case it's elliptic curve so it's a point)
type Point = kyber.Point

// Setup contains the setup information for pedersen commitments
type Setup struct {
	Gs []Point
	S  Point
	G  Group
}

// NewSetup returns the setup necessary to compute a commitment
func NewSetup(g Group, l int) Setup {
	gs := make([]Point, l)
	for i := 0; i < l; i++ {
		gs[i] = g.Point().Pick(oracle(g, l, i))
	}
	s := g.Point().Pick(oracle(g, l, l+1))
	return Setup{
		Gs: gs,
		S:  s,
		G:  g,
	}
}

// returns an oracle for deriving a point for pedersen commitment
func oracle(g Group, l, pos int) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-commit"))
	h.Write([]byte(g.String()))
	binary.Write(h, binary.LittleEndian, int32(l))
	binary.Write(h, binary.LittleEndian, int32(pos))
	xof, err := blake2s.NewXOF(0, h.Sum(nil))
	if err != nil {
		panic(err)
	}
	return random.New(xof)
}

// Commit compute the commitment of the points given the setup and the
// randomness to use. If r == nil, then r is set 0, so the commitment is not
// hiding anymore.
func Commit(s Setup, msgs []Scalar, r Scalar) Point {
	if r == nil {
		r = s.G.Scalar().Zero()
	}
	var acc = s.G.Point().Mul(r, s.S)
	for i, m := range msgs {
		gxi := s.G.Point().Mul(m, s.Gs[i])
		acc.Add(acc, gxi)
	}
	return acc
}
