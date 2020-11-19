package pedersen

import (
	"bytes"
	"crypto/cipher"
	"crypto/sha256"
	"encoding/binary"

	"github.com/drand/kyber"
	"github.com/drand/kyber/util/random"
	"github.com/nikkolasg/playsnark/pedersen/commit"
)

// Group has points on it and can create scalar from the scalar fields
type Group = kyber.Group

// Scalar of the field of the curve
type Scalar = kyber.Scalar

// Point in the group (in our case it's elliptic curve so it's a point)
type Point = kyber.Point

type Setup struct {
	commit.Setup
	h      Point
	degree int
}

// NewSetup returns the information necessary for using the pedersen polynomail
// commitment
func NewSetup(g Group, degree int) Setup {
	s := commit.NewSetup(g, degree+1)
	h := g.Point().Pick(oracle(g, degree))
	return Setup{
		Setup:  s,
		h:      h,
		degree: degree,
	}
}

type Polynomial = []Scalar

type Commitment = Point

// Commit the pedersen commitment of the coefficients
func Commit(s Setup, p Polynomial, r Scalar) Commitment {
	return commit.Commit(s.Setup, p, r)
}

func Open(s Setup, p Polynomial, c Commitment, z, r Scalar) Proof {
	// v = p(z)
	v := evalPoly(s.g, p, eval)
	// random poly p2 such that  p2(z) = 0
	// naive strategy: create a random one, then eval p2(z) = a then substract
	// a from the free coefficient to p2' such  that p2'(z) = p2(z) - a = 0
	p2 := randomPoly(s.g, s.degree)
	a := evalPoly(s.g, p2, z)
	p2[0] = p2[0].Sub(p2[0], a)
	// sample randomness w
	w := s.g.Scalar().Pick(random.New())
	// commit to p2 using w
	c2 := commit.Commit(s.Setup, p2, w)
	// compute challenge from random oracle
	alpha := oracleFromScalars(c, z, v, c2)

}

func evalPoly(g Group, p Polynomial, x Scalar) Scalar {
	xi := x.Clone()
	v := g.Scalar().Zero()
	for j := len(p) - 1; j >= 0; j-- {
		v.Mul(v, xi)
		v.Add(v, p[j])
	}
	return v
}

func randomPoly(g Group, degree int) Polynomial {
	p := make(Polynomial, degree+1)
	for i := 0; i <= degree; i++ {
		p[i] = g.Scalar().Pick(random.New())
	}
	return p
}

// returns an oracle for deriving a point for pedersen polynomial commitment
func oracle(g Group, degree int) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-polycommit"))
	h.Write([]byte(g.String()))
	binary.Write(h, binary.LittleEndian, degree)
	var b bytes.Buffer
	b.Write(h.Sum(nil))
	return random.New(&b)
}

func oracleFromScalar(g Group, scalars ...Scalar) cipher.Stream {
	var h = sha256.New()
	h.Write([]byte("pedersen-polycommit-scalar"))
	h.Write([]byte(g.String()))
	for _, s := range scalars {
		s.WriteTo(h)
	}
	var b bytes.Buffer
	b.Write(h.Sum(nil))
	return random.New(&b)
}
