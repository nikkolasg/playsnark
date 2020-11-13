package bproof

import (
	"github.com/drand/kyber"
	"github.com/drand/kyber/group/edwards25519"
)

type Scalar = kyber.Scalar

var Curve = edwards25519.NewBlakeSHA256Ed25519()

func NewScalar() Scalar {
	return Curve.Scalar()
}

func zero() Scalar {
	return NewScalar().Zero()
}
