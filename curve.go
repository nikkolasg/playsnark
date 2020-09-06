package playsnark

import (
	"github.com/drand/kyber"
	bls "github.com/drand/kyber-bls12381"
)

type Element = kyber.Scalar
type Commit = kyber.Point

var Group = bls.NewBLS12381Suite().G1()

func (v Value) ToFieldElement() Element {
	return Group.Scalar().SetInt64(int64(v))
}

func NewElement() Element {
	return Group.Scalar().Zero()
}
