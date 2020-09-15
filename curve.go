package playsnark

import (
	"github.com/drand/kyber"
	bls "github.com/drand/kyber-bls12381"
)

type Element = kyber.Scalar
type Commit = kyber.Point
type G1 = kyber.Point
type G2 = kyber.Point

var Suite = bls.NewBLS12381Suite()
var Group = Suite.G1()
var G2Group = Suite.G2()

func (v Value) ToFieldElement() Element {
	return Group.Scalar().SetInt64(int64(v))
}

func NewElement() Element {
	return Group.Scalar().Zero()
}

func NewG1() G1 {
	return Group.Point().Base()
}

func NewG2() G2 {
	return G2Group.Point().Base()
}

type Target = kyber.Point

// Pair returns e(a,b)
func Pair(a G1, b G2) Target {
	return Suite.Pair(a, b)
}

var zero = NewElement()
var one = NewElement().SetInt64(1)

var zeroG1 = NewG1().Null()
var zeroG2 = NewG2().Null()
var zeroGT = Suite.GT().Point().Null()
