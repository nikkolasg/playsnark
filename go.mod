module github.com/nikkolasg/playsnark

go 1.15

require (
	github.com/drand/kyber v1.1.4
	github.com/drand/kyber-bls12381 v0.2.1
	github.com/kilic/bls12-381 v0.0.0-20200820230200-6b2c19996391 // indirect
	github.com/stretchr/testify v1.4.0
	go.dedis.ch/kyber/v3 v3.0.9
	golang.org/x/crypto v0.0.0-20200820211705-5c72a883971a
)

//replace github.com/drand/kyber => ../drand/kyber
