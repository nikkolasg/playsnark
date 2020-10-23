[![Build Status](https://travis-ci.com/nikkolasg/playsnark.svg?branch=master)](https://travis-ci.com/nikkolasg/playsnark)
[![codecov](https://codecov.io/gh/nikkolasg/playsnark/branch/master/graph/badge.svg?token=F29VQZ22KO)](undefined)

**WARNING: This is a toy implementation with many pitfalls and bugs and is NOT
intented to be securely deployed, do ABSOLUTELY NOT use it in any production
setting, you've been warned!** For example, this implementation keeps around 
very sensitive information during trusted setup or proof generation. The goal is
to learn how the different proofs system work.

# Playsnark: a playground to learn proofs systems

The goal of this code is to learn different proofs systems from the
ground up. The code is a _very minimal_ implementation of the [Pinocchio proof
system](https://eprint.iacr.org/2013/879.pdf), over a toy function. 
The code is tested at each of these levels and to highlight what is the
satisfying "equation" that the prover is trying to prove.

## R1CS

The code in `r1cs.go` shows how to translate any function into R1CS format including the 
different variables and gates.  
For example, to construct a circuit for the equation `x^3 + x + 5 = 35`, you
need to decompose the operation with addition and multiplication gates.
```go
c := NewR1CS()
c.NewInput("x")
c.NewOutput("out")
// u = x * x
c.NewVar("u")
// v = u * x
c.NewVar("v")
// w = v + x
c.NewVar("w")

// u = x * x
c.Mul("x", "x", "u")
// v = u * x
c.Mul("u", "x", "v")
// w = v + x
c.Add("v", "x", "w")
// out = w + 5
c.AddConst("w", 5, "out")
```

Then you can construct the solution vector, giving one value to each variable,
as the following, where `r` is the R1CS struct. The intermediate values are the
one constructed when first giving the value 3 to the variable x.
**Note**: in R1CS, we have the first variable as being the variable "const"
which just takes the value 1, this is useful when doing addition. Then we have
the input variables ("x") and then the output variables ("out"), then the
intermediate variables ("v","w","y").
```go
var solution = make(Vector, len(r.vars))
solution[r.vars.IndexOf("const")] = 1
solution[r.vars.IndexOf("x")] = 3
solution[r.vars.IndexOf("out")] = 35
// v = x * x = 3 * 3
solution[r.vars.IndexOf("u")] = 9
// v = u * x = 9 * 3
solution[r.vars.IndexOf("v")] = 27
// w = v + x = 27 + 3
solution[r.vars.IndexOf("w")] = 30
```

## Quadratic Arithmetic Programs (QAP)

Now that we got the R1CS part, we need to translate it to equations involving
polynomials. To do so, we create one polynomial for each variable, which
evaluates to the value of the variable for each gate, with one polynomial for
the left input, the right input and the output of each gate. For example the fourth
polynomial (representing "u") for the output of the gate 1 evaluates to 1
because "u = x * x", so u is the output.

You can create the QAP from the R1CS as following:
```go
r1cs := createR1CS()
qap := ToQAP(r1cs)
```
and verify if the QAP is formed correctly by giving a witness to the problem and
see if the QAP equation resolves (see the `qap.go` for more info):
```go
s := createWitness(r1cs)
require.True(t, qap.IsValid(s))
```

## SNARK part

### Trusted Setup

Currently, playsnark implements the Pinocchio proof system which requires a
trusted setup for the prover and verifier. It enables the prover to blindly
evaluates the polynomials of the QAP above, where blindly means it evaluates
them at a unknown point s.

You can simply generate a trusted setup for the current circuit like so:
```go
r1cs := createR1CS()
qap := ToQAP(r1cs)
setup := NewTrustedSetup(qap)
```

The trusted setup has two parts, the **evaluation key** required by the prover
to evaluate its circuit with its witness and produce a valid proof, and the
**verification key** that is required from the verifier to verify a proof.
Note that this implementation keeps other materials that we call "toxic waste"
that should not be kept around after the trusted setup, but for the sake of
comprehension and testing, we keep it. For example, using the toxic waste, you
could generate a malicious proof passing the verification routine but does not
guarantee that the prover knows a witness to the circuit.
```go
type TrustedSetup struct {
	EK EvalKey         // required by the prover
	VK VerificationKey // required by the verifier
	t  ToxicWaste      // to be deleted - left for testing
}
```

### Generating Proof

The prover then needs the evaluation key, the QAP polynomials and the solution
vector:
```go
r1cs := createR1CS()
s := createWitness(r1cs)
qap := ToQAP(r1cs)
setup := NewTrustedSetup(qap)
proof := GenProof(setup.EK, qap, s)

```

### Verifiying proof

The verifier then runs three basic check with the proof, the verification key of
the trusted setup and the QAP, which are all explained in depth in
`pinocchio.go:VerifyProof()`:
* Division check: looks if the QAP equation resolves "in the exponent"
* Correct polynomial check: We check that the prover correctly used the
  polynomials of the QAP, the ones that were blindly evaluated in CRS
  _individually_.
* Linear combination check: We check that the prover used the same value for the
  variables for the three polynomials 

You can do so as the following: 
```go
r1cs := createR1CS()
s := createWitness(r1cs)
qap := ToQAP(r1cs)
diff := qap.nbVars - qap.nbIO
setup := NewTrustedSetup(qap)
proof := GenProof(setup.EK, qap, s)
fmt.Println(VerifyProof(setup.VK, qap, proof, s[:diff]))
```

Note the `diff` variable is just because in this test we have access to
everything. In reality, the verifier only needs the values of the input /
outputs so we restrict the solution vector to these variable when giving it to
the verifier.


## Resources

Well the first one I used is the series of [Vitalik blog post](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649), then I looked at this more technical small [paper](https://chriseth.github.io/notes/articles/zksnarks/zksnarks.pdf) and finally to implement correctly the Pinocchio proof system I used the original [paper](https://eprint.iacr.org/2013/879.pdf) as well as the [paper](https://eprint.iacr.org/2013/879.pdf) derived after that succintly describes the algorithm using an asymmetric pairing from Ben-Sasson, Chiesa, Tromer and Virza.
