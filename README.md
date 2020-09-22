[![Build Status](https://travis-ci.com/nikkolasg/playsnark.svg?branch=master)](https://travis-ci.com/nikkolasg/playsnark)
[![codecov](https://codecov.io/gh/nikkolasg/playsnark/branch/master/graph/badge.svg?token=F29VQZ22KO)](undefined)

**WARNING: This is a toy implementation with potentially many pitfalls and bugs, do ABSOLUTELY NOT use it in any production setting, you've been warned!**

# Playsnark: a playground to learn proofs systems

The goal of this code is to learn different proofs systems from the
ground up. The code is a _very minimal_ implementation of the [Pinocchio proof
system](https://eprint.iacr.org/2013/879.pdf), over a toy function. 
The code is tested at each of these levels and to highlight what is the
satisfying "equation" that the prover is trying to prove.

Currently, the code allows a prover to prove he knows an x such that `f(x) = x^3 + x + 5 = 35`.

You can follow the code by first looking at:
* `r1cs.go`: how does the function gets translated into R1CS format including the
  different variables and gates
* `qap.go`: how does the r1cs format gets translated to polynomials and what the
  check of validity becomes.
* `pinocchio.go`: what is the trusted setup in Pinocchio, what the prover is
  doing and what the verifier must compute in order to make sure the prover
  didn't cheat.

For each of these files, you can look at the corresponding tests that check whether the transformation is correct or not.

There are two additional files as well:
* `curve.go` contains the code defining the curve we use, the points and the generators
* `algebra.go` contains some logic related to matrices, vectors, polynomials needed in the r1cs and QAP modules.

## Future

I started this project merely as a self learning tool but I might continue to develop other proof systems in the future as time permits.

## Resources

Well the first one I used is the series of [Vitalik blog post](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649), then I looked at this more technical small [paper](https://chriseth.github.io/notes/articles/zksnarks/zksnarks.pdf) and finally to implement correctly the Pinocchio proof system I used the original [paper](https://eprint.iacr.org/2013/879.pdf) as well as the [paper](https://eprint.iacr.org/2013/879.pdf) derived after that succintly describes the algorithm using an asymmetric pairing from Ben-Sasson, Chiesa, Tromer and Virza.
