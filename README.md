[![Build Status](https://travis-ci.com/nikkolasg/playsnark.svg?branch=master)](https://travis-ci.com/nikkolasg/playsnark)

# Playsnark: a playground to learn proofs systems

The goal of this code is simply to learn different proofs systems from the
ground up. It includes how to interpret any function into circuits, then into
polynomials and how to prove that you know a witness. 
The code is tested at each of these levels and to highlight what is the
satisfying "equation" that the prover is trying to prove.

Currently, this code implements the [Pinocchio proof
system](https://eprint.iacr.org/2013/879.pdf), over a toy example.  Namely, it
allows a prover to prove he knows an x such that `f(x) = x^3 + x + 5 = 35`.
You can follow the code by first looking at:
* `r1cs.go`: how does the function gets translated into R1CS format including the
  different variables and gates
* `qap.go`: how does the r1cs format gets translated to polynomials and what the
  check of validity becomes.
* `pinocchio.go`: what is the trusted setup in Pinocchio, what the prover is
  doing and what the verifier must compute in order to make sure the prover
  didn't cheat.


