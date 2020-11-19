package playsnark

// Implements the Inner Product Argument from bulletproof
// paper https://eprint.iacr.org/2017/1066.pdf
//
// Goal is to prove, from input g,h \in G^n, scalar c, and P and G \in G that
// prover knows a and b \in Z^n such that P = g^a * h^b and c = < a b >
