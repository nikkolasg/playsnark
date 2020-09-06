# Snark tuto

Let's try to make a full snark to prove the following equation:
x^3 + x + 5 = 35

## R1CS

### Flattening

We need to flatten out this equation so we express it in a series of gate
(multiplication or addition):

u = x * x
v = u * x
w = v + x
out = w + 5

### Vector of solutions

Then we express these list of variables in a vector of solutions:
s = (const, x, out, u, v, w)

A few observations, in R1CS:
* we express "const" a variable that contains the value 1.
* we first express the input variable(s) then the output variable(s) then all
  intermediate variables in order afterwards

### R1CS equation

In RC1S, we must define three vectors A_l, A_r, A_o such that 
<A_l, s> * <A_r, s> - <A_o, s> = 0

Each vector A is of the same length as f obviously. 
Each vector A is filled such that the equation matches.

We will define such triplets (A_l,A_r,A_o) for each gates:

* Gate u = x * x
    + Remind the order of the solution vector: s = (const, x, out, u, v, w)
    + A_l = (0, 1, 0, 0, 0, 0)
    + A_r = (0, 1, 0, 0, 0, 0)
    + A_o = (0, 0, 0, 1, 0, 0)
    + We check: <A_l, s> * <A_r, s> - <A_o, s> = 0
        + indeed, we have (1 * x) * (1 * x) - (1 * u) = 0 <=> u = x * x
v = u * x
w = v + x
out = w + 5


