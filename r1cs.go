package playsnark

// let's construct the r1cs matrix A_l, A_r A_o for the equation
// x^3 + x + 5 = 35

type Var struct {
	Index int
	Name  string
}

func newVar(i int, name string) Var {
	return Var{
		Index: i, Name: name,
	}
}

// createVariables returns all the variables for the R1CS
// circuit we want to prove with their associated index
// u = x * x
// v = u * x
// w = v + x
// out = w + 5
type Variables []Var

// createVariables returns the indexes of the variables and the number of input
// output variables
func createVariables() (Variables, int) {
	var ordering []Var
	for i, name := range []string{"const", "x", "out", "u", "v", "w"} {
		ordering = append(ordering, newVar(i, name))
	}
	return ordering, 3
}

func (v *Variables) IndexOf(name string) int {
	for _, va := range *v {
		if va.Name == name {
			return va.Index
		}
	}
	panic("aie")
}

func (v *Variables) ConstraintOn(names ...string) Vector {
	var constraints Vector
	for _, variable := range *v {
		var found bool
		for _, name := range names {
			if variable.Name == name {
				found = true
				break
			}
		}
		if found {
			constraints = append(constraints, 1)
		} else {
			constraints = append(constraints, 0)
		}
	}
	return constraints
}

// createWitness returns the witness vector that satisfies the R1CS equation -
// see tests for more info on the equation.
// The vector is created by assigning the value of each variable in the circuit
// to the value that satisfies the whole circuit.
// Remember the flattened circuit
// u = x * x
// v = u * x
// w = v + x
// out = w + 5

// variable "const" has 1, always
// variable x is 3 because 3 is a solution to the equation  x^3 + x + 5 = 35
// variable "out" is 35 because that's the right side of the equation
// variable "u" is 9 because it's 3x3
// variable "v" is 27 because it's ux3
// variable "w" is 30 because it v + x = 27 + 3
func createWitness() Vector {
	vars, _ := createVariables()
	nbVars := len(vars)
	solution := make(Vector, nbVars)
	solution[vars.IndexOf("const")] = 1
	solution[vars.IndexOf("x")] = 3
	solution[vars.IndexOf("out")] = 35
	solution[vars.IndexOf("u")] = 9
	solution[vars.IndexOf("v")] = 27
	solution[vars.IndexOf("w")] = 30
	return solution
}

// R1CS describe the circuit: it contains the 3 matrixes left right and output
// such that when given a valid vector of solution "s", the equation
// <left,s> * <right,s> - <out,s> = 0 is satisfied.
type R1CS struct {
	vars  Variables
	nbIO  int
	left  Matrix
	right Matrix
	out   Matrix
}

// createR1CS returns the R1CS for the problem we consider
// x^3 + x + 5 = 35 using the variables described as in "createVariables"
func createR1CS() R1CS {
	variables, _ := createVariables()
	// We create a list of vector, i.e. a matrix
	// for all the left, right and out output
	var left, right, out Matrix

	// u = x * x
	left = append(left, variables.ConstraintOn("x"))
	right = append(right, variables.ConstraintOn("x"))
	out = append(out, variables.ConstraintOn("u"))

	// v = u * x
	left = append(left, variables.ConstraintOn("u"))
	right = append(right, variables.ConstraintOn("x"))
	out = append(out, variables.ConstraintOn("v"))

	// w = v + x
	// here since it is an addition, we simply mark the two in the left
	// variable for example that will get added togeter during the dot product
	// and only mark as one the constant in the right input so it gives
	// w = (v + x) * 1
	left = append(left, variables.ConstraintOn("v", "x"))
	right = append(right, variables.ConstraintOn("const"))
	out = append(out, variables.ConstraintOn("w"))

	// out = w + 5
	// Here it's same as before, so we consider
	// out = (w + 5) * 1
	l := variables.ConstraintOn("const", "w")
	// here I'm cheating because I know the first value is the constant one
	// TODO make API to change constant value in vector
	l[0] = l[0] * 5
	left = append(left, l)
	right = append(right, variables.ConstraintOn("const"))
	out = append(out, variables.ConstraintOn("out"))
	return R1CS{
		vars:  variables,
		left:  left,
		right: right,
		out:   out,
		// "const", "x", and "out"
		nbIO: 3,
	}
}
