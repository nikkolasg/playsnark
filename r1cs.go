package playsnark

// let's construct the r1cs matrix A_l, A_r A_o for the equation
// x^3 + x + 5 = 35

// Var holds the name of the variable and its index in the vector of variables
type Var struct {
	Index int
	Name  string
}

func newVar(i int, name string) Var {
	return Var{
		Index: i, Name: name,
	}
}

type Variables []Var

// IndexOf returns the index of a variable
func (v *Variables) IndexOf(name string) int {
	for _, va := range *v {
		if va.Name == name {
			return va.Index
		}
	}
	panic("plouf")
}

// ConstraintOn returns a vector where the i-th entry is set to 1 if there is a
// variable given whose name has the i-th index.
// For example, `ConstraintOn("x","out")` will return [0,1,1,0,0,0]
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
func createWitness(r R1CS) Vector {
	var solution = make(Vector, len(r.vars))
	solution[r.vars.IndexOf("const")] = 1
	solution[r.vars.IndexOf("x")] = 3
	solution[r.vars.IndexOf("out")] = 35
	solution[r.vars.IndexOf("u")] = 9
	solution[r.vars.IndexOf("v")] = 27
	solution[r.vars.IndexOf("w")] = 30
	return solution
}

// R1CS describe the circuit: it contains the 3 matrixes left right and output
// such that when given a valid vector of solution "s", the equation
// <left,s> * <right,s> - <out,s> = 0 is satisfied.
type R1CS struct {
	// slice of name of the input variables
	inputs []string
	// slice of name of the outputs variables
	outputs []string
	// slice of name of the intermediate variable in the circuit - it is
	// necessary to separate those as the verifier knows already the input and
	// outputs values, but doesn't know the intermediate ones, because the
	// prover is doing the computation, not the verifier. The prover then shows
	// correctness of the computations using these variables.
	intermediates []string
	// the concatatenated variables [const, inputs..., ouputs..., intermediates]
	// in order
	vars Variables
	// the matrices representing the wire of each variable for each gate:
	// the rows of the matrices represent the gate and the columns the variable
	// so if left(row=2,column=3) = 1, then it means the 3rd variable in `vars`
	// is wired up as the left input to the 2nd gate constraint.
	left  Matrix
	right Matrix
	out   Matrix
}

func NewR1CS() R1CS {
	return R1CS{}
}

func (r *R1CS) nbIO() int {
	return 1 + len(r.inputs) + len(r.outputs)
}

func (r *R1CS) NewInput(name string) {
	r.inputs = append(r.inputs, name)
	r.mergeVars()
}

func (r *R1CS) NewOutput(name string) {
	r.outputs = append(r.outputs, name)
	r.mergeVars()
}

func (r *R1CS) NewVar(name string) {
	r.intermediates = append(r.intermediates, name)
	r.mergeVars()
}

func (r *R1CS) mergeVars() {
	var vars = []Var{newVar(0, "const")}
	for _, n := range r.inputs {
		vars = append(vars, newVar(len(vars), n))
	}
	for _, n := range r.outputs {
		vars = append(vars, newVar(len(vars), n))
	}
	for _, n := range r.intermediates {
		vars = append(vars, newVar(len(vars), n))
	}
	r.vars = vars
}

// Mul takes the name of the left variable, right variable and the output
// variable and wires them such that left * right = out
func (r *R1CS) Mul(left, right, out string) {
	r.left = append(r.left, r.vars.ConstraintOn(left))
	r.right = append(r.right, r.vars.ConstraintOn(right))
	r.out = append(r.out, r.vars.ConstraintOn(out))
}

// Add takes the name of the first variable, second variable and the output
// variable and wires them such that var1 + var2 = out
func (r *R1CS) Add(var1, var2, out string) {
	// here since it is an addition, we simply mark the two in the left
	// variable for example that will get added togeter during the dot product
	// and only mark as one the constant in the right input so it gives
	// w = (v + x) * 1
	r.left = append(r.left, r.vars.ConstraintOn(var1, var2))
	r.right = append(r.right, r.vars.ConstraintOn("const"))
	r.out = append(r.out, r.vars.ConstraintOn(out))
}

// AddConst takes the name of the first variable, the value of the constant to
// add and the output variable name and wires them such that var1 + const = out
func (r *R1CS) AddConst(var1 string, add int, out string) {
	row := r.vars.ConstraintOn("const", var1)
	row[0] = row[0] * Value(add)
	r.left = append(r.left, row)
	r.right = append(r.right, r.vars.ConstraintOn("const"))
	r.out = append(r.out, r.vars.ConstraintOn(out))
}

// createR1CS returns the R1CS for the problem we consider
// x^3 + x + 5 = 35 using the variables described as in "createVariables"
func createR1CS() R1CS {
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
	return c
}
