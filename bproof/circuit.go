package bproof

import (
	"bytes"
	"fmt"
)

// In Bulletproof style of proof, we need to represent the circuit as:
// 1. three vectors of witnesses a,b,c  such that a_i * b_i = c_i
// 2. A _list_ of three vectors and one element such that for each such tuple we
// have SUM(i:0 -> N) w_l * a + SUM w_r * b + SUM w_o * c = k
// where N is the number of elements in a (or b or c). There are Q such
// linear constraints.
// implementation of circuit heavily inspired by the dalek bulletproof impl. but
// maybe simplified and less efficient ;)

// LinearCombination is a linear combination of different variables represented
// by a list of products. In this system, we can Add different linear
// combination together
type LinearCombination interface {
	Products() []product
}

// product represents the multiplication of a variable by a coefficient. It is a
// term of a linear combination.
type product struct {
	Witness variable
	Coeff   Scalar
}

func newOneProduct(v variable) product {
	return product{Witness: v, Coeff: NewScalar().One()}
}
func newMinusProduct(v variable) product {
	a := NewScalar().One()
	return product{Witness: v, Coeff: a.Neg(a)}
}

// kind represents the type of variable that is possible: Left Right OUT or
// Constant (i.e. a "k" coefficient in a linear combination).
type kind int

const (
	VLEFT kind = iota
	VRIGHT
	VOUT
	CST
)

type variable struct {
	Index int
	Kind  kind
}

func (v *variable) Products() []product {
	return []product{newOneProduct(*v)}
}

type linearCombination []product

func (l *linearCombination) Products() []product {
	return *l
}

type Constraint []Scalar

type Circuit struct {
	nbVars      int
	constraints []LinearCombination

	// filled later by "Flatten"
	wl []Constraint
	wr []Constraint
	wo []Constraint
	wk Constraint
}

type ProverCircuit struct {
	Circuit
	a []Scalar
	b []Scalar
	c []Scalar
}

func NewCircuit() *Circuit {
	return &Circuit{}
}

func NewProverCircuit() *ProverCircuit {
	return &ProverCircuit{Circuit: *NewCircuit()}
}

// Mul creates three new variables in the circuit such that a * b =
// variable in this constraint systems, and then we create a linear
// combination to say that the "previous" variable equals the "new
// variable"- this operation is called "wiring"
func (c *Circuit) Mul(a, b LinearCombination) (left variable, right variable, out variable) {
	left, right, out = c.newMultiplier()
	c.wire(a, left)
	c.wire(b, right)
	return left, right, out
}

func (c *Circuit) Add(a, b LinearCombination) LinearCombination {
	lc := linearCombination(append(a.Products(), b.Products()...))
	return &lc
}

func (c *Circuit) MulCst(a LinearCombination, cst Scalar) LinearCombination {
	for _, p := range a.Products() {
		p.Coeff = p.Coeff.Mul(p.Coeff, cst)
	}
	return a
}

func (c *Circuit) AddCst(a LinearCombination, cst Scalar) LinearCombination {
	// No need to specify index since the constant in the linear constraint is
	// not a vector it's a single value
	p := product{Coeff: cst, Witness: variable{Kind: CST}}
	lc := linearCombination(append(a.Products(), p))
	return &lc
}

// NewMultiplier creates a left and output variable to use for multiplication
// and additions
func (c *Circuit) newMultiplier() (left variable, right variable, out variable) {
	left = variable{Index: c.nbVars, Kind: VLEFT}
	right = variable{Index: c.nbVars, Kind: VRIGHT}
	out = variable{Index: c.nbVars, Kind: VOUT}
	c.nbVars += 1
	return left, right, out
}

// Constraint saves this constraints such that when this linear combination is
// evaluated with the witness vectors a b c then the equation is true.
func (c *Circuit) Constraint(lc LinearCombination) {
	c.constraints = append(c.constraints, lc)
}

// wire takes a linear combination and sets the variable v as equal to the
// result of that linear combination. In practice it means
// LC + (-1) * v = 0
func (c *Circuit) wire(input LinearCombination, v variable) {
	lc := linearCombination(append(input.Products(), newMinusProduct(v)))
	fmt.Println("New Wire Added: ", lc)
	c.Constraint(&lc)
}

// Allocate creates a left right and output variable where left is set to a,
// right is set to b, and out is set to a * b.
// At the beginnin of a circuit, one can create a "one" variable .
func (c *ProverCircuit) Allocate(a, b Scalar) (left LinearCombination, right LinearCombination, out LinearCombination) {
	o := NewScalar().Mul(a, b)
	c.pushWitness(a, b, o)
	av, bv, outv := c.newMultiplier()
	return &av, &bv, &outv
}

// Eval takes a linear combination and evaluates its value with the witness
// vector. The result can then be used an input in sbusequent gates.
func (c *ProverCircuit) eval(a LinearCombination) Scalar {
	s := zero()
	for _, p := range a.Products() {
		switch p.Witness.Kind {
		case VLEFT:
			s = s.Add(s, NewScalar().Mul(c.a[p.Witness.Index], p.Coeff))
		case VRIGHT:
			s = s.Add(s, NewScalar().Mul(c.b[p.Witness.Index], p.Coeff))
		case VOUT:
			s = s.Add(s, NewScalar().Mul(c.c[p.Witness.Index], p.Coeff))
		case CST:
			s = s.Add(s, p.Coeff)
		}
	}
	return s
}

func (c *ProverCircuit) Mul(a, b LinearCombination) (l variable, r variable, o variable) {
	va := c.eval(a)
	vb := c.eval(b)
	out := NewScalar().Mul(va, vb)
	fmt.Println("Prover.Mul(): a", va, "  - b", vb, " --> ", out)
	c.pushWitness(va, vb, out)
	l, r, o = c.Circuit.Mul(a, b)
	return
}
func (c *ProverCircuit) pushWitness(a, b, out Scalar) {
	c.a = append(c.a, a)
	c.b = append(c.b, b)
	c.c = append(c.c, out)
}

// Flatten will compute the fully fledged arithemtic circuit representation
// w_l w_r w_o are matrices of size N x Q and w_k \in F^Q where N is the number
// of variable and Q is the number of linear constraints.
// for constraint K
// SUM a_i * w_l_K_i + SUM b_i * w_r_i + SUM c_i * w_o_i = w_k
// Note that this implementation is not optimized at all: it generates full
// length vector while they are supposedly very sparse.
func (c *Circuit) Flatten() {
	n := c.nbVars
	q := len(c.constraints)

	wl := make([]Constraint, q)
	wr := make([]Constraint, q)
	wo := make([]Constraint, q)
	wk := make(Constraint, q)
	for i := 0; i < q; i++ {
		wl[i] = make(Constraint, n)
		wr[i] = make(Constraint, n)
		wo[i] = make(Constraint, n)
		for j := 0; j < n; j++ {
			wl[i][j] = zero()
			wr[i][j] = zero()
			wo[i][j] = zero()
		}
		wk[i] = zero()

		for _, product := range c.constraints[i].Products() {
			v := product.Witness
			switch product.Witness.Kind {
			case VLEFT:
				wl[i][v.Index] = product.Coeff
			case VRIGHT:
				wr[i][v.Index] = product.Coeff
			case VOUT:
				wo[i][v.Index] = product.Coeff
			case CST:
				// wl + wr + wo = wk
				// When this circuit adds constant it just adds them to wk, so
				// to satisfy the constraint we negate wk
				wk[i] = NewScalar().Neg(product.Coeff)
			}
		}
	}
	c.wl = wl
	c.wr = wr
	c.wo = wo
	c.wk = wk
}

func (c *Circuit) String() string {
	var b bytes.Buffer
	for k := 0; k < len(c.wl); k++ {
		b.WriteString(fmt.Sprintf("Constraint %d: \n", k))
		b.WriteString(fmt.Sprintf("\t- w_l: %v\n", c.wl[k]))
		b.WriteString(fmt.Sprintf("\t- w_r: %v\n", c.wr[k]))
		b.WriteString(fmt.Sprintf("\t- w_o: %v\n", c.wo[k]))
		b.WriteString(fmt.Sprintf("\t- w_k: %v\n", c.wk[k]))
	}
	return b.String()
}

func (c *ProverCircuit) String() string {
	var b bytes.Buffer
	b.WriteString(c.Circuit.String())
	b.WriteString(fmt.Sprintf("Witness A: %v\n", c.a))
	b.WriteString(fmt.Sprintf("Witness B: %v\n", c.b))
	b.WriteString(fmt.Sprintf("Witness C: %v\n", c.c))
	return b.String()
}
