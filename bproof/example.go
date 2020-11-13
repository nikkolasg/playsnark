package bproof

// GenerateCircuit generates a sample circuit for the computation
// x^2 + y^2 + 5 = z <=>
// x^2 + y^2 + 5 - z = 0
func GenerateCircuit() *ProverCircuit {
	c := NewProverCircuit()
	// 2^2 + 3^2 + 5 = 18 = z
	x := NewScalar().SetInt64(2)
	y := NewScalar().SetInt64(3)
	_, _, x2 := c.Allocate(x, x)
	_, _, y2 := c.Allocate(y, y)
	x2y2 := c.Add(x2, y2)
	left := c.AddCst(x2y2, NewScalar().SetInt64(5))

	// Solution here is to use the "left" Linear combination as a variable and
	// constraint on the following
	//        left     * right = z
	// (x^2 + y^2 + 5) * 1     = z
	// For that, it's a hack but I need declare a "one" variable such that
	// 1*1=1, it's then a var i can re-use everywhere in my circuit
	// I could change the circuit to have a "one" variable though.
	one, _, _ := c.Allocate(NewScalar().SetInt64(1), NewScalar().SetInt64(1))
	_, _, _ = c.Mul(left, one)
	return c
}
