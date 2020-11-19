package playsnark

type Value int

type Vector []Value

type Matrix []Vector

func NewMatrix(rows []Vector) Matrix {
	return Matrix(rows)
}

func (m Matrix) Column(i int) Vector {
	v := make(Vector, 0, len(m))
	for _, row := range m {
		v = append(v, row[i])
	}
	return v
}

func (m Matrix) Transpose() Matrix {
	nbRow := len(m[0])
	m2 := make(Matrix, 0, nbRow)
	for i := 0; i < nbRow; i++ {
		m2 = append(m2, m.Column(i))
	}
	return m2
}

func (m Matrix) Mul(v Vector) Vector {
	var out Vector
	for _, row := range m {
		var acc Value
		for i := range row {
			acc += row[i] * v[i]
		}
		out = append(out, acc)
	}
	return out
}

func (v Vector) Hadamard(v2 Vector) Vector {
	out := make(Vector, 0, len(v))
	for i := range v {
		out = append(out, v[i]*v2[i])
	}
	return out
}

func (v Vector) Sub(v2 Vector) Vector {
	out := make(Vector, len(v))
	for i := range v {
		out = append(out, v[i]-v2[i])
	}
	return out
}

func (v Vector) IsZero() bool {
	for _, e := range v {
		if e != 0 {
			return false
		}
	}
	return true
}

// GeneratePowersCommit returns { g^shift * s^i} for i=0...power included
func GeneratePowersCommit(base Commit, e Element, shift Element, power int) []Commit {
	var gi = make([]Commit, 0, power+1)
	gi = append(gi, base.Clone().Mul(shift, nil))
	var si = one.Clone()
	var tmp = NewElement()
	for i := 0; i < power; i++ {
		// s * (tmp) = s * ( s * ( .. ) )
		si = si.Mul(si, e)
		// s^i * shift
		tmp = tmp.Mul(si, shift)
		gi = append(gi, base.Clone().Mul(tmp, nil))
	}
	return gi
}
