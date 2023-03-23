package gravitree

import (
	"math"
)

func BindingEnergy(x, v [][3]float64, mp, eps float64, E []float64) {
	IterativeBindingEnergy(x, v, mp, eps, 1, E)
}

func IterativeBindingEnergy(
	x, v [][3]float64, mp, eps float64, iters int, E []float64,
) {
	ok := make([]bool, len(x))
	nPrev := len(x)
	for i := range ok { ok[i] = true }

	for j := 0; j < iters || nPrev == 0; j++ {
		bindingEnergy(x, v, ok, mp, eps, E)

		nCurr := 0
		for i := range ok {
			ok[i] = ok[i] && E[i] < 0
			if ok[i] { nCurr++ }
		}

		if nCurr == nPrev { break }
		nPrev = nCurr
	}

	inf := math.Inf(+1)
	for i := range ok {
		if !ok[i] { E[i] = inf }
	}
}

func bindingEnergy(
	x, v [][3]float64, ok []bool, mp, eps float64, E []float64,
) {
	n0 := 0
	for i := range ok {
		if ok[i] { n0++ }
	}
	
	x0 := make([][3]float64, len(x))
	v0 := make([][3]float64, len(v))
	E0 := make([]float64, len(x))

	j := 0
	for i := range x {
		if ok[i] {
			x0[j], v0[j], E0[j] = x[i], v[i], E[i]
			j++
		}
	}

	tree := NewTree(x0)
	tree.Potential(eps, E0)

	j = 0
	for i := range ok {
		if ok[i] {
			v2 := 0.0
			for dim := 0; dim < 3; dim++ {
				v2 += v0[j][dim]*v0[j][dim]
			}
			E[i] = E0[j]*mp*4.301e-6 + v2/2
			j++
		}
	}
}

func PotentialEnergy(x [][3]float64, mp, eps float64, E []float64) {
	tree := NewTree(x)
    tree.Potential(eps, E)

	for i := range E {
		E[i] *= mp*4.301e-6
	}
}
