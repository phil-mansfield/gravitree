package gravitree

import (
	"math"
)

/* NOTE: Many of the loops in this file are very very hot, so some
optimizations are a bit aggressive. */

// Potential computes the potential for each point in the tree. Potentials are
// written to the array phi and are in units where G * mp = 1.

//*-- interface --*////

func (phi Potential) Len() int {
	return len(phi)
}

func (phi Potential) AddMonopoleOneSided(
	i int, xi, xj *[3]float64, weight, eps2 float64) {

	dr2 := 0.0
	for k := 0; k < 3; k++ {
		dx := xi[k] - xj[k]
		dr2 += dx * dx
	}
	dr2 += eps2

	phi[i] -= weight / math.Sqrt(dr2)
}

func (phi Potential) AddMonopoleTwoSided(
	i, j int, xi, xj *[3]float64, weight, eps2 float64) {

	dr2 := 0.0
	for k := 0; k < 3; k++ {
		dx := xi[k] - xj[k]
		dr2 += dx * dx
	}
	dr2 += eps2

	val := weight / math.Sqrt(dr2)
	phi[i] -= val
	phi[j] -= val
}
