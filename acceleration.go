package gravitree

import (
	"math"
)

// Implementations of acceleration vectors

// Acceleration computed for each point in the tree. Vectors are
// written to the 3-vector acc and are in units where G * mp = 1.

func (acc Acceleration) Len() int {
	return len(acc)
}

func (acc Acceleration) AddMonopoleOneSided(
	i int, xi, xj *[3]float64, weight, eps2 float64) {

	x := []float64{0, 0, 0}
	dr2 := 0.0
	for k := 0; k < 3; k++ {
		x[k] = xi[k] - xj[k]
		dr2 += x[k] * x[k]
	}
	dr2 += eps2

	for j := 0; j < 3; j++ {
		acc[i][j] -= weight * x[j] / (dr2 * math.Sqrt(dr2))
	}
}

func (acc Acceleration) AddMonopoleTwoSided(
	i, j int, xi, xj *[3]float64, weight, eps2 float64) {

	x := []float64{0, 0, 0}
	dr2 := 0.0
	for k := 0; k < 3; k++ {
		x[k] = xi[k] - xj[k]
		dr2 += x[k] * x[k]
	}
	dr2 += eps2

	for k := 0; k < 3; k++ {
		acc[i][k] -= weight * x[k] / (dr2 * math.Sqrt(dr2))
		acc[j][k] += weight * x[k] / (dr2 * math.Sqrt(dr2))
	}
}
