package gravitree

import (
	"math"
)

//*-- interface --*////

func (tensor TidalTensor) Len() int {
	return len(tensor)
}

func (tensor TidalTensor) AddMonopoleOneSided(
	i int, xi, xj *[3]float64, weight, eps2 float64) {

	x := []float64{0, 0, 0}
	dr2 := 0.0
	for k := 0; k < 3; k++ {
		x[k] = xi[k] - xj[k]
		dr2 += x[k] * x[k]
	}
	dr2 += eps2
	denom := dr2 * dr2 * math.Sqrt(dr2)

	x2 := (xi[0] - xj[0]) * (xi[0] - xj[0])
	y2 := (xi[1] - xj[1]) * (xi[1] - xj[1])
	z2 := (xi[2] - xj[2]) * (xi[2] - xj[2])

	dtensor := [3][3]float64{
		{y2 + z2 - 2*x2, -3 * x[0] * x[1], -3 * x[0] * x[2]},
		{-3 * x[0] * x[1], x2 + z2 - 2*y2, -3 * x[1] * x[2]},
		{-3 * x[0] * x[2], -3 * x[1] * x[2], x2 + y2 - 2*z2},
	}

	// rescale
	for k := 0; k < 3; k++ {
		for l := 0; l < 3; l++ {
			tensor[i][k][l] += weight * dtensor[k][l] / denom
		}
	}
}

func (tensor TidalTensor) AddMonopoleTwoSided(
	i, j int, xi, xj *[3]float64, weight, eps2 float64) {

	x := []float64{0, 0, 0}
	dr2 := 0.0
	for k := 0; k < 3; k++ {
		x[k] = xi[k] - xj[k]
		dr2 += x[k] * x[k]
	}
	dr2 += eps2
	denom := dr2 * dr2 * math.Sqrt(dr2)

	x2 := (xi[0] - xj[0]) * (xi[0] - xj[0])
	y2 := (xi[1] - xj[1]) * (xi[1] - xj[1])
	z2 := (xi[2] - xj[2]) * (xi[2] - xj[2])

	dtensor := [3][3]float64{
		{y2 + z2 - 2*x2, -3 * x[0] * x[1], -3 * x[0] * x[2]},
		{-3 * x[0] * x[1], x2 + z2 - 2*y2, -3 * x[1] * x[2]},
		{-3 * x[0] * x[2], -3 * x[1] * x[2], x2 + y2 - 2*z2},
	}

	// rescale
	for k := 0; k < 3; k++ {
		for l := 0; l < 3; l++ {
			tensor[i][k][l] += weight * dtensor[k][l] / denom
			tensor[j][k][l] += weight * dtensor[k][l] / denom
		}
	}
}
