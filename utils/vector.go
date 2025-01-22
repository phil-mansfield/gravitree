package utils

import "math"

func RescalePoints(x [][3]float64, a float64) [][3]float64 {
	res := make([][3]float64, len(x))
	for i := range x {
		for k := 0; k < 3; k++ {
			res[i][k] = a * x[i][k]
		}
	}
	return res
}

func PointwiseAdd(a [][3]float64, b [][3]float64) [][3]float64 {

	if len(a) != len(b) {
		panic("Input vectors are not of the same length")
	}

	res := make([][3]float64, len(a))

	for i := 0; i < len(a); i++ {
		for k := 0; k < 3; k++ {
			res[i][k] = a[i][k] + b[i][k]
		}
	}
	return res
}

func GetNorm(x [3]float64) float64 {
	r := 0.0
	for i := 0; i < 3; i++ {
		r += x[i] * x[i]
	}
	return math.Sqrt(r)
}
