package old_gravitree

import (
	"fmt"
	"math"
	"math/rand"
)

// Generate points randomly within the rectangular prism defined by span.
// Span[0] gives the lower corner and span[1] gives the upper corner.
func generateUniformPoints(span [2][3]float64, out [][3]float64) {
	width := [3]float64{ }
	for d := 0; d < 3; d++ { width[d] = span[1][d] - span[0][d] }
	
	for i := range out {
		for d := 0; d < 3; d++ {
			out[i][d] = span[0][d] + rand.Float64()*width[d]
		}
	}
}


// generatePowerLawPoints generates points from a distribution with
// a power law density profile rho(r) ~ r^alpha. Mass profiles diverge if
// alpha <= -3.0.
//
// This function works by applying inverse transfer sampling to the mass profile
// and again to the spherical coordinates to get random cartesian coordinates.
func generatePowerLawPoints(alpha, rMax float64, out [][3]float64) {
	if alpha <= -3.0 {
		panic(fmt.Sprintf(
			"alpha = %f leads to a divergence mass profile.", alpha))
	}

	invMEnc := func (r float64) float64 { return math.Pow(r, alpha + 3) }
	mMax := math.Pow(rMax, alpha + 3)
	
	for i := range out {
		r := randomRadius(invMEnc, mMax)
		theta := math.Acos(2*rand.Float64() - 1)
		phi := 2*math.Pi * rand.Float64()
		
		sinTheta, cosTheta := math.Sincos(theta)
		sinPhi, cosPhi := math.Sincos(phi)
		out[i] = [3]float64 { r*sinTheta*cosPhi, r*sinTheta*sinPhi, r*cosTheta }
	}
}

type func1D func(float64) float64

// randomRadius returns a random radius <= rMax for a mass distribution M(<r).
// The user must supply M(<r)^-1, invMEnc as well as M(<rMax) = mMax.
func randomRadius(invMEnc func1D, mMax float64) float64 {
	M := rand.Float64() * mMax
	return invMEnc(M)
}
