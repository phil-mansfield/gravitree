package gravitree

import "github.com/phil-mansfield/gravitree/utils"

// Re-exporting for exposure to testing
var ReadPointFile = utils.ReadPointFile
var ReadPhaseSpaceFile = utils.ReadPhaseSpaceFile
var RescalePoints = utils.RescalePoints

type Profile interface {
	EnclosedMass(r float64) float64
	GetCircularVelocity(r float64) float64
	GetAcceleration(x, acc [][3]float64)
	GetPotential(r float64) float64
	GetEnergy(x, v [][3]float64)
	LeapfrogStep(pos, vel [][3]float64, dt float64)
	RK4Step(pos, vel [][3]float64, dt float64)
}

type Einasto struct {
	Rs    float64
	Alpha float64
}

type Plummer struct {
	B float64
}

type PointMass struct {
	Mass float64
}
