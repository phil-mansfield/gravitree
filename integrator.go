package gravitree

import (
	"math"

	"github.com/phil-mansfield/gravitree/utils"
)

type SimulationOptions struct {
	Steps         int
	Dt            float64
	SaveEvery     int
	Eps           float64
	Brute         bool
	SaveDirectory string
	TreeOptions   TreeOptions
	Integrator    string
}

// Code units:
// mass     : M_vir (i.e., each point has mass 1/npoints)
// length   : R_vir
// velocity : V_circ(R_vir)
// time     : T_circ(R_vir) / 2 pi

func LeapfrogStep(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	tree *Tree, // tree loaded with massive particles
	tracer *ArrayTree, // tree loaded with tracers to evaluate forces on
	ok []bool, // flag to evaluate forces on certain particles
	dt float64, // timestep (see units below)
	eps float64, // force softening (Plummer kernel)
	brute bool, // calculate exact forces?
) {

	npts := float64(len(tree.Points))

	acc := Acceleration(make([][3]float64, len(pos)))

	if brute {
		BruteForceAccelerationAt(eps, tree.Points, pos, acc)
	} else {
		tracer.Tree.Points = pos
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	// Perform first half-kick
	for i := 0; i < len(pos); i++ {
		if ok[i] {
			for j := 0; j < 3; j++ {
				// mass = 1 / npts
				vel[i][j] += (dt / 2.0) * (acc)[i][j] / npts
			}
		}
	}

	// Perform drift step
	for i := 0; i < len(pos); i++ {
		if ok[i] {
			for j := 0; j < 3; j++ {
				pos[i][j] += dt * vel[i][j]
			}
		}
	}

	// clear acc before recalculating!
	acc = Acceleration(make([][3]float64, len(pos)))
	if brute {
		BruteForceAccelerationAt(eps, tree.Points, pos, acc)
	} else {

		tracer.Tree.Points = pos
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	// Perform second half-kick
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				vel[i][j] += (dt / 2.0) * (acc)[i][j] / npts
			}
		}
	}
}

func RK4Step(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	tree *Tree, // tree loaded with massive particles
	tracer *ArrayTree, // tree loaded with tracers to evaluate forces on
	ok []bool, // flag to evaluate forces on certain particles
	dt float64, // timestep (see units below)
	eps float64, // force softening (Plummer kernel)
	brute bool, // calculate exact forces?
) {

	npts := float64(len(tree.Points))

	// create acceleration in memory and evaluate force
	acc := Acceleration(make([][3]float64, len(pos)))

	if brute {
		BruteForceAccelerationAt(eps, tree.Points, pos, acc)
	} else {
		tracer.Tree.Points = pos
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	//
	kx1 := utils.RescalePoints(vel, dt)
	kv1 := utils.RescalePoints(acc, dt/npts)

	// phase space position half a step forward
	wx1 := utils.PointwiseAdd(pos, utils.RescalePoints(kx1, 0.5))
	wv1 := utils.PointwiseAdd(vel, utils.RescalePoints(kv1, 0.5))

	// clear accelerations
	acc = Acceleration(make([][3]float64, len(pos)))

	// calculate accelerations for second step
	if brute {
		BruteForceAccelerationAt(eps, tree.Points, wx1, acc)
	} else {
		tracer.Tree.Points = wx1
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	kx2 := utils.RescalePoints(wv1, dt)
	kv2 := utils.RescalePoints(acc, dt/npts)

	// move half-step forward
	wx2 := utils.PointwiseAdd(pos, utils.RescalePoints(kx2, 0.5))
	wv2 := utils.PointwiseAdd(vel, utils.RescalePoints(kv2, 0.5))

	// calculate forces for the third step
	acc = Acceleration(make([][3]float64, len(pos)))

	if brute {
		BruteForceAccelerationAt(eps, tree.Points, wx2, acc)
	} else {
		tracer.Tree.Points = wx2
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	kx3 := utils.RescalePoints(wv2, dt)
	kv3 := utils.RescalePoints(acc, dt/npts)

	// step 4
	wx3 := utils.PointwiseAdd(pos, kx3)
	wv3 := utils.PointwiseAdd(vel, kv3)

	acc = Acceleration(make([][3]float64, len(pos)))

	if brute {
		BruteForceAccelerationAt(eps, tree.Points, wx3, acc)
	} else {
		tracer.Tree.Points = wx3
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	kx4 := utils.RescalePoints(wv3, dt)
	kv4 := utils.RescalePoints(acc, dt/npts) // note this is full step

	//

	for i := range pos {
		if ok[i] {
			for j := range 3 {
				pos[i][j] += (1. / 6.) * (kx1[i][j] + 2*kx2[i][j] + 2*kx3[i][j] + kx4[i][j])
				vel[i][j] += (1. / 6.) * (kv1[i][j] + 2*kv2[i][j] + 2*kv3[i][j] + kv4[i][j])
			}
		}
	}
}

func getParticleCount(x [][3]float64, r float64) int {

	count := 0

	for i := 0; i < len(x); i++ {
		ri := 0.0

		for k := 0; k < 3; k++ {
			ri += x[i][k] * x[i][k]
		}

		ri = math.Sqrt(ri)

		if ri <= r {
			count++
		}
	}
	return count
}

func CalculateEnergy(
	pos [][3]float64,
	vel [][3]float64,
	extPos [][3]float64,
	eps float64,
	ok []bool,
) []float64 {

	// Calculates the energy of a tracer orbiting
	// a distribution of points `extPos`

	res := make([]float64, len(pos))
	pot := Potential(make([]float64, len(pos)))
	BruteForcePotentialAt(eps, extPos, pos, pot)

	for i := 0; i < len(pos); i++ {
		if ok[i] {
			// kinetic term
			res[i] = utils.GetNorm(vel[i]) * utils.GetNorm(vel[i]) / 2.0
			// potential term
			res[i] += pot[i]
		}
	}

	return res
}

func GetCircularVelocity(x [][3]float64, r float64) float64 {
	npart := float64(getParticleCount(x, r))
	ntot := float64(len(x))
	return math.Sqrt(1. * npart / (r * ntot))
}
