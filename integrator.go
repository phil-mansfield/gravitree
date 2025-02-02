package gravitree

import (
	"math"

	"github.com/phil-mansfield/gravitree/utils"
	"gonum.org/v1/gonum/mathext"
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
	ParticleMass  float64
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
	mp float64,
) {

	// pos[step, :] = pos[step-1, :] + 0.5 * dt * vel[step-1, :]
	// acc = get_accelerations(pos, step)

	// if _acc is not None:
	//     _acc[step, :] = acc

	// vel[step, :] = vel[step-1, :] + dt * acc
	// pos[step, :] += 0.5 * dt * vel[step, :]

	acc := Acceleration(make([][3]float64, len(pos)))

	// Perform first half-kick
	for i := 0; i < len(pos); i++ {
		if ok[i] {
			for j := 0; j < 3; j++ {
				// mass = 1 / npts
				pos[i][j] += (dt / 2.0) * vel[i][j]
			}
		}
	}

	if brute {
		BruteForceAccelerationAt(eps, tree.Points, pos, acc)
	} else {
		tracer.Tree.Points = pos
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	// Perform kick step
	for i := 0; i < len(pos); i++ {
		if ok[i] {
			for j := 0; j < 3; j++ {
				vel[i][j] += dt * acc[i][j] * mp
			}
		}
	}

	// Perform second half-drift
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				pos[i][j] += (dt / 2.0) * vel[i][j]
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
	mp float64,
) {

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
	kv1 := utils.RescalePoints(acc, dt*mp)

	// phase space position half a step forward
	wx1 := utils.PointwiseAdd(pos, utils.RescalePoints(kx1, 0.5))
	wv1 := utils.PointwiseAdd(vel, utils.RescalePoints(kv1, 0.5))

	// clear accelerations
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}

	// calculate accelerations for second step
	if brute {
		BruteForceAccelerationAt(eps, tree.Points, wx1, acc)
	} else {
		tracer.Tree.Points = wx1
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	kx2 := utils.RescalePoints(wv1, dt)
	kv2 := utils.RescalePoints(acc, dt*mp)

	// move half-step forward
	wx2 := utils.PointwiseAdd(pos, utils.RescalePoints(kx2, 0.5))
	wv2 := utils.PointwiseAdd(vel, utils.RescalePoints(kv2, 0.5))

	// calculate forces for the third step
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	if brute {
		BruteForceAccelerationAt(eps, tree.Points, wx2, acc)
	} else {
		tracer.Tree.Points = wx2
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	kx3 := utils.RescalePoints(wv2, dt)
	kv3 := utils.RescalePoints(acc, dt*mp)

	// step 4
	wx3 := utils.PointwiseAdd(pos, kx3)
	wv3 := utils.PointwiseAdd(vel, kv3)

	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}

	if brute {
		BruteForceAccelerationAt(eps, tree.Points, wx3, acc)
	} else {
		tracer.Tree.Points = wx3
		tracer.Update()
		tree.EvaluateAt(&tracer.Tree, eps, acc)
	}

	kx4 := utils.RescalePoints(wv3, dt)
	kv4 := utils.RescalePoints(acc, dt*mp) // note this is full step

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

func GetParticleCount(x [][3]float64, r float64) int {

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

func GetCircularVelocity(x [][3]float64, r, particleMass float64) float64 {
	npart := float64(GetParticleCount(x, r))
	return math.Sqrt(1. * npart * particleMass / r)
}

// Smooth potential

func (e *Einasto) EnclosedMass(r float64) float64 {
	// Mtot is the total mass
	// since code units are [m] = 1/npts, then
	// Mtot = 1 when all points are enclosed.
	// return lower incomplete gamma function

	res := 0.0
	a := 3. / e.Alpha
	x := (2. / e.Alpha) * math.Pow(r/e.Rs, e.Alpha)
	res = (1 - mathext.GammaIncRegComp(a, x))
	return res
}

func (e *Einasto) GetCircularVelocity(r float64) float64 {
	res := math.Sqrt(e.EnclosedMass(r) / r)
	return res
}

func (e *Einasto) GetAcceleration(x, acc [][3]float64) {
	for i := range x {
		r := utils.GetNorm(x[i])
		for k := 0; k < 3; k++ {
			acc[i][k] -= x[i][k] * e.EnclosedMass(r) / (r * r * r)
		}
	}
}

func (e *Einasto) GetDConstant() float64 {
	// eqn 16 from https://www.aanda.org/articles/aa/full_html/2012/04/aa18543-11/aa18543-11.html
	n := 1 / e.Alpha
	d := (3. * n)
	d += -(1. / 3.)
	d += (8. / (1215. * n))
	d += (184. / (229635. * n * n))
	d += (1084. / (31000725 * n * n * n))
	d += -(17557576. / (1242974068875. * n * n * n * n))
	return d
}

func (e *Einasto) MathCalF(r float64) float64 {
	x := r / e.Rs
	t1 := (1. / x) * (1. - mathext.GammaIncRegComp(3./e.Alpha, (2./e.Alpha)*math.Pow(x, e.Alpha)))
	t2a := math.Gamma(2. / e.Alpha)
	t2b := mathext.GammaIncRegComp(2./e.Alpha, (2./e.Alpha)*math.Pow(x, e.Alpha)) * math.Gamma(2./e.Alpha) / math.Gamma(3./e.Alpha)
	return t1 + (math.Pow((2./e.Alpha), 1./e.Alpha)/math.Gamma(3./e.Alpha))*(t2a+t2b)
}

func (e *Einasto) GetPotential(r float64) float64 {
	// eqn 19 from https://www.aanda.org/articles/aa/full_html/2012/04/aa18543-11/aa18543-11.html
	h := e.Rs / math.Pow(2./e.Alpha, 1./e.Alpha)
	x := math.Pow(r/h, e.Alpha)
	return -(1. / h) * (mathext.GammaIncReg(3./e.Alpha, x)/math.Pow(x, 1./e.Alpha) + (math.Gamma(2./e.Alpha) * mathext.GammaIncRegComp(2./e.Alpha, x) / math.Gamma(3./e.Alpha)))

}

func (e *Einasto) GetEnergy(x, v [][3]float64) []float64 {
	// Calculates the energy of a tracer orbiting
	// the einasto distribution

	res := make([]float64, len(x))

	for i := 0; i < len(x); i++ {
		// kinetic term
		res[i] = utils.GetNorm(v[i]) * utils.GetNorm(v[i]) / 2.0
		// potential term
		res[i] += e.GetPotential(utils.GetNorm(x[i]))
	}

	return res
}

func (prof *Einasto) LeapfrogStep(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	dt float64, // timestep (see units below)
) {

	acc := Acceleration(make([][3]float64, len(pos)))
	prof.GetAcceleration(pos, acc)

	// Perform first half-kick
	for i := 0; i < len(pos); i++ {
		for j := 0; j < 3; j++ {
			// mass = 1 / npts
			vel[i][j] += (dt / 2.0) * (acc)[i][j]
		}
	}

	// Perform drift step
	for i := 0; i < len(pos); i++ {
		for j := 0; j < 3; j++ {
			pos[i][j] += dt * vel[i][j]
		}
	}

	// clear acc before recalculating!
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(pos, acc)

	// Perform second half-kick
	for i := range pos {
		for j := range 3 {
			vel[i][j] += (dt / 2.0) * (acc)[i][j]
		}
	}
}

func (prof *Einasto) RK4Step(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	dt float64, // timestep (see units below)
) {

	// create acceleration in memory and evaluate force
	acc := Acceleration(make([][3]float64, len(pos)))
	prof.GetAcceleration(pos, acc)
	//
	kx1 := utils.RescalePoints(vel, dt)
	kv1 := utils.RescalePoints(acc, dt)

	// phase space position half a step forward
	wx1 := utils.PointwiseAdd(pos, utils.RescalePoints(kx1, 0.5))
	wv1 := utils.PointwiseAdd(vel, utils.RescalePoints(kv1, 0.5))

	// clear accelerations
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(wx1, acc)

	kx2 := utils.RescalePoints(wv1, dt)
	kv2 := utils.RescalePoints(acc, dt)

	// move half-step forward
	wx2 := utils.PointwiseAdd(pos, utils.RescalePoints(kx2, 0.5))
	wv2 := utils.PointwiseAdd(vel, utils.RescalePoints(kv2, 0.5))

	// calculate forces for the third step
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(wx2, acc)

	kx3 := utils.RescalePoints(wv2, dt)
	kv3 := utils.RescalePoints(acc, dt)

	// step 4
	wx3 := utils.PointwiseAdd(pos, kx3)
	wv3 := utils.PointwiseAdd(vel, kv3)

	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(wx3, acc)

	kx4 := utils.RescalePoints(wv3, dt)
	kv4 := utils.RescalePoints(acc, dt) // note this is full step

	//

	for i := range pos {
		for j := range 3 {
			pos[i][j] += (1. / 6.) * (kx1[i][j] + 2*kx2[i][j] + 2*kx3[i][j] + kx4[i][j])
			vel[i][j] += (1. / 6.) * (kv1[i][j] + 2*kv2[i][j] + 2*kv3[i][j] + kv4[i][j])
		}
	}
}

func (prof *Plummer) RK4Step(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	dt float64, // timestep (see units below)
) {

	// create acceleration in memory and evaluate force
	acc := Acceleration(make([][3]float64, len(pos)))
	prof.GetAcceleration(pos, acc)
	//
	kx1 := utils.RescalePoints(vel, dt)
	kv1 := utils.RescalePoints(acc, dt)

	// phase space position half a step forward
	wx1 := utils.PointwiseAdd(pos, utils.RescalePoints(kx1, 0.5))
	wv1 := utils.PointwiseAdd(vel, utils.RescalePoints(kv1, 0.5))

	// clear accelerations
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(wx1, acc)

	kx2 := utils.RescalePoints(wv1, dt)
	kv2 := utils.RescalePoints(acc, dt)

	// move half-step forward
	wx2 := utils.PointwiseAdd(pos, utils.RescalePoints(kx2, 0.5))
	wv2 := utils.PointwiseAdd(vel, utils.RescalePoints(kv2, 0.5))

	// calculate forces for the third step
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(wx2, acc)

	kx3 := utils.RescalePoints(wv2, dt)
	kv3 := utils.RescalePoints(acc, dt)

	// step 4
	wx3 := utils.PointwiseAdd(pos, kx3)
	wv3 := utils.PointwiseAdd(vel, kv3)

	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	prof.GetAcceleration(wx3, acc)

	kx4 := utils.RescalePoints(wv3, dt)
	kv4 := utils.RescalePoints(acc, dt) // note this is full step

	//

	for i := range pos {
		for j := range 3 {
			pos[i][j] += (1. / 6.) * (kx1[i][j] + 2*kx2[i][j] + 2*kx3[i][j] + kx4[i][j])
			vel[i][j] += (1. / 6.) * (kv1[i][j] + 2*kv2[i][j] + 2*kv3[i][j] + kv4[i][j])
		}
	}
}

// Plummer

func (p *Plummer) EnclosedMass(r float64) float64 {
	// M0 is the total mass, in this package it is currently Mtot = 1.
	return (r * r * r) / math.Pow(p.B*p.B+r*r, 1.5)
}

func (p *Plummer) GetCircularVelocity(r float64) float64 {
	res := math.Sqrt(p.EnclosedMass(r) / r)
	return res
}

func (p *Plummer) GetAcceleration(x, acc [][3]float64) {
	for i := range x {
		r := utils.GetNorm(x[i])
		for k := 0; k < 3; k++ {
			acc[i][k] -= x[i][k] * p.EnclosedMass(r) / (r * r * r)
		}
	}
}

func (p *Plummer) GetPotential(r float64) float64 {
	return -1 / math.Sqrt(p.B*p.B+r*r)
}

func (p *Plummer) GetEnergy(x, v [][3]float64) []float64 {
	// Calculates the energy of a tracer orbiting
	// the Plummer distribution

	res := make([]float64, len(x))

	for i := 0; i < len(x); i++ {
		// kinetic term
		res[i] = utils.GetNorm(v[i]) * utils.GetNorm(v[i]) / 2.0
		// potential term
		res[i] += p.GetPotential(utils.GetNorm(x[i]))
	}

	return res
}

func (prof *Plummer) LeapfrogStep(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	dt float64, // timestep (see units below)
) {

	acc := Acceleration(make([][3]float64, len(pos)))
	prof.GetAcceleration(pos, acc)

	// Perform first half-kick
	for i := 0; i < len(pos); i++ {
		for j := 0; j < 3; j++ {
			// mass = 1 / npts
			vel[i][j] += (dt / 2.0) * (acc)[i][j]
		}
	}

	// Perform drift step
	for i := 0; i < len(pos); i++ {
		for j := 0; j < 3; j++ {
			pos[i][j] += dt * vel[i][j]
		}
	}

	// clear acc before recalculating!
	acc = Acceleration(make([][3]float64, len(pos)))
	prof.GetAcceleration(pos, acc)

	// Perform second half-kick
	for i := range pos {
		for j := range 3 {
			vel[i][j] += (dt / 2.0) * (acc)[i][j]
		}
	}
}

// Point source

func (p *PointMass) EnclosedMass(r float64) float64 {
	return p.Mass
}

func (p *PointMass) GetCircularVelocity(r float64) float64 {
	return math.Sqrt(p.Mass / r)
}

func (p *PointMass) GetAcceleration(x, acc [][3]float64) {
	for i := range x {
		r := utils.GetNorm(x[i])
		for k := 0; k < 3; k++ {
			acc[i][k] -= x[i][k] * p.Mass / (r * r * r)
		}
	}
}

func (p *PointMass) GetPotential(r float64) float64 {
	return -p.Mass / r
}

func (p *PointMass) GetEnergy(x, v [][3]float64) []float64 {
	res := make([]float64, len(x))
	for i := range x {
		// kinetic term
		res[i] = utils.GetNorm(v[i]) * utils.GetNorm(v[i]) / 2.0
		// potential term
		res[i] += p.GetPotential(utils.GetNorm(x[i]))
	}
	return res
}

func (p *PointMass) LeapfrogStep(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	dt float64, // timestep (see units below)
) {
	acc := Acceleration(make([][3]float64, len(pos)))
	p.GetAcceleration(pos, acc)

	// Perform first half-kick
	for i := 0; i < len(pos); i++ {
		for j := 0; j < 3; j++ {
			vel[i][j] += (dt / 2.0) * acc[i][j]
		}
	}

	// Perform drift step
	for i := 0; i < len(pos); i++ {
		for j := 0; j < 3; j++ {
			pos[i][j] += dt * vel[i][j]
		}
	}

	// clear acc before recalculating!
	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	p.GetAcceleration(pos, acc)

	// Perform second half-kick
	for i := range pos {
		for j := range 3 {
			vel[i][j] += (dt / 2.0) * acc[i][j]
		}
	}
}

func (p *PointMass) RK4Step(
	pos [][3]float64, // positions
	vel [][3]float64, // velocities
	dt float64, // timestep (see units below)
) {
	acc := Acceleration(make([][3]float64, len(pos)))
	p.GetAcceleration(pos, acc)

	kx1 := utils.RescalePoints(vel, dt)
	kv1 := utils.RescalePoints(acc, dt)

	wx1 := utils.PointwiseAdd(pos, utils.RescalePoints(kx1, 0.5))
	wv1 := utils.PointwiseAdd(vel, utils.RescalePoints(kv1, 0.5))

	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	p.GetAcceleration(wx1, acc)

	kx2 := utils.RescalePoints(wv1, dt)
	kv2 := utils.RescalePoints(acc, dt)

	wx2 := utils.PointwiseAdd(pos, utils.RescalePoints(kx2, 0.5))
	wv2 := utils.PointwiseAdd(vel, utils.RescalePoints(kv2, 0.5))

	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	p.GetAcceleration(wx2, acc)

	kx3 := utils.RescalePoints(wv2, dt)
	kv3 := utils.RescalePoints(acc, dt)

	wx3 := utils.PointwiseAdd(pos, kx3)
	wv3 := utils.PointwiseAdd(vel, kv3)

	for i := range acc {
		for j := range acc[i] {
			acc[i][j] = 0.0
		}
	}
	p.GetAcceleration(wx3, acc)

	kx4 := utils.RescalePoints(wv3, dt)
	kv4 := utils.RescalePoints(acc, dt)

	for i := range pos {
		for j := range 3 {
			pos[i][j] += (1. / 6.) * (kx1[i][j] + 2*kx2[i][j] + 2*kx3[i][j] + kx4[i][j])
			vel[i][j] += (1. / 6.) * (kv1[i][j] + 2*kv2[i][j] + 2*kv3[i][j] + kv4[i][j])
		}
	}
}
