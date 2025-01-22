package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strconv"

	"github.com/phil-mansfield/gravitree"
)

func writeOrbit(filename string, x [][3]float64, v [][3]float64, ok []bool) {
	f, err := os.Create(filename)
	if err != nil {
		panic(err.Error())
	}
	for i := range x {
		if ok[i] {
			fmt.Fprintf(f, "%.10g %.10g %.10g %.10g %.10g %.10g \n",
				x[i][0], x[i][1], x[i][2], v[i][0], v[i][1], v[i][2])
		}
	}
}

func writeQuantity(filename string, q []float64) {
	// writes a float64 quantity `q` to a file
	// source: https://stackoverflow.com/questions/50205448/whats-the-difference-between-os-o-append-and-os-modeappend
	f, err := os.OpenFile(filename, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		panic(err.Error())
	}
	// write to file
	for i := 0; i < len(q); i++ {
		fmt.Fprintf(f, "%v ", q[i])
	}
	// demarcate
	fmt.Fprintf(f, "\n")
}

func MakeDirectory(path string) {
	_, err := os.Stat(path)

	if os.IsNotExist(err) {
		err := os.MkdirAll(path, os.ModePerm)
		if err != nil {
			panic("Unable to create snapshot folder.")
		}
	}
}

func SimulateCircularOrbits(
	npart int,
	r0 float64,
	dr float64,
	extPos [][3]float64,
	opt gravitree.SimulationOptions,
) {

	orbitFile := filepath.Join(opt.SaveDirectory, "circ_t_%d.dat")
	energyFile := filepath.Join(opt.SaveDirectory, "circ_e.dat")
	timeFile := filepath.Join(opt.SaveDirectory, "circ_t.dat")

	if opt.Brute {
		orbitFile = filepath.Join(opt.SaveDirectory, "circ_bf_t_%d.dat") //
		energyFile = filepath.Join(opt.SaveDirectory, "circ_bf_e.dat")
		fmt.Printf("Calculating brute force... \n")
	}

	pos := make([][3]float64, npart)
	vel := make([][3]float64, npart)
	tCirc := make([]float64, npart)
	ok := make([]bool, npart)

	for i := 0; i < npart; i++ {
		ok[i] = true
		pos[i] = [3]float64{r0 + float64(i)*dr, 0, 0}
		vel[i] = [3]float64{0, gravitree.GetCircularVelocity(extPos, r0+float64(i)*dr), 0}
		tCirc[i] = 2 * math.Pi * (r0 + float64(i)*dr) / gravitree.GetCircularVelocity(extPos, r0+float64(i)*dr)
	}

	tree := gravitree.NewTree(extPos)

	tracer := gravitree.NewArrayTree(pos)

	t := 0.0

	_, err := os.Create(energyFile)
	if err != nil {
		panic(err.Error())
	}

	_, err = os.Create(timeFile)
	if err != nil {
		panic(err.Error())
	}

	for i := range opt.Steps {

		if i%opt.SaveEvery == 0 {
			fmt.Printf("(circ) step %v; brute: %t \n", i, opt.Brute)
			writeOrbit(fmt.Sprintf(orbitFile, i/opt.SaveEvery), pos, vel, ok)
			e := gravitree.CalculateEnergy(pos, vel, extPos, opt.Eps, ok)

			writeQuantity(energyFile, e)

			ti := make([]float64, npart)

			for k := 0; k < npart; k++ {
				ti[k] = t / tCirc[k]
			}

			writeQuantity(timeFile, ti)

		}

		if opt.Integrator == "leapfrog" {
			gravitree.LeapfrogStep(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else if opt.Integrator == "rk4" {
			gravitree.RK4Step(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt
	}

}

func SimulatePointOrbit(opt gravitree.SimulationOptions) {
	orbitFile := filepath.Join(opt.SaveDirectory, "point_t_%d.dat")
	energyFile := filepath.Join(opt.SaveDirectory, "point_e.dat")
	timeFile := filepath.Join(opt.SaveDirectory, "point_t.dat")

	if opt.Brute {
		orbitFile = filepath.Join(opt.SaveDirectory, "point_bf_t_%d.dat")
		energyFile = filepath.Join(opt.SaveDirectory, "point_bf_e.dat")
		fmt.Printf("Calculating using brute force. \n")
	}

	pos := make([][3]float64, 1)
	vel := make([][3]float64, 1)
	ok := make([]bool, 1)

	extPos := [][3]float64{{0., 0., 0.}}
	tree := gravitree.NewTree(extPos)

	ok[0] = true
	pos[0] = [3]float64{1., 0, 0}
	vel[0] = [3]float64{0, 1., 0}

	tracer := gravitree.NewArrayTree(pos)

	t := 0.0

	_, err := os.Create(energyFile)
	if err != nil {
		panic(err.Error())
	}

	_, err = os.Create(timeFile)
	if err != nil {
		panic(err.Error())
	}

	for i := range opt.Steps {

		if i%opt.SaveEvery == 0 {
			fmt.Printf("(single) step %v; brute: %t \n", i, opt.Brute)
			writeOrbit(fmt.Sprintf(orbitFile, i/opt.SaveEvery), pos, vel, ok)
			e := gravitree.CalculateEnergy(pos, vel, extPos, opt.Eps, ok)
			writeQuantity(energyFile, e)
			ti := make([]float64, 1)
			ti[0] = t / (2 * math.Pi) // 2 pi
			writeQuantity(timeFile, ti)
		}

		if opt.Integrator == "leapfrog" {
			gravitree.LeapfrogStep(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else if opt.Integrator == "rk4" {
			gravitree.RK4Step(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt
	}
}

func SimulateDipoleOrbit(r, dFraction float64, opt gravitree.SimulationOptions) {
	orbitFile := "snapshots/dipole_t_%d.dat"
	energyFile := "snapshots/dipole_e.dat"

	if opt.Brute {
		orbitFile = "snapshots/dipole_bf_t_%d.dat"
		energyFile = "snapshots/dipole_bf_e.dat"
		fmt.Printf("Calculating using brute force. \n")
	}

	pos := make([][3]float64, 1)
	vel := make([][3]float64, 1)
	ok := make([]bool, 1)

	extPos := [][3]float64{{r * dFraction, 0., 0.}, {-r * dFraction, 0., 0.}}
	tree := gravitree.NewTree(extPos)

	ok[0] = true
	pos[0] = [3]float64{r, 0, 0}
	vel[0] = [3]float64{0, gravitree.GetCircularVelocity(extPos, r), 0}

	tracer := gravitree.NewArrayTree(pos)

	t := 0.0

	_, err := os.Create(energyFile)
	if err != nil {
		panic(err.Error())
	}

	for i := range opt.Steps {

		if i%opt.SaveEvery == 0 {
			fmt.Printf("(single) step %v; brute: %t \n", i, opt.Brute)
			writeOrbit(fmt.Sprintf(orbitFile, i/opt.SaveEvery), pos, vel, ok)
			e := gravitree.CalculateEnergy(pos, vel, extPos, opt.Eps, ok)
			writeQuantity(energyFile, e)
		}

		if opt.Integrator == "leapfrog" {
			gravitree.LeapfrogStep(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else if opt.Integrator == "rk4" {
			gravitree.RK4Step(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt
	}
}

func cartToSph(pos [3]float64) [3]float64 {
	res := [3]float64{0, 0, 0}
	x, y, z := pos[0], pos[1], pos[2]
	res[0] = math.Sqrt(x*x + y*y + z*z) // this is the radius
	res[1] = math.Acos(z / res[0])
	res[2] = math.Atan2(y, x)
	return res
}

func sphToCart(sph [3]float64) [3]float64 {
	res := [3]float64{0, 0, 0}
	r, theta, phi := sph[0], sph[1], sph[2]

	res[0] = r * math.Sin(theta) * math.Cos(phi)
	res[1] = r * math.Sin(theta) * math.Sin(phi)
	res[2] = r * math.Cos(theta)

	return res
}

func sampleNorm(mu, sigma float64) float64 {
	return mu + sigma*rand.NormFloat64()
}

func getSprayPair(x0, v0 [3]float64, rt float64) [][3]float64 {
	q_sph := cartToSph(x0)

	q1 := [3]float64{0., 0., 0.}
	q2 := [3]float64{0., 0., 0.}

	p1 := [3]float64{0., 0., 0.}
	p2 := [3]float64{0., 0., 0.}

	r0 := 0.0

	for k := range 3 {
		q1[k] = q_sph[k]
		q2[k] = q_sph[k]
		p1[k] = v0[k]
		p2[k] = v0[k]
		r0 += x0[k] * x0[k]
	}

	// r0 = math.Sqrt(r0)

	// modify tangential velocity
	// vt ~
	// v_t := (v_circ * rt / r0)

	q1[0] -= rt * sampleNorm(2.0, 0.4)
	q2[0] += rt * sampleNorm(2.0, 0.4)

	q1 = sphToCart(q1)
	q2 = sphToCart(q2)

	return [][3]float64{q1, p1, q2, p2}
}

func SimulateStream(
	x0 [3]float64,
	v0 [3]float64,
	rt float64,
	extPos [][3]float64,
	spawnEvery int,
	opt gravitree.SimulationOptions,
) {
	orbit_fn := "snapshots/stream_t_%d.dat"

	if opt.Brute {
		orbit_fn = "snapshots/stream_bf_t_%d.dat"
		fmt.Printf("Calculating using brute force. \n")
	}

	nStream := int(opt.Steps/spawnEvery)*2 + 1

	fmt.Printf("stream count %v \n", nStream)

	pos := make([][3]float64, nStream)
	vel := make([][3]float64, nStream)
	ok := make([]bool, nStream)

	tree := gravitree.NewTree(extPos)
	// fix this

	t := 0.0

	for i := 0; i < nStream; i++ {
		ok[i] = false
		pos[i] = [3]float64{0, 0, 0}
		vel[i] = [3]float64{0, 0, 0}
	}

	for k := 1; k < nStream; k++ {
		pos[k][0] += rand.NormFloat64() * 1000
		pos[k][1] += rand.NormFloat64() * 1000
		pos[k][2] += rand.NormFloat64() * 1000
	}

	ok[0] = true

	// set initial conditions of stream particle
	pos[0] = x0
	vel[0] = v0

	tracer := gravitree.NewArrayTree(pos)

	for i := range opt.Steps {
		fmt.Printf("(stream) step %v; brute: %t \n", i, opt.Brute)

		if int(t/opt.Dt)%opt.SaveEvery == 0 {
			writeOrbit(fmt.Sprintf(orbit_fn, i/opt.SaveEvery), pos, vel, ok)
		}

		if int(t/opt.Dt)%spawnEvery == 0 {
			ok[i/spawnEvery+1] = true
			ok[i/spawnEvery+2] = true

			// v_circ := getCircularVelocity(ext_pos, vecNorm(pos[0]))
			spray := getSprayPair(pos[0], vel[0], rt)

			pos[i/spawnEvery+1] = spray[0]
			vel[i/spawnEvery+1] = spray[1]

			pos[i/spawnEvery+2] = spray[2]
			vel[i/spawnEvery+2] = spray[3]
		}

		if opt.Integrator == "leapfrog" {
			gravitree.LeapfrogStep(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else if opt.Integrator == "rk4" {
			gravitree.RK4Step(pos, vel, tree, tracer, ok, opt.Dt, opt.Eps, opt.Brute)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt
	}

}

func main() {

	method := flag.String("method", "approx", "approx or brute (exact) soln")
	gen := flag.String("gen", "stream", "sim configuration")
	npts := flag.String("npts", "3", "einasto profile to use")
	integrator := flag.String("int", "leapfrog", "integrator")
	feps := flag.String("feps", "1", "softening multiplier")
	dt := flag.String("dt", "1e-4", "timestep")
	theta := flag.String("theta", "1e-2", "opening angle")
	crit := flag.String("crit", "pkd", "bh scheme")
	saveDir := flag.String("saveto", "", "directory to store snapshots")

	// default is pkd
	var opCrit gravitree.OpeningCriteria

	if *crit == "pkd" {
		opCrit = gravitree.PKDGRAV3
	} else if *crit == "bh" {
		opCrit = gravitree.BarnesHut
	} else if *crit == "sw" {
		opCrit = gravitree.SalmonWarren
	} else {
		panic("Invalid opening criteria specified. Supported criteria are `pkd`, `bh`, and `sw`.")
	}

	flag.Parse()

	path := "snapshots"

	brute := false

	if *method == "brute" {
		brute = true
	}

	// obtain static distribution to orbit around
	particle_fn := fmt.Sprintf("../einasto_n=%s_a=18.dat", *npts)
	extPos := gravitree.ReadPointFile(particle_fn)
	extPos = gravitree.RescalePoints(extPos, 1./(0.1))

	x0 := [3]float64{1., 0, 0}
	v0 := [3]float64{0, gravitree.GetCircularVelocity(extPos, 1.), 0}

	// eps—Mimics Symphony softening scale
	// given that the particle distribution
	// represents the DM potential

	// eps := 0.096 / math.Pow(float64(len(extPos)), 1./3.)
	// eps := 0.004
	eps := 0.004 / math.Pow(float64(len(extPos)), 1./3.)

	fepsValue, _ := strconv.ParseFloat(*feps, 64)
	thetaValue, _ := strconv.ParseFloat(*theta, 64)
	dtVal, _ := strconv.ParseFloat(*dt, 64)

	if *saveDir != "" {
		MakeDirectory(*saveDir)
	}

	if *gen == "stream" {

		path = fmt.Sprintf("snapshots_n=%s", *npts)
		path = filepath.Join(*saveDir, path)
		MakeDirectory(path)

		spawn_every := 10000

		opt := gravitree.SimulationOptions{
			Brute: brute, Steps: int(1e6), Dt: 0.0001, SaveEvery: 1000,
			Eps: eps, SaveDirectory: path, Integrator: *integrator,
			TreeOptions: gravitree.TreeOptions{Criteria: opCrit, Theta: 1e-2},
		}

		SimulateStream(x0, v0, 1e-2, extPos, spawn_every, opt)

	} else if *gen == "circ" {

		// TEST:
		// Compare energies against smooth Einasto profile

		path = fmt.Sprintf("snapshots_n=%s_feps=%s_int=%s_dt=%s_th=%s", *npts, *feps, *integrator, *dt, *theta)
		path = filepath.Join(*saveDir, path)

		MakeDirectory(path)

		time := 40.
		steps := int(time / dtVal)

		opt := gravitree.SimulationOptions{
			Brute: brute, Steps: steps, Dt: dtVal, SaveEvery: 10,
			Eps: fepsValue * eps, SaveDirectory: path, Integrator: *integrator,
			TreeOptions: gravitree.TreeOptions{Criteria: opCrit, Theta: thetaValue},
		}

		SimulateCircularOrbits(10, .1, .1, extPos, opt)

	} else if *gen == "single" {

		// TEST:
		// Need to simulate at least 10 orbits
		// 3 orbits ~ 10 time units

		path = fmt.Sprintf("snapshots_int=%s_dt=%s", *integrator, *dt)
		path = filepath.Join(*saveDir, path)

		eps = 0

		dtVal, _ := strconv.ParseFloat(*dt, 64)
		time := 40.

		MakeDirectory(path)
		opt := gravitree.SimulationOptions{
			// ten orbits
			// sample 1/10th of an orbit per snapshot
			Brute: brute, Steps: int(time / dtVal), Dt: dtVal, SaveEvery: 10,
			Eps: eps, SaveDirectory: path, Integrator: *integrator,
			TreeOptions: gravitree.TreeOptions{Criteria: opCrit, Theta: thetaValue},
		}

		SimulatePointOrbit(opt)

	} else if *gen == "dipole" {

		// TEST:
		// See accuracy as a function of dipole length

		path = fmt.Sprintf("snapshots_n=%s_feps=%s_int=%s_dt=%s_th=%s", *npts, *feps, *integrator, *dt, *theta)

		path = filepath.Join(*saveDir, path)
		MakeDirectory(path)

		opt := gravitree.SimulationOptions{
			Brute: brute, Steps: int(1e5), Dt: 0.001, SaveEvery: 100,
			Eps: fepsValue, SaveDirectory: path, Integrator: *integrator,
			TreeOptions: gravitree.TreeOptions{Criteria: gravitree.PKDGRAV3, Theta: thetaValue},
		}

		SimulateDipoleOrbit(1., 1e-2, opt)

	} else {
		panic("Unknown tracer configuration! (`gen` can be `stream` or `circ`)")
	}
}

// test to do (actually)
// analytic potential
// with an orbit we believe

// change values on force softening scale
// pericenter > softening length
