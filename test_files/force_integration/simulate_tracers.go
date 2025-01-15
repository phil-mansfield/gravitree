package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand/v2"
	"os"

	"path/filepath"

	"github.com/phil-mansfield/gravitree"
	"github.com/phil-mansfield/symtable"
)

type SimulationOptions struct {
	steps      int
	dt         float64
	save_every int
	eps        float64
	brute      bool
	save_dir   string
	snap_fn    string
}

func readPointFile(filename string) [][3]float64 {

	t := symtable.TextFile(filename)
	cols := t.ReadFloat64s([]int{0, 1, 2}) // column indices
	xs := cols[0]
	ys := cols[1]
	zs := cols[2]

	var result [][3]float64

	for i := 0; i < len(xs); i++ {
		var point [3]float64
		point[0] = xs[i]
		point[1] = ys[i]
		point[2] = zs[i]
		result = append(result, point)
	}

	return result
}

func rescalePoints(x [][3]float64, a float64) [][3]float64 {
	res := make([][3]float64, len(x))
	for i := range x {
		for k := 0; k < 3; k++ {
			res[i][k] = a * x[i][k]
		}
	}
	return res
}

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
	f, err := os.Create(filename)
	if err != nil {
		panic(err.Error())
	}
	for i := 0; i < len(q); i++ {
		fmt.Fprintf(f, "%.10g \n", q[i])

	}
}

func cart_to_sph(pos [3]float64) [3]float64 {
	res := [3]float64{0, 0, 0}
	x, y, z := pos[0], pos[1], pos[2]
	res[0] = math.Sqrt(x*x + y*y + z*z) // this is the radius
	res[1] = math.Acos(z / res[0])
	res[2] = math.Atan2(y, x)
	return res
}

func sph_to_cart(sph [3]float64) [3]float64 {
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

func LeapfrogStep(
	pos [][3]float64,
	vel [][3]float64,
	tree *gravitree.Tree,
	dt float64,
	ok []bool,
	w_tree *gravitree.ArrayTree,
	eps float64,
) {

	w_acc := gravitree.Acceleration(make([][3]float64, len(pos)))

	// calculate current acceleration
	w_tree.Tree.Points = pos
	w_tree.Update()
	tree.EvaluateAt(&w_tree.Tree, eps, w_acc)

	// perform half kick
	for i := range pos {

		if ok[i] {
			for j := range 3 {
				vel[i][j] += (dt / 2.0) * (w_acc)[i][j]
			}
		}
	}

	// perform drift with new velocity
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				pos[i][j] += dt * vel[i][j]
			}
		}
	}

	// rewrite tree and get new accelerations
	// Set tree points to new positions
	w_tree.Tree.Points = pos
	w_tree.Update()
	tree.EvaluateAt(&w_tree.Tree, eps, w_acc)

	// half kick
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				vel[i][j] += (dt / 2.0) * (w_acc)[i][j]
			}
		}
	}
}

func LeapfrogStepBF(
	pos [][3]float64,
	vel [][3]float64,
	halo_pos [][3]float64,
	dt float64,
	ok []bool,
	eps float64,
) {

	w_acc := gravitree.Acceleration(make([][3]float64, len(pos)))

	// calculate current acceleration
	gravitree.BruteForceAccelerationAt(eps, halo_pos, pos, w_acc)

	// perform half kick
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				vel[i][j] += (dt / 2.0) * w_acc[i][j]
			}
		}
	}

	// perform drift with new velocity
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				pos[i][j] += dt * vel[i][j]
			}
		}
	}

	// rewrite tree and get new accelerations
	gravitree.BruteForceAccelerationAt(eps, halo_pos, pos, w_acc)

	// half kick
	for i := range pos {
		if ok[i] {
			for j := range 3 {
				vel[i][j] += (dt / 2.0) * w_acc[i][j]
			}
		}
	}
}

func particleCount(x [][3]float64, r float64) int {

	count := 0

	for i := range x {
		r_i := 0.0

		for k := 0; k < 3; k++ {
			r_i += x[i][k] * x[i][k]
		}

		r_i = math.Sqrt(r_i)

		if r_i <= r {
			count++
		}
	}
	return count
}

func getCircularVelocity(x [][3]float64, r float64) float64 {
	n_part := particleCount(x, r)
	return math.Sqrt(1. * float64(n_part) / r)
}

func vecNorm(x [3]float64) float64 {
	r := 0.0
	for i := 0; i < 3; i++ {
		r += x[i] * x[i]
	}
	return math.Sqrt(r)
}

func stream(
	x0 [3]float64,
	v0 [3]float64,
	rt float64,
	steps int,
	dt float64,
	spawn_every int,
	ext_pos [][3]float64,
	snap_cadence int,
	brute_force bool,
	eps float64,
) {
	orbit_fn := "snapshots/stream_t_%d.dat"

	if brute_force {
		orbit_fn = "snapshots/stream_bf_t_%d.dat"
		fmt.Printf("Calculating using brute force. \n")
	}

	stream_count := int(steps/spawn_every)*2 + 1

	fmt.Printf("stream count %v \n", stream_count)

	pos := make([][3]float64, stream_count)
	vel := make([][3]float64, stream_count)
	ok := make([]bool, stream_count)

	tree := gravitree.NewTree(ext_pos)
	// fix this

	t := 0.0

	for i := 0; i < stream_count; i++ {
		ok[i] = false
		pos[i] = [3]float64{0, 0, 0}
		vel[i] = [3]float64{0, 0, 0}
	}

	for k := 1; k < stream_count; k++ {
		pos[k][0] += rand.NormFloat64() * 1000
		pos[k][1] += rand.NormFloat64() * 1000
		pos[k][2] += rand.NormFloat64() * 1000
	}

	ok[0] = true

	// set initial conditions of stream particle
	pos[0] = x0
	vel[0] = v0

	w_tree := gravitree.NewArrayTree(pos)

	for i := range steps {
		fmt.Printf("(stream) step %v; brute: %t \n", i, brute_force)

		if int(t/dt)%snap_cadence == 0 {
			writeOrbit(fmt.Sprintf(orbit_fn, i/snap_cadence), pos, vel, ok)
		}

		if int(t/dt)%spawn_every == 0 {
			ok[i/spawn_every+1] = true
			ok[i/spawn_every+2] = true

			// v_circ := getCircularVelocity(ext_pos, vecNorm(pos[0]))
			spray := getSprayPair(pos[0], vel[0], rt)

			pos[i/spawn_every+1] = spray[0]
			vel[i/spawn_every+1] = spray[1]

			pos[i/spawn_every+2] = spray[2]
			vel[i/spawn_every+2] = spray[3]
		}

		if brute_force {
			LeapfrogStepBF(pos, vel, ext_pos, dt, ok, eps)

		} else {
			LeapfrogStep(pos, vel, tree, dt, ok, w_tree, eps)
		}

		t += dt
	}

}

func calculateEnergy(
	pos [][3]float64,
	vel [][3]float64,
	halo_pos [][3]float64,
	ok []bool,
	eps float64,
	brute bool,
) []float64 {

	res := make([]float64, len(pos))
	pot := gravitree.Potential(make([]float64, len(pos)))
	// potential
	if !brute {
		tree := gravitree.NewTree(halo_pos)
		w_tree := gravitree.NewArrayTree(pos)
		w_tree.SetEvaluateFlag(ok)
		w_tree.Update()
		tree.EvaluateAt(&w_tree.Tree, eps, pot)
	} else {
		gravitree.BruteForcePotentialAt(eps, halo_pos, pos, pot)
	}

	for i := range pos {
		// kinetic
		if ok[i] {
			res[i] += vecNorm(vel[i])
			res[i] *= res[i] / 2.
			res[i] += pot[i]
		}
	}

	return res
}

func single_particle(
	steps int,
	dt float64,
	snap_cadence int,
	brute_force bool,
	eps float64,
) {
	orbit_fn := "snapshots/single_t_%d.dat"
	energy_fn := "snapshots/single_e_%d.dat"

	if brute_force {
		orbit_fn = "snapshots/single_bf_t_%d.dat"
		energy_fn = "snapshots/single_bf_e_%d.dat"
		fmt.Printf("Calculating using brute force. \n")
	}

	pos := make([][3]float64, 1)
	vel := make([][3]float64, 1)
	ok := make([]bool, 1)

	ext_pos := [][3]float64{
		{0., 0., 0.},
	}

	tree := gravitree.NewTree(ext_pos)
	// fix this

	t := 0.0

	ok[0] = true
	pos[0] = [3]float64{.1, 0, 0}
	vel[0] = [3]float64{0, 1. / math.Sqrt(.1), 0}

	w_tree := gravitree.NewArrayTree(pos)

	for i := range steps {
		fmt.Printf("(single) step %v; brute: %t \n", i, brute_force)

		if int(t/dt)%snap_cadence == 0 {
			writeOrbit(fmt.Sprintf(orbit_fn, i/snap_cadence), pos, vel, ok)
			e := calculateEnergy(
				pos,
				vel,
				ext_pos,
				ok,
				eps,
				brute_force,
			)
			writeQuantity(fmt.Sprintf(energy_fn, i/snap_cadence), e)
		}

		if brute_force {
			LeapfrogStepBF(pos, vel, ext_pos, dt, ok, eps)

		} else {
			LeapfrogStep(pos, vel, tree, dt, ok, w_tree, eps)

		}

		t += dt
	}

}

func circ(
	npart int,
	r0 float64,
	dr float64,
	ext_pos [][3]float64,
	opt SimulationOptions,
) {

	orbit_fn := filepath.Join(opt.save_dir, "circ_t_%d.dat")
	energy_fn := filepath.Join(opt.save_dir, "circ_e_%d.dat")

	if opt.brute {
		orbit_fn = filepath.Join(opt.save_dir, "circ_bf_t_%d.dat")  //
		energy_fn = filepath.Join(opt.save_dir, "circ_bf_e_%d.dat") //
		fmt.Printf("Calculating brute force... \n")
	}

	pos := make([][3]float64, npart)
	vel := make([][3]float64, npart)
	ok := make([]bool, npart)

	tree := gravitree.NewTree(ext_pos)

	t := 0.0

	for i := 0; i < npart; i++ {
		ok[i] = true
		pos[i] = [3]float64{r0 + float64(i)*dr, 0, 0}
		vel[i] = [3]float64{0, getCircularVelocity(ext_pos, r0+float64(i)*dr), 0}
	}

	w_tree := gravitree.NewArrayTree(pos)

	for i := range opt.steps {
		fmt.Printf("(circ) step %v; brute: %t \n", i, opt.brute)

		if int(t/opt.dt)%opt.save_every == 0 {
			writeOrbit(fmt.Sprintf(orbit_fn, i/opt.save_every), pos, vel, ok)
			e := calculateEnergy(
				pos,
				vel,
				ext_pos,
				ok,
				opt.eps,
				opt.brute,
			)
			writeQuantity(fmt.Sprintf(energy_fn, i/opt.save_every), e)
		}

		if opt.brute {
			LeapfrogStepBF(pos, vel, ext_pos, opt.dt, ok, opt.eps)

		} else {
			LeapfrogStep(pos, vel, tree, opt.dt, ok, w_tree, opt.eps)
		}

		t += opt.dt
	}

}

func getSprayPair(x0, v0 [3]float64, rt float64) [][3]float64 {
	q_sph := cart_to_sph(x0)

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

	q1 = sph_to_cart(q1)
	q2 = sph_to_cart(q2)

	return [][3]float64{q1, p1, q2, p2}
}

func main() {

	method := flag.String("method", "approx", "Brute force or tree approx.")
	gen := flag.String("gen", "stream", "generate a stream or Gaussian cloud")
	npts := flag.String("npts", "3", "einasto profile")

	flag.Parse()

	// check if snapshot folder exists

	save_dir := fmt.Sprintf("snapshots_n=%s", *npts)

	_, err := os.Stat(save_dir)

	if os.IsNotExist(err) {
		err := os.MkdirAll(save_dir, os.ModePerm)
		if err != nil {
			panic("Unable to create snapshot folder.")
		}
	}

	brute := false

	if *method == "brute" {
		brute = true
	}

	steps := int(1e6)
	// obtain static distribution to orbit around
	particle_fn := fmt.Sprintf("../einasto_n=%s_a=18.dat", *npts)
	ext_pos := readPointFile(particle_fn)
	ext_pos = rescalePoints(ext_pos, 1./(0.1))

	dt := .0001

	save_every := 1000
	spawn_every := 10000

	x0 := [3]float64{40., 0, 0}
	v0 := [3]float64{0, getCircularVelocity(ext_pos, 40.), 0}

	// fix this
	eps := 0.096 / math.Pow(float64(len(ext_pos)), 1./3.)

	if *gen == "stream" {
		stream(
			x0,
			v0,
			1e-2,
			steps,
			dt,
			spawn_every,
			ext_pos,
			save_every,
			brute,
			eps,
		)
	} else if *gen == "circ" {

		steps := int(5e5)

		opt := SimulationOptions{
			brute:      brute,
			steps:      steps,
			dt:         1e-6,
			save_every: int(steps / 1000),
			eps:        eps,
			save_dir:   save_dir,
		}

		circ(
			10,
			.5,
			.1,
			ext_pos,
			opt,
		)
	} else if *gen == "sing" {
		single_particle(
			int(1e4),
			dt,
			100,
			brute,
			eps,
		)
	} else {
		panic("Unknown tracer configuration! (`gen` can be `stream` or `circ`)")
	}
}
