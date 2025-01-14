package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand/v2"
	"os"

	"github.com/phil-mansfield/gravitree"
	"github.com/phil-mansfield/symtable"
)

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
	w_acc *gravitree.Acceleration,
	eps float64,
) {

	// calculate current acceleration
	tree.EvaluateAt(&w_tree.Tree, eps, w_acc)

	// perform half kick
	for i := range pos {

		if ok[i] {
			for j := range 3 {
				vel[i][j] += (dt / 2.0) * (*w_acc)[i][j]
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
				vel[i][j] += (dt / 2.0) * (*w_acc)[i][j]
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

// func vecNorm(x [3]float64) float64 {
// 	r := 0.0
// 	for i := 0; i < 3; i++ {
// 		r += x[i] * x[i]
// 	}
// 	return math.Sqrt(r)
// }

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

	w_acc := gravitree.Acceleration(make([][3]float64, len(pos)))
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
			particle_fn := "../einasto_n=3_a=18.dat"
			halo_x := readPointFile(particle_fn)
			LeapfrogStepBF(pos, vel, halo_x, dt, ok, eps)

		} else {
			LeapfrogStep(pos, vel, tree, dt, ok, w_tree, &w_acc, eps)
		}

		t += dt
	}

}

func circ(
	npart int,
	r0 float64,
	dr float64,
	steps int,
	dt float64,
	ext_pos [][3]float64,
	snap_cadence int,
	brute_force bool,
	eps float64,
) {
	orbit_fn := "snapshots/circ_t_%d.dat"

	if brute_force {
		orbit_fn = "snapshots/circ_bf_t_%d.dat"
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

	w_acc := gravitree.Acceleration(make([][3]float64, len(pos)))
	w_tree := gravitree.NewArrayTree(pos)

	for i := range steps {
		fmt.Printf("(circ) step %v; brute: %t \n", i, brute_force)

		if int(t/dt)%snap_cadence == 0 {
			writeOrbit(fmt.Sprintf(orbit_fn, i/snap_cadence), pos, vel, ok)
		}

		if brute_force {
			LeapfrogStepBF(pos, vel, ext_pos, dt, ok, eps)

		} else {
			LeapfrogStep(pos, vel, tree, dt, ok, w_tree, &w_acc, eps)
		}

		t += dt
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

	flag.Parse()

	// check if snapshot folder exists

	_, err := os.Stat("snapshots")

	if os.IsNotExist(err) {
		err := os.MkdirAll("snapshots", os.ModePerm)
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
	particle_fn := "../einasto_n=3_a=18.dat"
	ext_pos := readPointFile(particle_fn)

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
		circ(
			10,
			30.,
			1.,
			int(2e6),
			dt,
			ext_pos,
			save_every,
			brute,
			eps,
		)
	} else {
		panic("Unknown tracer configuration! (`gen` can be `stream` or `circ`)")
	}
}
