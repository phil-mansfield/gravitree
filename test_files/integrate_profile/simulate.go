package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"

	"github.com/phil-mansfield/gravitree"
	"github.com/phil-mansfield/gravitree/utils"
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

func SimulateEinasto(
	pos [][3]float64,
	vel [][3]float64,
	opt gravitree.SimulationOptions,
	prof *gravitree.Einasto,
) {

	orbitFile := filepath.Join(opt.SaveDirectory, "t_%d.dat")
	energyFile := filepath.Join(opt.SaveDirectory, "energy.dat")
	ngFile := filepath.Join(opt.SaveDirectory, "ng.dat")
	timeFile := filepath.Join(opt.SaveDirectory, "time.dat")
	accFile := filepath.Join(opt.SaveDirectory, "acc.dat")

	ok := make([]bool, len(pos))
	for i := range ok {
		ok[i] = true
	}

	t := 0.0

	_, err := os.Create(energyFile)
	if err != nil {
		panic(err.Error())
	}

	_, err = os.Create(timeFile)
	if err != nil {
		panic(err.Error())
	}

	endTime := float64(opt.Steps) * opt.Dt

	for i := range opt.Steps {

		accNorm := 0.0

		if i%opt.SaveEvery == 0 {

			// fmt.Printf("vc: %v, pot: %v, mr: %v \n", prof.GetCircularVelocity(1.), prof.GetPotential(1.), prof.EnclosedMass(1.))
			fmt.Printf("(einasto) step %v \n", i)

			// energy
			e := prof.GetEnergy(pos, vel)

			// acceleration
			acc := make([][3]float64, len(pos))
			prof.GetAcceleration(pos, acc)
			a := make([]float64, len(pos))
			for k := 0; k < len(pos); k++ {
				a[k] = utils.GetNorm(acc[k])
			}
			ti := []float64{t}

			// orbit parameters
			writeOrbit(fmt.Sprintf(orbitFile, i/opt.SaveEvery), pos, vel, ok)

			//
			ng := make([]float64, len(pos))
			for k := 0; k < len(pos); k++ {
				ng[k] = utils.GetNorm(utils.Cross(pos[k], vel[k]))
			}
			writeQuantity(ngFile, ng)

			writeQuantity(energyFile, e)
			writeQuantity(accFile, a)
			writeQuantity(timeFile, ti)
		}
		if opt.Integrator == "leapfrog" || opt.Integrator == "lfadp" {
			prof.LeapfrogStep(pos, vel, opt.Dt)
			if opt.Integrator == "lfadp" {
				// acc := gravitree.Acceleration(make([][3]float64, 1))
				// gravitree.BruteForceAccelerationAt(opt.Eps, tree.Points, pos, acc)

				acc := make([][3]float64, len(pos))
				prof.GetAcceleration(pos, acc)

				for k := 0; k < len(acc); k++ {
					accNorm = math.Max(accNorm, utils.GetNorm(acc[k]))
				}

				opt.Dt = gravitree.AdaptTimeStep(opt.Eps, accNorm, 0.025)
			}
		} else if opt.Integrator == "rk4" {
			prof.RK4Step(pos, vel, opt.Dt)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt

		if t >= endTime {
			break
		}
	}
}

func SimulatePlummer(
	pos [][3]float64,
	vel [][3]float64,
	opt gravitree.SimulationOptions,
	prof *gravitree.Plummer,
) {

	orbitFile := filepath.Join(opt.SaveDirectory, "t_%d.dat")
	energyFile := filepath.Join(opt.SaveDirectory, "energy.dat")
	timeFile := filepath.Join(opt.SaveDirectory, "time.dat")
	accFile := filepath.Join(opt.SaveDirectory, "acc.dat")

	ok := make([]bool, len(pos))

	for i := range ok {
		ok[i] = true
	}

	t := 0.0

	_, err := os.Create(energyFile)
	if err != nil {
		panic(err.Error())
	}

	_, err = os.Create(timeFile)
	if err != nil {
		panic(err.Error())
	}

	endTime := float64(opt.Steps) * opt.Dt

	for i := range opt.Steps {

		accNorm := 0.0

		if i%opt.SaveEvery == 0 {
			fmt.Printf("(plummer) step %v \n", i)

			// energy
			e := prof.GetEnergy(pos, vel)

			// acceleration
			acc := make([][3]float64, len(pos))
			prof.GetAcceleration(pos, acc)
			a := make([]float64, len(pos))
			for k := 0; k < len(pos); k++ {
				a[k] = utils.GetNorm(acc[k])
			}

			// time
			ti := []float64{t}

			// orbit parameters
			writeOrbit(fmt.Sprintf(orbitFile, i/opt.SaveEvery), pos, vel, ok)

			writeQuantity(energyFile, e)
			writeQuantity(accFile, a)
			writeQuantity(timeFile, ti)
		}

		if opt.Integrator == "leapfrog" || opt.Integrator == "lfadp" {
			prof.LeapfrogStep(pos, vel, opt.Dt)
			// Might be slow.
			if opt.Integrator == "lfadp" {
				// acc := gravitree.Acceleration(make([][3]float64, 1))
				// gravitree.BruteForceAccelerationAt(opt.Eps, tree.Points, pos, acc)

				acc := make([][3]float64, len(pos))
				prof.GetAcceleration(pos, acc)

				for k := 0; k < len(acc); k++ {
					accNorm = math.Max(accNorm, utils.GetNorm(acc[k]))
				}

				opt.Dt = gravitree.AdaptTimeStep(opt.Eps, accNorm, 0.025)

				if opt.Dt == 0 {
					panic("Dt is zero.")
				}
			}
		} else if opt.Integrator == "rk4" {
			prof.RK4Step(pos, vel, opt.Dt)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt

		if t >= endTime {
			break
		}

	}
}

func SimulatePoint(
	pos [][3]float64,
	vel [][3]float64,
	opt gravitree.SimulationOptions,
	prof *gravitree.PointMass,
) {

	orbitFile := filepath.Join(opt.SaveDirectory, "t_%d.dat")
	energyFile := filepath.Join(opt.SaveDirectory, "energy.dat")
	timeFile := filepath.Join(opt.SaveDirectory, "time.dat")
	accFile := filepath.Join(opt.SaveDirectory, "acc.dat")

	ok := make([]bool, len(pos))

	for i := range ok {
		ok[i] = true
	}

	t := 0.0

	_, err := os.Create(energyFile)
	if err != nil {
		panic(err.Error())
	}

	_, err = os.Create(timeFile)
	if err != nil {
		panic(err.Error())
	}

	endTime := float64(opt.Steps) * opt.Dt

	for i := range opt.Steps {

		accNorm := 0.0

		if i%opt.SaveEvery == 0 {
			fmt.Printf("(point) step %v \n", i)

			// energy
			e := prof.GetEnergy(pos, vel)

			// acceleration
			acc := make([][3]float64, len(pos))
			prof.GetAcceleration(pos, acc)
			a := make([]float64, len(pos))
			for k := 0; k < len(pos); k++ {
				a[k] = utils.GetNorm(acc[k])
			}
			ti := []float64{t}

			// orbit parameters
			writeOrbit(fmt.Sprintf(orbitFile, i/opt.SaveEvery), pos, vel, ok)

			writeQuantity(energyFile, e)
			writeQuantity(accFile, a)
			writeQuantity(timeFile, ti)
		}
		if opt.Integrator == "leapfrog" || opt.Integrator == "lfadp" {
			prof.LeapfrogStep(pos, vel, opt.Dt)
			if opt.Integrator == "lfadp" {
				// acc := gravitree.Acceleration(make([][3]float64, 1))
				// gravitree.BruteForceAccelerationAt(opt.Eps, tree.Points, pos, acc)

				acc := make([][3]float64, len(pos))
				prof.GetAcceleration(pos, acc)

				for k := 0; k < len(acc); k++ {
					accNorm = math.Max(accNorm, utils.GetNorm(acc[k]))
				}

				opt.Dt = gravitree.AdaptTimeStep(opt.Eps, accNorm, 0.025)
			}
		} else if opt.Integrator == "rk4" {
			prof.RK4Step(pos, vel, opt.Dt)
		} else {
			panic("Unknown integrator specified.")
		}

		t += opt.Dt

		if t >= endTime {
			break
		}
	}
}

func main() {

	gen := flag.String("gen", "einasto", "sim configuration")
	integrator := flag.String("int", "leapfrog", "integrator")
	dt := flag.String("dt", "1e-4", "timestep")
	// icFile := flag.String("ic", "ic.dat", "initial phase space file")
	vcFrac := flag.Float64("vc", 1., "fraction of circular speed")
	saveDir := flag.String("saveto", "", "directory to store snapshots")
	adaptive := flag.Bool("adp", false, "adaptive timestepping for leapfrog")

	flag.Parse()

	path := "snapshots"

	dtVal, _ := strconv.ParseFloat(*dt, 64)
	// ic := gravitree.ReadPhaseSpaceFile(*icFile)
	pos := [][3]float64{{1., 0., 0.}}
	r0 := utils.GetNorm(pos[0])

	time := 3 * math.Pi
	steps := int(time / dtVal)

	// steps := 1000

	integ_type := *integrator

	if *adaptive {
		integ_type = "lfadp"
	}

	eps := 1e-5

	if *gen == "einasto" {
		path = fmt.Sprintf("einasto_int=%s_dt=%s", *integrator, *dt)
		path = filepath.Join(*saveDir, path)
		MakeDirectory(path)

		opt := gravitree.SimulationOptions{
			Steps: steps, Dt: dtVal, SaveEvery: 10,
			SaveDirectory: path, Integrator: integ_type,
			Eps: eps,
		}

		prof := gravitree.Einasto{
			Alpha: 0.18,
			Rs:    0.1,
		}

		vel := [][3]float64{{0., (*vcFrac) * prof.GetCircularVelocity(r0), 0.}}
		SimulateEinasto(pos, vel, opt, &prof)

	} else if *gen == "plummer" {
		path = fmt.Sprintf("plummer_int=%s_dt=%s", integ_type, *dt)
		path = filepath.Join(*saveDir, path)
		MakeDirectory(path)

		opt := gravitree.SimulationOptions{
			Steps: steps, Dt: dtVal, SaveEvery: 10,
			SaveDirectory: path, Integrator: integ_type,
			Eps: eps,
		}

		prof := gravitree.Plummer{
			B: .05,
		}
		vel := [][3]float64{{0., (*vcFrac) * prof.GetCircularVelocity(r0), 0.}}
		SimulatePlummer(pos, vel, opt, &prof)
	} else if *gen == "point" {
		path = fmt.Sprintf("point_int=%s_dt=%s", integ_type, *dt)
		path = filepath.Join(*saveDir, path)
		MakeDirectory(path)

		opt := gravitree.SimulationOptions{
			Steps: steps, Dt: dtVal, SaveEvery: 10,
			SaveDirectory: path, Integrator: integ_type,
			Eps: eps,
		}

		prof := gravitree.PointMass{
			Mass: 1.,
		}

		vel := [][3]float64{{0., (*vcFrac) * prof.GetCircularVelocity(r0), 0.}}
		SimulatePoint(pos, vel, opt, &prof)
	} else {
		panic("Unknown tracer configuration! (`gen` can be `einasto`, `plummer`, or `point`)")
	}
}
