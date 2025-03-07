package main

import (
	"fmt"
	"math"
	"os"

	"github.com/phil-mansfield/gravitree"
)

func writeAcceleration(filename string, acc [][3]float64) {
	f, err := os.Create(filename)
	if err != nil {
		panic(err.Error())
	}
	for i := range acc {
		fmt.Fprintf(f, "%.10g %.10g %.10g\n", acc[i][0], acc[i][1], acc[i][2])
	}
}

func writePotential(filename string, phi []float64) {
	f, err := os.Create(filename)
	if err != nil {
		panic(err.Error())
	}
	for i := range phi {
		fmt.Fprintf(f, "%.10g\n", phi[i])
	}
}

func main() {

	thetas := []float64{-1, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0}
	criterias := []gravitree.OpeningCriteria{
		gravitree.BarnesHut, gravitree.PKDGRAV3,
		gravitree.SalmonWarren,
	}

	// number of particles in the halo tree
	// to compute quantities over
	nStrs := []string{"4"}

	for in := range nStrs {
		t1_fn := fmt.Sprintf("einasto_n=%s_a=18.dat", nStrs[in])
		// t2_fn := "shell_n=5_r=20.dat"
		t2_fn := "einasto_n=2_a=18.dat"

		x1 := gravitree.ReadPointFile(t1_fn)
		x2 := gravitree.ReadPointFile(t2_fn)

		phi := gravitree.Potential(make([]float64, len(x2)))
		acc := gravitree.Acceleration(make([][3]float64, len(x2)))

		for ic := range criterias {
			for it := range thetas {
				for i := range x2 {
					phi[i] = 0.0
					acc[i] = [3]float64{}
				}

				if it == 0 {
					gravitree.BruteForcePotentialAt(0.0, x1, x2, phi)
					gravitree.BruteForceAccelerationAt(0.0, x1, x2, acc)
				} else {
					opt1 := gravitree.TreeOptions{}
					opt1.Criteria = criterias[ic]
					opt1.Theta = thetas[it]

					opt2 := gravitree.TreeOptions{LeafSize: 1}
					opt2.Criteria = criterias[ic]
					opt2.Theta = thetas[it]

					t1 := gravitree.NewTree(x1, opt1)
					t2 := gravitree.NewTree(x2, opt2)

					t1.EvaluateAt(t2, 0.0, acc)
					t1.EvaluateAt(t2, 0.0, phi)
				}

				phi_fn := fmt.Sprintf("pot_table_at_n=%s_ic=%d_it=%d.dat", nStrs[in], ic, it)
				acc_fn := fmt.Sprintf("force_table_at_n=%s_ic=%d_it=%d.dat", nStrs[in], ic, it)

				writePotential(phi_fn, phi)
				writeAcceleration(acc_fn, acc)
			}
		}

	}

	// AccelerationEvaluateAt("einasto_n=3.5_a=18.dat", "shell_n=5_r=20.dat", 1e-3)
	// AccelerationEvaluateAt("einasto_n=2_a=18.dat", "shell_n=10_r=20.dat", 1e-3)

	// PotentialEvaluateAt("einasto_n=3.5_a=18.dat", "shell_n=5_r=20.dat", 1e-3)
	// PotentialEvaluateAt("einasto_n=2_a=18.dat", "shell_n=10_r=20.dat", 1e-3)
}

// ignore everything below this line for now
// gonna remove it in the future but I'd like it for reference]
// - jay

func PotentialEvaluateAt(t1_fn, t2_fn string, eps float64) {

	x1 := readPointFile(t1_fn)
	x2 := readPointFile(t2_fn)

	t1 := gravitree.NewTree(x1)
	t2 := gravitree.NewTree(x2, gravitree.TreeOptions{LeafSize: 1})

	eps2 := eps * eps

	phi := gravitree.Potential(make([]float64, len(x2)))
	phi_bf := gravitree.Potential(make([]float64, len(x2)))

	t1.EvaluateAt(t2, eps, phi)
	gravitree.BruteForcePotentialAt(eps2, x1, x2, phi_bf)

	if !multArrayAlmostEq(phi, 1., phi_bf, eps) {
		fmt.Printf("expected phi = %.4f, but phi = %.4f", phi_bf, phi)
	} else {
		fmt.Printf("PASS")
	}
}

func AccelerationEvaluateAt(t1_fn, t2_fn string, eps float64) {
	x1 := readPointFile(t1_fn)
	x2 := readPointFile(t2_fn)

	t1 := gravitree.NewTree(x1)
	t2 := gravitree.NewTree(x2, gravitree.TreeOptions{LeafSize: 1})

	eps2 := eps * eps

	acc := gravitree.Acceleration(make([][3]float64, len(x2)))
	acc_bf := gravitree.Acceleration(make([][3]float64, len(x2)))

	t1.EvaluateAt(t2, eps, acc)
	gravitree.BruteForceAccelerationAt(eps2, x1, x2, acc_bf)

	for k := range x2 {
		acc_pt := acc[k]
		acc_bf_pt := acc_bf[k]

		var acc_mag, acc_mag_bf float64
		for i := 0; i < 3; i++ {
			acc_mag += acc_pt[i] * acc_pt[i]
			acc_mag_bf += acc_bf_pt[i] * acc_bf_pt[i]
		}
		acc_mag = math.Sqrt(acc_mag)
		acc_mag_bf = math.Sqrt(acc_mag_bf)

		// check if it matches
		if !almostEq(acc_mag, acc_mag_bf, 1e-3) {
			fmt.Printf("(%d.0) expected acc = %.4f, got acc = %.4f \n", k, acc_mag_bf, acc_mag)
		}
	}

}

func almostEq(x, y, eps float64) bool {
	return x+eps > y && x-eps < y
}

func multArrayAlmostEq(
	x []float64, mult float64, y []float64, eps float64,
) bool {
	if len(x) != len(y) {
		return false
	}

	for i := range x {
		if !almostEq(x[i], mult*y[i], eps) {
			return false
		}
	}
	return true
}
