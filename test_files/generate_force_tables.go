package main

import (
	"fmt"
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

func main() {
	thetas := []float64{-1, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0}
	criterias := []gravitree.OpeningCriteria{
		gravitree.BarnesHut, gravitree.PKDGRAV3,
		gravitree.SalmonWarren,
	}
	nStrs := []string{"4"}

	for in := range nStrs {
		filename := fmt.Sprintf("einasto_n=%s_a=18.dat", nStrs[in])
		x := readPointFile(filename)
		acc := gravitree.Acceleration(make([][3]float64, len(x)))

		for ic := range criterias {
			for it := range thetas {
				for i := range acc {
					acc[i] = [3]float64{}
				}

				if it == 0 {
					gravitree.BruteForceAcceleration(0.0, x, acc)
				} else {
					opt := gravitree.TreeOptions{}
					opt.Criteria = criterias[ic]
					opt.Theta = thetas[it]

					tree := gravitree.NewTree(x, opt)
					tree.Evaluate(0.0, acc)
				}
				filename = fmt.Sprintf("force_table_n=%s_ic=%d_it=%d.dat",
					nStrs[in], ic, it)
				writeAcceleration(filename, acc)
			}
		}
	}
}
