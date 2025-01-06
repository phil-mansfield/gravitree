package main

import (
	"fmt"
	"os"

	"github.com/phil-mansfield/symtable"
	"github.com/phil-mansfield/gravitree"
)

func readPointFile(filename string) [][3]float64{

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

func writeAcceleration(filename string, acc [][3]float64) {
	f, err := os.Create(filename)
	if err != nil { panic(err.Error()) }
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
	nStrs := []string{ "4" }

	for in := range nStrs {
		filename := fmt.Sprintf("einasto_n=%s_a=18.dat", nStrs[in])
		x := readPointFile(filename)
		acc := gravitree.Acceleration(make([][3]float64, len(x)))

		for ic := range criterias {
			for it := range thetas {
				for i := range acc { acc[i] = [3]float64{ } }

				if it == 0 {
					gravitree.BruteForceAcceleration(0.0, x, acc)
				} else {
					opt := gravitree.TreeOptions{ }
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