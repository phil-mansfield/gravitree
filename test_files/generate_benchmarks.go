package main

import (
	"fmt"
	"os"

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

func writeBenchmarks(filename, comment string, times [][]float64) {
	// (The comment should be a a string with '#' in front of each line which explains
	// the file contents.)
	f, err := os.Create(filename)
	if err != nil {
		panic(err.Error())
	}
	fmt.Fprintf(f, "%s\n", comment)
	for i := range times {
		for j := range times[i] {
			fmt.Fprintf(f, "%.10g ", acc[i][0], acc[i][1], acc[i][2])
		}
		fmt.Fprintf(f, "\n")
	}
}

func main() {
	// TODO: Write something similar to the main funciton of generate_force_tables.go
	// which instead runs benchmarks and writes corresponding .dat files to test_files/.
	// I think the idea would be something like  this - make a nsted for loop over all
	// the meta parameters you want to test (e.g., opening angle, opening criteria, leaf
	// node size, whatever), then each combination generates a new .dat file. Inside the
	// .dat file is a table where each row is a different particle count and each column
	// is a different function you're benchmarking. Then you make Python scritps which
	// read in these tables and display their outputs.
	//
	// One set of tests that I think is borderline essential is one where we put a shell
	// of test points at different radii and compute the time to compute potentials/
	// forces/whatever at that radius using brute force (i.e. just a double for loop
	// over the the points in the shell, then all the points in the halo) versus
	// a tree code calculation. We know at very small particle counts, brute force should
	// win and at very large counts the opposite should be true. We know that how much
	// speed up we get from the tree should degrade as the opening angle decreases.
	panic("NYI")

	// metaparameters
	// thetas := []float64{-1, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0}
	// criterias := []gravitree.OpeningCriteria{
	// 	gravitree.BarnesHut, gravitree.PKDGRAV3,
	// 	gravitree.SalmonWarren,
	// }
	// nStrs := []string{"4"}

	// for in := range nStrs {
	// 	filename := fmt.Sprintf("einasto_n=%s_a=18.dat", nStrs[in])
	// 	x := readPointFile(filename)
	// 	acc := gravitree.Acceleration(make([][3]float64, len(x)))

	// 	for ic := range criterias {
	// 		for it := range thetas {
	// 			for i := range acc {
	// 				acc[i] = [3]float64{}
	// 			}

	// 			if it == 0 {
	// 				gravitree.BruteForceAcceleration(0.0, x, acc)
	// 			} else {
	// 				opt := gravitree.TreeOptions{}
	// 				opt.Criteria = criterias[ic]
	// 				opt.Theta = thetas[it]

	// 				tree := gravitree.NewTree(x, opt)
	// 				tree.Evaluate(0.0, acc)
	// 			}
	// 			filename = fmt.Sprintf("force_table_n=%s_ic=%d_it=%d.dat",
	// 				nStrs[in], ic, it)
	// 			writeAcceleration(filename, acc)
	// 		}
	// 	}
	// }
}
