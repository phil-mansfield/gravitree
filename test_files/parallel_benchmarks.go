package main

import (
	"fmt"
	"time"
	"runtime"
	
	"github.com/phil-mansfield/gravitree"
	"github.com/phil-mansfield/symtable"
)

func main() {
	files := []string{
		"einasto_n=2_a=18.dat",
		"einasto_n=2.5_a=18.dat",
		"einasto_n=3_a=18.dat",
		"einasto_n=3.5_a=18.dat",
		"einasto_n=4_a=18.dat",
		"einasto_n=4.5_a=18.dat",
		"einasto_n=5_a=18.dat",
	}

	exp := []string{
		"2", "2.5", "3", "3.5", "4", "4.5", "5",
	}
	trials := []int{ 100, 100, 30, 10, 10, 10, 10 }
	threads := []int{ 1, 2, 4, 8, 16, 32, 64, 128, 256 }

	fmt.Printf("Cores: %d\n", runtime.GOMAXPROCS(-1))
	
	for k := range files {
		dt := make([]float64, len(threads))
		
		x := readPointFile(files[k])
		tree := gravitree.NewTree(x)
		phi := gravitree.Potential(make([]float64, len(x)))
		
		for i := range threads {
			gravitree.SetThreads(threads[i])
			
			for j := 0; j < trials[k]; j++ {
			
				t0 := time.Now()
				
				tree.Evaluate(0.001, phi)
				
				t1 := time.Now()
				dt[i] += t1.Sub(t0).Seconds() / float64(len(x))
			}

			dt[i] /= float64(trials[k])
			dt[i] *= 1e6
		}

		fmt.Printf("np = 10^%s\n", exp[k])
		fmt.Printf("Threads:                 %7d\n", threads)
		fmt.Printf("wall time/particle (µs): %7.3f\n", dt)
		for i := range dt { dt[i] *= float64(threads[i]) }
		fmt.Printf("cpu time/particle (µs):  %7.3f\n", dt)
	}
}

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
