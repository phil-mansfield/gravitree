package main

import (
	"os"
	"strconv"
	"fmt"

	"github.com/phil-mansfield/symtable"
	"github.com/phil-mansfield/gravitree"
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

func main() {
	filename := os.Args[1]
	theta, err := strconv.ParseFloat(os.Args[2], 64)
	if err != nil { panic(err.Error()) }

	x := readPointFile(filename)
	n := len(x)
	t := gravitree.NewTree(x, gravitree.TreeOptions{ Theta: theta })

	phi := gravitree.Potential(make([]float64, n))
	cs := gravitree.NewCallStructure(t, t, phi)
	t.Evaluate(0, cs)

	fmt.Printf("Naive force evaluations: %.3g\n", float64(n*(n-1)/2))
	fmt.Printf("Points: %.3g\n", float64(len(t.Points)))
	fmt.Printf("Nodes: %.3g\n", float64(len(t.Nodes)))

	nl := 0
	for i := range t.Nodes {
		if t.Nodes[i].Left == -1 { nl++ }
	}
	fmt.Printf("Leafs: %.3g\n", float64(nl))


	nc, nf, nftot := 0, 0, 0
	for i := range cs.ApproximateCalls {
		nc += len(cs.ApproximateCalls[i])

		for j := range cs.ApproximateCalls[i] {
			n := cs.Sizes2[cs.ApproximateCalls[i][j].I2]
			
			nf += n
		}		
	}
	nftot += nf

	fmt.Printf("Approximate() calls: %.3g force evaluations: %.3g\n",
		float64(nc), float64(nf))

	nc, nf = 0, 0
	for i := range cs.OneSidedLeafCalls {
		nc += len(cs.OneSidedLeafCalls[i])
		for j := range cs.OneSidedLeafCalls[i] {
			n1 := cs.Sizes1[cs.OneSidedLeafCalls[i][j].I1]
			n2 := cs.Sizes2[cs.OneSidedLeafCalls[i][j].I2]

			nf += n1*n2
		}
	}
	nftot += nf

	fmt.Printf("OneSidedLeaf() calls: %.3g force evaluations: %.3g\n",
		float64(nc), float64(nf))

	nc, nf = 0, 0


	nc, nf = 0, 0
	for i := range cs.TwoSidedLeafCalls {
		nc += len(cs.TwoSidedLeafCalls[i])
		for j := range cs.TwoSidedLeafCalls[i] {
			n := cs.Sizes1[cs.TwoSidedLeafCalls[i][j].I1]

			nf += n*(n-1)/2
		}
	}

	nftot += nf

	fmt.Printf("TwoSidedLeaf() calls:%.3g force evaluations: %.3g\n",
		float64(nc), float64(nf))

	fmt.Printf("Total tree force evaluations: %.3g\n", float64(nftot))
	fmt.Printf("Tree efficiency: %.3g\n",
		float64(nftot)/float64(len(x)*(len(x)-1)/2))
}