package gravitree

import (
	"fmt"
	"math"
)

// Potential computes the potential for each point in the tree.
func (t *Tree) Potential(eps float64, phi []float64) {
	if len(phi) != len(t.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(phi) = %d",
			len(t.Nodes), len(phi)))
	}

	t.eps2 = eps*eps
	
	for i := range t.Nodes {
		// Only compute the potential at the leaf nodes.
		if t.Nodes[i].Left == -1 {
			t.walkNodePotential(i, 0, phi)
		}
	}
}

// useMonopole returns true if node i is close enough to node j that a monopole
// approximation can be used.
func (t *Tree) useMonopole(i, j int) bool {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	r2 := dist2(nodei.Center, nodej.Center)
	return r2 > nodei.RMax2 + nodej.ROpen2 || r2 < t.eps2
}

// i is the index of the node that phi corresponds to, j is the index of the
// node that is being walked.
func (t *Tree) walkNodePotential(i, j int, phi []float64) {
	target := &t.Nodes[j]
	
	if i == j {
		t.pairwisePotential(i, phi)
	} else if t.useMonopole(i, j) {
		t.monopolePotential(i, j, phi)
	} else if target.Left == -1 {
		t.oneSidedPotential(i, j, phi)
	} else {
		t.walkNodePotential(i, target.Left, phi)
		t.walkNodePotential(i, target.Right, phi)
	}
}

// pairwisePotential computes the potential at every point in node i using the
// other points in node i.
func (t *Tree) pairwisePotential(i int, phi []float64) {
	node := &t.Nodes[i]
	for i := node.Start; i < node.End; i++ {
		xi, idxi := t.Points[i], t.Index[i]
		for j := i+1; j < node.End; j++ {
			xj, idxj := t.Points[j], t.Index[j]

			phiij := pointPotential(dist2(xi, xj), t.eps2) 
			phi[idxi] += phiij
			phi[idxj] += phiij
		}
	}
}

// monopolePotential computes the potential 
func (t *Tree) monopolePotential(i, j int, phi []float64) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	xj := nodej.Center
	massj := float64(nodej.End - nodej.Start)
	
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := t.Points[i], t.Index[i]
		phi[idxi] += pointPotential(dist2(xi, xj), t.eps2)*massj
	}
}

// oneSidedPotential computes the potential at every point in i using every
// point in j.
func (t *Tree) oneSidedPotential(i, j int, phi []float64) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]

	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := t.Points[i], t.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := t.Points[j]

			phi[idxi] += pointPotential(dist2(xi, xj), t.eps2)
		}
	}
}

func dist2(x1, x2 [3]float64) float64 {
	dx2 := 0.0
	for k := 0; k < 3; k++ {
		dx := x2[k] - x1[k]
		dx2 += dx*dx
	}
	return dx2
}

func pointPotential(r2, eps2 float64) float64 {
	return -1.0 / math.Sqrt(r2 + eps2)
}
