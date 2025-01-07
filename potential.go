package gravitree

import (
	"fmt"
	"math"
)

type Potential []float64

var _ Quantity = Potential{}

func (phi Potential) Len() int { return len(phi) }

// pairwisePotential computes the potential at every point in node i using the
// other points in node i.
func (phi Potential) TwoSidedLeaf(t *Tree, i int) {
	node := &t.Nodes[i]
	for i := node.Start; i < node.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := i + 1; j < node.End; j++ {
			xj, idxj := &t.Points[j], t.Index[j]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]
			dx2 := dx*dx + dy*dy + dz*dz

			phiij := pointPotential(dx2, t.eps2)
			phi[idxi] += phiij
			phi[idxj] += phiij
		}
	}
}

func (phi Potential) Approximate(t1, t2 *Tree, i, j int) {
	// writes the approximated potential for nodes in
	// t2 from nodes in t1
	switch t2.Order {
	case Monopole:
		// Loops over the t2 nodes, calculating the
		// contributions from the t1 nodes.
		nodei, nodej := &t2.Nodes[i], &t1.Nodes[j]
		xj := &nodej.Center
		massj := float64(nodej.End - nodej.Start)

		for i := nodei.Start; i < nodei.End; i++ {
			xi, idxi := &t2.Points[i], t2.Index[i]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]
			dx2 := dx*dx + dy*dy + dz*dz
			phi[idxi] += pointPotential(dx2, t2.eps2) * massj
		}
	case Quadrupole:
		// TODO: Implement (Plummer) quadrupole appoximation of potential
		panic("NYI")
	default:
		panic(fmt.Sprintf("Unrecognized approximation order code, %d", t2.Order))
	}
}

func (phi Potential) OneSidedLeaf(t1, t2 *Tree, i, j int) {
	nodei, nodej := &t2.Nodes[i], &t1.Nodes[j]

	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t2.Points[i], t2.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := &t1.Points[j]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]
			dx2 := dx*dx + dy*dy + dz*dz
			phi[idxi] += pointPotential(dx2, t2.eps2)
		}
	}
}

func BruteForcePotential(eps float64, x [][3]float64, phi []float64) {
	eps2 := eps * eps
	for i := range x {
		for j := i + 1; j < len(x); j++ {
			dx := x[j][0] - x[i][0]
			dy := x[j][1] - x[i][1]
			dz := x[j][2] - x[i][2]
			dx2 := dx*dx + dy*dy + dz*dz

			phiij := 1.0 / math.Sqrt(dx2+eps2)
			phi[i] -= phiij
			phi[j] -= phiij
		}
	}
}

func dist2(x1, x2 [3]float64) float64 {
	dx2 := 0.0
	for k := 0; k < 3; k++ {
		dx := x2[k] - x1[k]
		dx2 += dx * dx
	}
	return dx2
}

func pointPotential(r2, eps2 float64) float64 {
	return -1.0 / math.Sqrt(r2+eps2)
}
