package gravitree

import (
	"fmt"
	"math"
)

// Implementations of acceleration vectors

// Acceleration computed for each point in the tree. Vectors are
// written to the 3-vector acc and are in units where G * mp = 1.

func (t *Tree) Acceleration(eps float64, acc [][3]float64) {
	if len(acc) != len(t.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(acc) = %d",
			len(t.Nodes), len(acc)))
	}

	t.eps2 = eps * eps

	for i := range t.Nodes {
		// Only compute the potential at the leaf nodes.
		if t.Nodes[i].Left == -1 {
			t.walkNodeAcceleration(i, 0, acc)
		}
	}
}

func (t *Tree) walkNodeAcceleration(i, j int, acc [][3]float64) {
	target := &t.Nodes[j]

	if i == j {
		t.pairwiseAcceleration(i, acc)
	} else if t.useMonopole(i, j) {
		t.monopoleAcceleration(i, j, acc)
	} else if target.Left == -1 {
		t.oneSidedAcceleration(i, j, acc)
	} else {
		t.walkNodeAcceleration(i, target.Left, acc)
		t.walkNodeAcceleration(i, target.Right, acc)
	}
}

func pointAcceleration(x [3]float64, eps2 float64) [3]float64 {
	r2 := x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + eps2
	r := math.Sqrt(r2)
	return [3]float64{x[0] / (r * r2), x[1] / (r * r2), x[2] / (r * r2)}
}

func (t *Tree) oneSidedAcceleration(i, j int, acc [][3]float64) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]

	// This loop is very hot.
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := &t.Points[j]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]

			da := [3]float64{dx, dy, dz}
			da = pointAcceleration(da, t.eps2)

			for k := 0; k < 3; k++ {
				acc[idxi][k] += da[k]
			}
		}
	}
}

func (t *Tree) monopoleAcceleration(i, j int, acc [][3]float64) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	xj := &nodej.Center
	massj := float64(nodej.End - nodej.Start)

	// This loop is very hot.
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]

		dx := xj[0] - xi[0]
		dy := xj[1] - xi[1]
		dz := xj[2] - xi[2]

		da := [3]float64{dx, dy, dz}
		da = pointAcceleration(da, t.eps2)

		for k := 0; k < 3; k++ {
			acc[idxi][k] += da[k] * massj
		}
	}
}

// pairwisePotential computes the potential at every point in node i using the
// other points in node i.
func (t *Tree) pairwiseAcceleration(i int, acc [][3]float64) {
	node := &t.Nodes[i]
	for i := node.Start; i < node.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := i + 1; j < node.End; j++ {
			xj, idxj := &t.Points[j], t.Index[j]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]

			da := [3]float64{dx, dy, dz}
			da = pointAcceleration(da, t.eps2)

			for k := 0; k < 3; k++ {
				acc[idxi][k] += da[k]
				acc[idxj][k] -= da[k] // reflexive
			}
		}
	}
}
