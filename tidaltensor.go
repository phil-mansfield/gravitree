package gravitree

import (
	"fmt"
	"math"
)

// Tidal Tensor
func (t *Tree) TidalTensor(eps float64, T [][3][3]float64) {
	if len(T) != len(t.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(T) = %d",
			len(t.Nodes), len(T)))
	}

	t.eps2 = eps * eps

	for i := range t.Nodes {
		// Only compute the tensor at the leaf nodes.
		if t.Nodes[i].Left == -1 {
			t.walkNodeTidalTensor(i, 0, T)
		}
	}
}

func (t *Tree) walkNodeTidalTensor(i, j int, T [][3][3]float64) {
	target := &t.Nodes[j]

	if i == j {
		t.pairwiseTidalTensor(i, T)
	} else if t.useMonopole(i, j) {
		t.monopoleTidalTensor(i, j, T)
	} else if target.Left == -1 {
		t.oneSidedTidalTensor(i, j, T)
	} else {
		t.walkNodeTidalTensor(i, target.Left, T)
		t.walkNodeTidalTensor(i, target.Right, T)
	}
}

func pointTidalTensor(x [3]float64, eps2 float64) [3][3]float64 {

	a := (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + eps2)
	a = a * a * math.Sqrt(a)

	x2 := x[0] * x[0]
	y2 := x[1] * x[1]
	z2 := x[2] * x[2]

	tensor := [3][3]float64{
		{y2 + z2 - 2*x2, -3 * x[0] * x[1], -3 * x[0] * x[2]},
		{-3 * x[0] * x[1], x2 + z2 - 2*y2, -3 * x[1] * x[2]},
		{-3 * x[0] * x[2], -3 * x[1] * x[2], x2 + y2 - 2*z2},
	}

	// rescale
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			tensor[i][j] /= a
		}
	}

	return tensor
}

func (t *Tree) oneSidedTidalTensor(i, j int, T [][3][3]float64) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]

	// This loop is very hot.
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := &t.Points[j]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]

			dt := pointTidalTensor([3]float64{dx, dy, dz}, t.eps2)

			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					T[idxi][k][l] += dt[k][l]
				}
			}
		}
	}
}

func (t *Tree) monopoleTidalTensor(i, j int, T [][3][3]float64) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	xj := &nodej.Center
	massj := float64(nodej.End - nodej.Start)

	// This loop is very hot.
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]

		dx := xj[0] - xi[0]
		dy := xj[1] - xi[1]
		dz := xj[2] - xi[2]

		dt := pointTidalTensor([3]float64{dx, dy, dz}, t.eps2)

		for k := 0; k < 3; k++ {
			for l := 0; l < 3; l++ {
				T[idxi][k][l] += dt[k][l] * massj
			}
		}
	}
}

func (t *Tree) pairwiseTidalTensor(i int, T [][3][3]float64) {
	node := &t.Nodes[i]
	for i := node.Start; i < node.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := i + 1; j < node.End; j++ {
			xj, idxj := &t.Points[j], t.Index[j]

			dx := xj[0] - xi[0]
			dy := xj[1] - xi[1]
			dz := xj[2] - xi[2]

			dt := pointTidalTensor([3]float64{dx, dy, dz}, t.eps2)

			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					T[idxi][k][l] += dt[k][l]
					T[idxj][k][l] += dt[k][l]
				}
			}
		}
	}
}
