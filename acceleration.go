package gravitree

import (
	"math"
)

// Implementations of acceleration vectors

// Acceleration computed for each point in the tree. Vectors are
// written to the 3-vector acc and are in units where G * mp = 1.


type Acceleration [][3]float64
var _ Quantity = Acceleration{ }

func (acc Acceleration) Len() int { return len(acc) }

func (acc Acceleration) TwoSidedLeaf(t *Tree, i int) {
	node := &t.Nodes[i]
	for i := node.Start; i < node.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := i+1; j < node.End; j++ {
			xj, idxj := &t.Points[j], t.Index[j]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = xi[k] - xj[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += t.eps2

			for k := 0; k < 3; k++ {
				acc[idxi][k] -= dx[k] / (dr2 * math.Sqrt(dr2))
				acc[idxj][k] += dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	}
}


func (acc Acceleration) Approximate(t *Tree, i, j int) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	xj := &nodej.Center
	massj := float64(nodej.End - nodej.Start)

	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]

		dx := []float64{0, 0, 0}
		dr2 := 0.0
		for k := 0; k < 3; k++ {
			dx[k] = xi[k] - xj[k]
			dr2 += dx[k] * dx[k]
		}
		dr2 += t.eps2

		for k := 0; k < 3; k++ {
			acc[idxi][k] -= massj*dx[k] / (dr2 * math.Sqrt(dr2))
		}
	}
}

func (acc Acceleration) OneSidedLeaf(t *Tree, i, j int) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]

	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := &t.Points[j]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = xi[k] - xj[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += t.eps2

			for k := 0; k < 3; k++ {
				acc[idxi][k] -= dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	}
}