package gravitree

import (
	"fmt"
	"math"
)

type Acceleration [][3]float64

var _ Quantity = Acceleration{}

func (acc Acceleration) Len() int { return len(acc) }

func (acc Acceleration) TwoSidedLeaf(t *Tree, i int) {
	node := &t.Nodes[i]
	for i := node.Start; i < node.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := i + 1; j < node.End; j++ {
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

func (acc Acceleration) Approximate(t1, t2 *Tree, i2, i1 int) {
	switch t2.Order {
	case Monopole:
		node_i2, node_i1 := &t2.Nodes[i2], &t1.Nodes[i1]
		x_i1 := &node_i1.Center
		mass_i1 := float64(node_i1.End - node_i1.Start)

		for i2 := node_i2.Start; i2 < node_i2.End; i2++ {
			x_i2, idx_i2 := &t2.Points[i2], t2.Index[i2]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = x_i2[k] - x_i1[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += t2.eps2

			for k := 0; k < 3; k++ {
				acc[idx_i2][k] -= mass_i1 * dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	case Quadrupole:
		// TODO: Implement (Plummer) quadrupole appoximation of acceleration
		panic("NYI")
	default:
		panic(fmt.Sprintf("Unrecognized approximaiton order code, %d", t2.Order))
	}
}

func (acc Acceleration) OneSidedLeaf(t1, t2 *Tree, i2, i1 int) {
	node_i2, node_i1 := &t2.Nodes[i2], &t1.Nodes[i1]

	for i2 := node_i2.Start; i2 < node_i2.End; i2++ {
		x_i2, idx_i2 := &t2.Points[i2], t2.Index[i2]
		for i1 := node_i1.Start; i1 < node_i1.End; i1++ {
			x_i1 := &t1.Points[i1]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = x_i2[k] - x_i1[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += t2.eps2

			for k := 0; k < 3; k++ {
				acc[idx_i2][k] -= dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	}
}

func BruteForceAcceleration(eps float64, x [][3]float64, acc [][3]float64) {
	eps2 := eps * eps
	for i := range x {
		xi := x[i]
		for j := i + 1; j < len(x); j++ {
			xj := x[j]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = xi[k] - xj[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += eps2

			for k := 0; k < 3; k++ {
				acc[i][k] -= dx[k] / (dr2 * math.Sqrt(dr2))
				acc[j][k] += dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	}
}

func BruteForceAccelerationAt(eps float64, x1, x2 [][3]float64, acc [][3]float64) {
	eps2 := eps * eps
	for i := range x1 {
		xi := x1[i]
		for j := range x2 {
			xj := x2[j]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = xi[k] - xj[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += eps2

			for k := 0; k < 3; k++ {
				acc[j][k] += dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	}
}
