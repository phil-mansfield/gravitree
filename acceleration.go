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

func (acc Acceleration) Approximate(t1, t2 *Tree, i, j int) {
	switch t2.Order {
	case Monopole:
		nodei, nodej := &t2.Nodes[i], &t1.Nodes[j]
		xj := &nodej.Center
		massj := float64(nodej.End - nodej.Start)

		for i := nodei.Start; i < nodei.End; i++ {
			xi, idxi := &t2.Points[i], t2.Index[i]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = xi[k] - xj[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += t2.eps2

			for k := 0; k < 3; k++ {
				acc[idxi][k] -= massj * dx[k] / (dr2 * math.Sqrt(dr2))
			}
		}
	case Quadrupole:
		// TODO: Implement (Plummer) quadrupole appoximation of acceleration
		panic("NYI")
	default:
		panic(fmt.Sprintf("Unrecognized approximaiton order code, %d", t2.Order))
	}
}

func (acc Acceleration) OneSidedLeaf(t1, t2 *Tree, i, j int) {
	nodei, nodej := &t2.Nodes[i], &t1.Nodes[j]

	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t2.Points[i], t2.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := &t1.Points[j]

			dx := []float64{0, 0, 0}
			dr2 := 0.0
			for k := 0; k < 3; k++ {
				dx[k] = xi[k] - xj[k]
				dr2 += dx[k] * dx[k]
			}
			dr2 += t2.eps2

			for k := 0; k < 3; k++ {
				acc[idxi][k] -= dx[k] / (dr2 * math.Sqrt(dr2))
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
