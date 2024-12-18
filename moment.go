package gravitree

import "fmt"

type GravityMoment interface {
	// You may want to change this so you pass the distance between xi and xj
	// in as a parameter to avoid needing to recompute it.
	Len() int
	AddMonopoleOneSided(i int, xi, xj *[3]float64, weight, eps2 float64)
	AddMonopoleTwoSided(i, j int, xi, xj *[3]float64, weight, eps2 float64)
}

type Potential []float64
type Acceleration [][3]float64
type TidalTensor [][3][3]float64

// useMonopole returns true if node i is close enough to node j that a monopole
// approximation can be used.
func (t *Tree) useMonopole(i, j int) bool {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]

	xi, xj := &nodei.Center, &nodej.Center
	dx := xj[0] - xi[0]
	dy := xj[1] - xi[1]
	dz := xj[2] - xi[2]
	dx2 := dx*dx + dy*dy + dz*dz

	return dx2 > nodei.RMax2+nodej.ROpen2 || dx2 < t.eps2
}

func (t *Tree) Quantity(eps float64, quant GravityMoment) {
	if quant.Len() != len(t.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(quant) = %d",
			len(t.Nodes), quant.Len()))
	}

	t.eps2 = eps * eps

	for i := range t.Nodes {
		// Only compute the potential at the leaf nodes.
		if t.Nodes[i].Left == -1 {
			t.walkNode(i, 0, quant)
		}
	}
}

func (t *Tree) walkNode(i, j int, quant GravityMoment) {
	target := &t.Nodes[j]

	if i == j {
		t.Pairwise(i, quant)
	} else if t.useMonopole(i, j) {
		t.Monopole(i, j, quant)
	} else if target.Left == -1 {
		t.OneSided(i, j, quant)
	} else {
		t.walkNode(i, target.Left, quant)
		t.walkNode(i, target.Right, quant)
	}
}

func (t *Tree) Pairwise(i int, quant GravityMoment) {
	node := &t.Nodes[i]

	for i := node.Start; i < node.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := i + 1; j < node.End; j++ {
			xj, idxj := &t.Points[j], t.Index[j]
			quant.AddMonopoleTwoSided(idxi, idxj, xi, xj, 1, t.eps2)
		}
	}
}

func (t *Tree) Monopole(i, j int, quant GravityMoment) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	xj := &nodej.Center
	massj := float64(nodej.End - nodej.Start)

	// This loop is very hot.
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		quant.AddMonopoleOneSided(idxi, xi, xj, massj, t.eps2)
	}
}

func (t *Tree) OneSided(i, j int, quant GravityMoment) {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]
	for i := nodei.Start; i < nodei.End; i++ {
		xi, idxi := &t.Points[i], t.Index[i]
		for j := nodej.Start; j < nodej.End; j++ {
			xj := &t.Points[j]
			quant.AddMonopoleOneSided(idxi, xi, xj, 1, t.eps2)
		}
	}
}
