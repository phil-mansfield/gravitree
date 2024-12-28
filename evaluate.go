package gravitree

import "fmt"

type Quantity interface {
	Len() int
	TwoSidedLeaf(t *Tree, i int) // Evaluate quantity within i
	Approximate(t *Tree, i, j int) // Approximate quantity of j at i
	OneSidedLeaf(t *Tree, i, j int) // Evaluate quantity of j at i
}

func (t *Tree) useApproximation(i, j int) bool {
	nodei, nodej := &t.Nodes[i], &t.Nodes[j]

	xi, xj := &nodei.Center, &nodej.Center
	dx := xj[0] - xi[0]
	dy := xj[1] - xi[1]
	dz := xj[2] - xi[2]
	dx2 := dx*dx + dy*dy + dz*dz

	return dx2 > nodei.RMax2+nodej.ROpen2 || dx2 < t.eps2
}

func (t *Tree) Evaluate(eps float64, q Quantity) {
	if q.Len() != len(t.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(q) = %d",
			len(t.Nodes), q.Len()))
	}

	t.eps2 = eps * eps

	for i := range t.Nodes {
		// Only compute the potential at the leaf nodes.
		if t.Nodes[i].Left == -1 {
			t.walkNodeEvaluate(i, 0, q)
		}
	}
}

func (t *Tree) walkNodeEvaluate(i, j int, q Quantity) {
	target := &t.Nodes[j]

	if i == j {
		q.TwoSidedLeaf(t, i)
	} else if t.useApproximation(i, j) {
		q.Approximate(t, i, j)
	} else if target.Left == -1 {
		q.OneSidedLeaf(t, i, j)
	} else {
		t.walkNodeEvaluate(i, target.Left, q)
		t.walkNodeEvaluate(i, target.Right, q)
	}
}
