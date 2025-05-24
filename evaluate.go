package gravitree

import (
	"fmt"
)

type Quantity interface {
	Len() int
	TwoSidedLeaf(t *Tree, i int)          // Evaluate quantity within i
	Approximate(t, t2 *Tree, i1, i2 int)  // Approximate quantity of j at i
	OneSidedLeaf(t, t2 *Tree, i1, i2 int) // Evaluate quantity of j at i
}

func (t1 *Tree) useApproximation(t2 *Tree, i1, i2 int) bool {
	node2, node1 := &t2.Nodes[i2], &t1.Nodes[i1]

	x_i1, x_i2 := &node1.Center, &node2.Center
	dx := x_i2[0] - x_i1[0]
	dy := x_i2[1] - x_i1[1]
	dz := x_i2[2] - x_i1[2]
	dx2 := dx*dx + dy*dy + dz*dz
	
	return dx2 > node2.RMax2+node1.ROpen2
}

func (t *Tree) Evaluate(eps float64, q Quantity) {
	if q.Len() != len(t.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(q) = %d",
			len(t.Nodes), q.Len()))
	}

	t.eps2 = eps * eps
	WorkerQueue(nWorkers, len(t.Nodes), func(worker, i int) {		
		if t.Nodes[i].Left == -1 {
			t.walkNodeEvaluate(0, i, q)
		}
	})
}


func (t *Tree) walkNodeEvaluate(i, j int, q Quantity) {
	target := &t.Nodes[i]
	
	if i == j {
		q.TwoSidedLeaf(t, i)
	} else if t.useApproximation(t, i, j) {
		q.Approximate(t, t, i, j) // passing in the same tree
	} else if target.Left == -1 {
		q.OneSidedLeaf(t, t, i, j) // passing in the same tree
	} else {
		t.walkNodeEvaluate(target.Left, j, q)
		t.walkNodeEvaluate(target.Right, j, q)
	}
}

// Same function as walkNodeEvaluate except it calculates quantities
// for a secondary tree.
func (t1 *Tree) walkNodeEvaluateAt(t2 *Tree, i1, i2 int, q Quantity) {
	target := &t1.Nodes[i1]

	if t1.useApproximation(t2, i1, i2) {
		q.Approximate(t1, t2, i1, i2) // passing in the secondary tree
	} else if target.Left == -1 {
		q.OneSidedLeaf(t1, t2, i1, i2)
	} else {
		t1.walkNodeEvaluateAt(t2, target.Left, i2, q)
		t1.walkNodeEvaluateAt(t2, target.Right, i2, q)
	}
}

func (t1 *Tree) EvaluateAt(t2 *Tree, eps float64, q Quantity) {
	if q.Len() != len(t2.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(q) = %d",
			len(t2.Nodes), q.Len()))
	}

	t1.eps2 = eps * eps

	WorkerQueue(nWorkers, len(t2.Nodes), func(worker, i2 int) {
		if t2.Nodes[i2].Left == -1 {
			t1.walkNodeEvaluateAt(t2, 0, i2, q)
		}
	})
}
