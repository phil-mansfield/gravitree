package gravitree

import "fmt"

type Quantity interface {
	Len() int
	TwoSidedLeaf(t *Tree, i int)        // Evaluate quantity within i
	Approximate(t, t2 *Tree, i, j int)  // Approximate quantity of j at i
	OneSidedLeaf(t, t2 *Tree, i, j int) // Evaluate quantity of j at i
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
		q.Approximate(t, t, i, j) // passing in the same tree
	} else if target.Left == -1 {
		q.OneSidedLeaf(t, t, i, j) // passing in the same tree
	} else {
		t.walkNodeEvaluate(i, target.Left, q)
		t.walkNodeEvaluate(i, target.Right, q)
	}
}

// Same function as walkNodeEvaluate except it calculates quantities
// for a secondary tree.
func (t1 *Tree) walkNodeEvaluateAt(t2 *Tree, i, j int, q Quantity) {
	target := &t2.Nodes[j]

	if t1.useApproximation(i, j) {
		q.Approximate(t1, t2, i, j) // passing in the secondary tree
	} else if target.Left == -1 {
		q.OneSidedLeaf(t1, t2, i, j)
	} else {
		t1.walkNodeEvaluateAt(t2, i, target.Left, q)
		t1.walkNodeEvaluateAt(t2, i, target.Right, q)
	}
}

func (t1 *Tree) EvaluateAt(t2 *Tree, eps float64, q Quantity) {
	// TODO: Write an alternate version of Evaluate which takes in a set of test
	// points and updates the Quantity accordingly. For now, let's have the other points
	// passed to us as a Tree. We'll want to do new benchmarking of how much refinement
	// this new tree should have. My guess is that for small point sets, we'll want to
	// set the minimum leaf size to 1.
	//
	// This new funciton will only need to call OneSidedLeaf and Approximate. Those
	// two methods of the Quantity interface will need to be changed so they take in
	// a separte tree for index i and and for index j (and then the vanilla
	// walkNodeEvaluate will need to be updated so that the same tree is passed twice.)
	//
	// (It's totally possible that we want this to work without the secondary tree.)
	// panic("NYI")

	if q.Len() != len(t2.Points) {
		panic(fmt.Sprintf("Tree has %d points, but len(q) = %d",
			len(t2.Nodes), q.Len()))
	}

	t1.eps2 = eps * eps

	// loop over the t2 nodes
	for i := range t2.Nodes {
		if t2.Nodes[i].Left == -1 {
			t1.walkNodeEvaluateAt(t2, i, 0, q)
		}
	}
}
