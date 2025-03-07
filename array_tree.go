package gravitree

// An ArrayTree is a wrapper around an array of particles which (a) allows the
// Quantity interface to interact with them and (b) allows particles to be
// turned off. Each particle will be put in its own isolated node and the
// array will be used directly by gravitree without copying. This tree can 
// only support being passed to EvaluateAt as a second argument.
//
// A typical usage of Array tree is when updating a set of test paritcles
// orbiting a set of massive paritlces. In this case, using ArrayTree
// would look something like this:
//
//  t1 := ... // Tree contining particles that the test points are orbiting
//  x, v, acc := ... // initial positions, velocities, and accelerations
//                  // of the test points
//  t2 := NewArrayTree(x)
//  t1.EvaluateAt(t2, acc)
//  makeTimestep(x, v, acc) // User function which updates x and v
//  t2.Update() // Updates tree based on values in x
//  t1.EvaluateAt(t2, acc) // t2 can be used again
//  // etc.
//  
type ArrayTree struct {
	Tree
}

// NewArrayTree creates an ArrayTree wrapper around an array of position
// vectors. The evaluate flag is set to true by default.
func NewArrayTree(x [][3]float64) *ArrayTree {
	idx := make([]int, len(x))
	for i := range idx { idx[i] = i }
	nodes := make([]Node, len(x))
	for i := range nodes {
		nodes[i].Start, nodes[i].End = i, i+1
		nodes[i].Left, nodes[i].Center = -1, x[i]
	}
	return  &ArrayTree{ Tree{ Points: x, Index: idx, Nodes: nodes } }
}

// Update update's an ArrayTree's internal information to be consistent with
// the current values of its input array.
func (t *ArrayTree) Update() {
	for i := range t.Nodes {
		t.Nodes[i].Center = t.Points[i]
	}
}

// SetEvaluateFlag sets the evaluation flag for each particle. Particles with
// ok = true will have quantities evaluated at their location and particles
// with ok = false will not.
func (t *ArrayTree) SetEvaluateFlag(ok []bool) {
	for i := range t.Nodes {
		if ok[i] {
			t.Nodes[i].Left = -1
		} else {
			t.Nodes[i].Left = 0
		}
	}
}