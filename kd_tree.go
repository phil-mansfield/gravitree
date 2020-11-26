package gravitree

import (
	"fmt"
)

// A k-d tree (https://en.wikipedia.org/wiki/K-d_tree). Since this tree is only
// every used for gravitational calculations and 3d bounds checking, the number
// of dimensions is fixed at 3.
type KDTree struct {
	Root *KDNode // The root node.
	Points [][3]float64 // The points in the k-d tree, reordered by node.
	Nodes []KDNode // Flat array of the KDNodes.
	LeafSize int
}

// NewKDTree creates a KDTree from a colleciton of vectors, x. The footprint of
// the tree is given by span (first element is the lowermost corner, and the
// uppermost corner is the second element). 
func NewKDTree(x [][3]float64, span [2][3]float64, leafSize int) *KDTree {
	if leafSize <= 0 {
		panic(fmt.Sprintf("leafSize is %d", leafSize))
	}
	
	buf := make([][3]float64, len(x))
	copy(buf, x)
	
	t:= &KDTree{ Points: buf, Nodes: []KDNode{ }, LeafSize: leafSize }
	if len(x) > 0 {
		depth = 0
		t.addNode(0, len(x), span)
		t.Root = &t.Nodes[0]
	}

	return t
}

var depth = 0

// addNode adds a single node to a KDTree instance. parent is the index of the
// parent node, start and end give the range of elements in t.Points that are
// inside the node, and span is the geometric range of the node.
func (t *KDTree) addNode(start, end int, span [2][3]float64) {
	
	depth++
	blankNode := KDNode{ -1, -1, -1, -1, start, end}
	i := len(t.Nodes)
	t.Nodes = append(t.Nodes, blankNode)
	node := &t.Nodes[i]
	
	if end - start <= t.LeafSize { return }
	
	node.Dim = chooseNodeDimension(span)
	node.Pivot = partition(t.Points[start: end], node.Dim) + start
	leftSpan, rightSpan := t.splitSpan(node, span)
	
	node.Left = len(t.Nodes)
	t.addNode(start, node.Pivot, leftSpan)
	depth--
	
	// Can't use "node" here because addNode performs an append.
	t.Nodes[i].Right = len(t.Nodes)
	t.addNode(node.Pivot+1, end, rightSpan)
	depth--
}

// splitSpan splits a span in two at a node's pivot point and returns the
// lower and upper portions.
func (t *KDTree) splitSpan(n *KDNode, span [2][3]float64) (
	left, right [2][3]float64,
) {
	pivot := t.Points[n.Pivot][n.Dim]
	left, right = span, span
	left[1][n.Dim], right[0][n.Dim] = pivot, pivot
	return left, right
}

// chooseNodeDimension returns the dimension that a node with a givne span
// should be split in.
func chooseNodeDimension(span [2][3]float64) int {
	width := [3]float64{ }
	for i := range width {
		width[i] = span[1][i] - span[0][i]
	}
	
	if width[1] > width[0] && width[1] > width[2] {
		return 1
	} else if width[2] > width[0] && width[2] > width[1] {
		return 2
	}
	return 0
}


// A single node in the k-d tree
type KDNode struct {
	Left, Right int // Indices of the left and right nodes in KDTree.Nodes.
	
	Dim int // The dimension that the split is being performed over.
	Pivot int // The index of this node's pivot point in KDTree.Points.
	Start, End int // The starting and ending index in KDTree.Points with this
	               // node's points.
}

// sort3 sorts three values from largest to smallest.
func sort3(x, y, z [3]float64, d int) (max, mid, min [3]float64) {
	if x[d] > y[d] {
		if x[d] > z[d] {
			if y[d] > z[d] {
				return x, y, z
			} else {
				return x, z, y
			}
		} else {
			return z, x, y
		}
	} else {
		if y[d] > z[d] {
			if x[d] > z[d] {
				return y, x, z
			} else {
				return y, z, x
			}
		} else {
			return z, y, x
		}
	}
}

// partition rearranges the elements of a slice, x, into two contiguous
// groups, such that the d'th element of every vector the first group is <=
// to the pivot and every vector in the second group is >= the pivot. This
// function returns the index of the pivot element. Everything before it is in
// the first group and everything after it is in the second group.
func partition(x [][3]float64, d int) int {
	switch len(x) {
	case 0, 1:
		return 0
	case 2:
		if x[0][d] > x[1][d] {
			x[1], x[0] = x[0], x[1]
		}
		return 0
	case 3:
		max, mid, min := sort3(x[0], x[1], x[2], d)
		x[0], x[1], x[2] = min, mid, max
		return 1
	}
	
	n, n2 := len(x), len(x)/2
	// Take three values. The median will be the pivot, the other two will
	// be sentinel values so that we can avoid additional bounds checks. This is
	// based on the Numerical Recipes median algorithm.
	max, mid, min := sort3(x[0], x[n2], x[n-1], d)
	x[0], x[n2], x[n-1] = min, mid, max
	x[1], x[n2] = x[n2], x[1]

	lo, hi := 1, n-1
	for {
		lo++
		for x[lo][d] < mid[d] {
			lo++
		}
		hi--
		for x[hi][d] > mid[d] {
			hi--
		}
		if hi < lo {
			break
		}
		x[lo], x[hi] = x[hi], x[lo]
	}

	// Swap the pivot into the middle
	x[1], x[hi] = x[hi], x[1]

	return hi
}
