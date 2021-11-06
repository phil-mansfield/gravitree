package gravitree

import (
	"fmt"
	"math"
)

// OpeningCriteria represents a type of criteria used to decide whether a
// monopole approximation can be used for given tree node.
type OpeningCriteria int

const (
	// The PKDGRAV3 opening criteria (Potter, Stadel, & Teyssier 2017; S 3.1).
	// This computes R_i for each cell: R_i = R_max/theta, where R_max is
	// measured relative to the center of mass. Cells are opened if
	// R < R_target + 1.5*R_i.
	PKDGRAV3 OpeningCriteria = iota
	// This criteria is based on maximum-error analysis in Salmon & Warren
	// (1994), although rederiving it from that paper is tedious. A concise
	// description can be found in Behroozi, Wechsler, & Wu (2013); Appendix B.
	SalmonWarren
	// BarnesHut uses the classical Barnes-Hut opening criteria (Barnes & Hut
	// 1986). Cells are opened if R > MaxWidth/theta.
	BarnesHut
)

// Tree is a gravitational KD-tree which can be used to compute gravitaional
// forces and potentials.
type Tree struct {
	Root *Node // The root node of the tree.
	Nodes []Node // Array containing the Tree's Nodes.
	
	Points [][3]float64 // The (re-arranged) points in the Tree.
	Index []int // The original indices of points in the input array.
	
	LeafSize int // The maximum number of points that can be stored in a leaf.
	Criteria OpeningCriteria // Flag indicating the opening criteria.
	Theta float64 // The critical opening angle of the chosen criteria.

	eps2 float64
}

// Node is KD-node in a gravitational tree.
type Node struct {
	Center [3]float64 // Center of mass in the cell.
	RMax2, ROpen2 float64 // Squared radii used to determine cell opening.
	Left, Right int // The index of the left and right nodes 
	Start, End int // The indices of the points within the node in Tree.Points.
}

// TreeOptions allows the user to specify advanced options to tune performance
// and accuracy.
type TreeOptions struct {
	LeafSize int // Default: 16
	Criteria OpeningCriteria // Default: PKDGRAV
	Theta float64 // Default: 0.7

	// Buffers which can be reused between trees to reduce allocation. They
	// will be resized if the provided buffers are too small, but this can be
	// prevented by setting PointsBuffer and IndexBuffer to have length len(x),
	// and setting NodeBuffer to have length ceil(2*len(x)/LeafSize) or larger.
	// In practice, this is only useful if you're running a simulation. In this
	// case, just use the arrays from the previous Tree incarnation.
	PointsBuffer [][3]float64
	IndexBuffer []int
	NodeBuffer []Node
}

// NewTree creates a Tree from a colleciton of vectors, x. The tree can be
// customized through the optional TreeOptions parameter. Only the first
// TreeOptions argument will be used. If fields in the TreeOptions argument
// Are set to zero/nil, they will be replaced with the default values.
func NewTree(x [][3]float64, opt ...TreeOptions) *Tree {
	// Use default
	if len(opt) == 0 {
		opt = []TreeOptions{ {} }
	}
	if opt[0].LeafSize == 0 {
		opt[0].LeafSize = 16
	}
	if opt[0].Theta == 0.0 {
		opt[0].Theta = 0.7
	}
	
	t := &Tree{ Nodes: []Node{ }, LeafSize: opt[0].LeafSize,
		Theta: opt[0].Theta, Criteria: opt[0].Criteria}

	// Initialize points and indices.
	n := len(x)
	t.Points = append(opt[0].PointsBuffer[:0], make([][3]float64, n)...)
	t.Index = append(opt[0].IndexBuffer[:0], make([]int, n)...)
	
	nodeEstimate := int(math.Ceil(2*float64(len(x))/float64(t.LeafSize)))
	t.Nodes = append(opt[0].NodeBuffer[:0], make([]Node, nodeEstimate)...)
	t.Nodes = t.Nodes[:0]
	
	copy(t.Points, x)
	for i := range t.Index { t.Index[i] = i }

	if len(x) == 0 { return t }


	t.addNode(0, 0, len(x))
	t.Root = &t.Nodes[0]
	
	return t
}

// addNode adds a node to a tree which corresponds to points in the range
// [start: end] and a given depth.
func (t *Tree) addNode(depth, start, end int) {
	blankNode := Node{ [3]float64{}, 0, 0, -1, -1, start, end }
	span := pointSpan(t.Points[start: end])
	
	i := len(t.Nodes)
	t.Nodes = append(t.Nodes, blankNode)
	node := &t.Nodes[i]

	node.ROpen2 = t.rOpen2(i, span)
	
	if end - start <= t.LeafSize { return }
	
	dim := chooseNodeDimension(span)
	pivot := choosePivot(t.Points[start: end], span, dim)
	mid := partition(t.Points[start: end], t.Index[start: end], dim, pivot)
	
	node.Left = len(t.Nodes)
	t.addNode(depth+1, start, mid + start)
	t.Nodes[i].Right = len(t.Nodes)
	t.addNode(depth+1, mid + start, end)
}

// rOpen2 computes r_open^2 for node i with the given span.
func (t *Tree) rOpen2(i int, span [2][3]float64) float64 {
	node := &t.Nodes[i]

	pts := t.Points[node.Start: node.End]
	node.Center = centerOfMass(pts)
	node.RMax2 = rMax2(node.Center, pts)
	
	switch t.Criteria {
	case SalmonWarren:
		return t.rOpen2SalmonWarren(i, span)
	case BarnesHut:
		return t.rOpen2BarnesHut(i, span)
	case PKDGRAV3:
		return t.rOpen2PKDGRAV3(i, span)
	}
	panic(fmt.Sprintf("Unknown tree criteria %d.", t.Criteria))
}

// rOpen2PKDGRAV computes r_open^2 for the node, i, with span, span, using the
// Salmon-Warren monopole criteria.
func (t *Tree) rOpen2SalmonWarren(i int, span [2][3]float64) float64 {
	node := &t.Nodes[i]
	rMax := math.Sqrt(node.RMax2)
	
	sigmaX2 := 0.0
	for i := node.Start; i < node.End; i++ {
		for k := 0; k < 3; k++ {
			dx := node.Center[k] - t.Points[i][k]
			sigmaX2 += dx*dx
		}
	}
	sigmaX2 /= float64(node.End - node.Start)

	rOpen := rMax/2 + math.Sqrt(rMax*rMax/4 + sigmaX2/t.Theta)
	return rOpen*rOpen
}

// rOpen2PKDGRAV computes r_open^2 for the node, i, with span, span, using the
// classic Barnes-Hut criteria
func (t *Tree) rOpen2BarnesHut(i int, span [2][3]float64) float64 {
	width := span[1][0] - span[0][0]
	for k := 1; k < 3; k++ {
		dx := span[1][k] - span[0][k]
		if width < dx { dx = width }
	}

	return width*width / (t.Theta*t.Theta)
}

// rOpen2PKDGRAV3 computes r_open^2 for the node, i, with span, span, using the
// PKDGRAV3 criteria.
func (t *Tree) rOpen2PKDGRAV3(i int, span [2][3]float64) float64 {
	node := &t.Nodes[i]
	return 1.5*1.5 * node.RMax2 / (t.Theta*t.Theta)
}

// centerOfMass returns the center of mass for a collection of points, x.
func centerOfMass(x [][3]float64) [3]float64 {
	sum := &[3]float64{ }
	// Hot loop:
	for i := range x {
		xi := &x[i]

		sum[0] += xi[0]
		sum[1] += xi[1]
		sum[2] += xi[2]
	}

	n := float64(len(x))
	for k := 0; k < 3; k++ { sum[k] /= n }

	return *sum
}

// rMax2 returns the maximum squared distance between any point in x and the
// center of mass.
func rMax2(center [3]float64, x [][3]float64) float64 {
	rMax2 := 0.0
	c := &center
	// Hot loop:
	for i := range x {
		xi := &x[i]
		dx := xi[0] - c[0]
		dy := xi[1] - c[1]
		dz := xi[2] - c[2]
		dx2 := dx*dx + dy*dy + dz*dz
		
		if dx2 > rMax2 { rMax2 = dx2 }
	}

	return rMax2
}

// pointSpan returns the span of a collection of points.
func pointSpan(x [][3]float64) [2][3]float64 {
	if len(x) == 0 { return [2][3]float64{ } }
	
	min, max := x[0], x[0]
	
	// Host loop:
	for i := 1; i < len(x); i++ {
		xi := x[i]
		if xi[0] < min[0] {
			min[0] = xi[0]
		} else if xi[0] > max[0] {
			max[0] = xi[0]
		}
		if xi[1] < min[1] {
			min[1] = xi[1]
		} else if xi[1] > max[1] {
			max[1] = xi[1]
		}
		if xi[2] < min[2] {
			min[2] = xi[2]
		} else if xi[2] > max[2] {
			max[2] = xi[2]
		}
	}

	return [2][3]float64{ min, max }
}

// chooseNodeDimension returns the dimension that a node with the given span
// should be split along.
func chooseNodeDimension(span [2][3]float64) int {
	width := [3]float64{ }
	for k := 0; k < 3; k++ {
		width[k] = span[1][k] - span[0][k]
	}

	if width[1] >= width[0] && width[1] >= width[2] {
		return 1
	} else if width[2] >= width[0] && width[2] >= width[1] {
		return 2
	}
	return 0
}

// choosePivot returns a pivot value in dimension dim for a node containing
// the points x with the span span.
func choosePivot(x [][3]float64, span [2][3]float64, dim int) float64 {
	// Midpoint of the span.
	width := span[1][dim] - span[0][dim]
	mid := span[0][dim] + width/2
	return mid
}

// partition partitions the array x into a "left" sub array where each element
// <= pivot and a "right" sub array > pivot. This is evaluated in the dim
// dimension.  The length of the left sub array is returned. Both x and x are
// reordered accordingly.
func partition(x [][3]float64, idx []int, dim int, pivot float64) int {
	// Small arrays need to be handled manually.
	switch len(x) {
	case 0:
		return 0
	case 1:
		if x[0][dim] <= pivot { return 1 }
		return 0
	}

	// Find the numebr of objects <= the pivot
	left := 0
	for i := range x {
		if x[i][dim] <= pivot { left++ }
	}

	// Ending early on this condition avoids an extra bounds check in the
	// main loop.
	if left == 0 { return 0 }	
	l, r := 0, len(x) - 1
	for {
		// Loop until you find a pair that needs to be switched.
		for x[l][dim] <= pivot {
			l++
			if l == left { return left }
		}
		for x[r][dim] > pivot {
			r--
		}

		// Switch the pair.
		x[l], x[r] = x[r], x[l]
		idx[l], idx[r] = idx[r], idx[l]
		l++
		r--
		if l == left { return left }
	}
}
