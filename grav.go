package gravitree

import (
	"fmt"
	"math"
)

// CenterOfMass calculates the center of mass for each Node in tree.
func (tree *KDTree) CenterOfMass(out [][3]float64) {
	N := len(tree.Nodes)
	if len(out) != N {
		panic(fmt.Sprintf("len(tree.Nodes) = %d, but len(out) = %d",
			len(tree.Nodes), len(out)))
	}
	
	// There's no need to iterate through the tree recursively due to the
	// memory layout of the nodes. Helper functions are inlined.
	for i := N - 1; i >= 0; i-- {
		node := &tree.Nodes[i]
		
		if node.Left == -1 {
			// Leaf. Manually find center of mass.
			
			if node.Start == node.End { continue } // Empty leaf.
			if node.Start + 1 == node.End {
				x := tree.Points[node.Start]
				out[i] = x
				continue
			}
			
			for j := node.Start; j < node.End; j++ {
				for d := 0; d < 3; d++ {
					out[i][d] += tree.Points[j][d]
				}
			}
			n := float64(node.End - node.Start)
			for d := 0; d < 3; d++ {
				out[i][d] /= n
			}
		} else {
			// Node. Compute center of mass from the pivot + the left and right
			// monopoles.
			
			out[i] = tree.Points[node.Pivot]
			nLeft := float64(node.Pivot - node.Start)
			nRight := float64(node.End - node.Pivot - 1)
			
			for d := 0; d < 3; d++ {
				out[i][d] += nLeft*out[node.Left][d]
				out[i][d] += nRight*out[node.Right][d]
				out[i][d] /= 1 + nLeft + nRight
			}
		}
	}
}

// Span computes the span of each Node in tree. span[i][0] is the lower corner
// of the node and span[i][1] is the upper corner.
func (tree *KDTree) Span(out [][2][3]float64) {
	N := len(tree.Nodes)

	if len(out) != N {
		panic(fmt.Sprintf("len(tree.Nodes) = %d, but len(out) = %d",
			len(tree.Nodes), len(out)))
	}

	for i := N - 1; i >= 0; i-- {
		node := &tree.Nodes[i]
		
		if node.Left == -1 {
			// Leaf. Manually compute span.
			
			if node.Start == node.End { continue } // Empty node.
			if node.Start + 1 == node.End {
				x := tree.Points[node.Start]
				out[i][0], out[i][1] = x, x
				continue
			}
			
			// Start with the first point in the range.
			
			low, high := tree.Points[node.Start], tree.Points[node.Start]
			for j := node.Start+1; j < node.End; j++ {
				x := tree.Points[j]
				for d := 0; d < 3; d++ {
					if x[d] < low[d] {
						low[d] = x[d]
					} else if x[d] > high[d] {
						high[d] = x[d] }
				}
			}
			out[i][0], out[i][1] = low, high
		} else {
			//  Node. Combine pivot with the left and right spans.
			low, high := tree.Points[node.Pivot], tree.Points[node.Pivot]
			left, right := &tree.Nodes[node.Left], &tree.Nodes[node.Right]
			// In principle, unrolling these d loops, skipping d = node.Dim,
			// and updating those ranges manually would speed this up by a
			// factor of 33%.
			if left.Start != left.End { // Non-empty node.
				leftLow, leftHigh  := out[node.Left][0], out[node.Left][1]
				for d := 0; d < 3; d++ {
					if leftLow[d] < low[d] {
						low[d] = leftLow[d]
					} else if leftHigh[d] > high[d] {
						high[d] = leftHigh[d]
					}
				}
			}
			if right.Start != right.End { // Non-empty node.
				rightLow, rightHigh  := out[node.Right][0], out[node.Right][1]
				for d := 0; d < 3; d++ {
					if rightLow[d] < low[d] {
						low[d] = rightLow[d]
					} else if rightHigh[d] > high[d] {
						high[d] = rightHigh[d]
					}
				}
			}
			out[i][0], out[i][1] = low, high
		}
	}
}

// SquaredMonopoleDistance computes the squared distances above which each
// node in the tree can be safely approximated by a monopole. Results are
// written to out. Errors in the potential are bounded to be smaller than
// relErr according to the monopole criteria in Salmon & Warren (1994).
func (tree *KDTree) SquaredMonopoleDistance(
	relErr float64, centers [][3]float64, out []float64,
) {	
	for i := range tree.Nodes {
		node := &tree.Nodes[i]
		if node.End - node.Start <= 1 {
			continue
		}
		
		rMax2, sigmaX2 := tree.rStats(centers, i)

		term1 := rMax2/4
		term3 := rMax2/4 + sigmaX2/relErr
		term2 := 2*math.Sqrt(term1*term3)
		
		out[i] = term1 + term2 + term3
	}
}

// rStats returns R_max^2 and sigma_x^2 for node i in a KDTree. These quantities
// are defined in Salmon & Warren (1994): R_max = |x - x_com|^2, and
// sigma_x = < |x - x_com|^2 >
func (tree *KDTree) rStats(centers [][3]float64, i int) (
	rMax2, sigmaX2 float64,
) {
	node := &tree.Nodes[i]
	x0 := centers[i]
	sigmaX2, rMax2 = 0.0, 0.0
	for j := node.Start; j < node.End; j++ {
		r2 := 0.0
		for k := 0; k < 3; k++ {
			dx := tree.Points[j][k] - x0[k]
			r2 += dx*dx
		}

		sigmaX2 += r2
		if r2 > rMax2 { rMax2 = r2 }
	}

	n := float64(node.End - node.Start)
		
	return rMax2, sigmaX2 / n
}

// Gravitree is a full gravitational KD-tree. Most users will only need to use
// this type and its methods.
type Gravitree struct {
	Tree *KDTree
	Center [][3]float64
	SMD []float64 // SquaredMonopoleDistance
	Mp float64
}


// NewGravitree creates a new Gravitree.
func NewGravitree(kt *KDTree, relAcc float64,) *Gravitree {
	gt := &Gravitree{ }
	gt.Tree = kt
	
	gt.Center = make([][3]float64, len(gt.Tree.Nodes))
	gt.SMD = make([]float64, len(gt.Tree.Nodes))

	gt.Tree.CenterOfMass(gt.Center)
	gt.Tree.SquaredMonopoleDistance(relAcc, gt.Center, gt.SMD)
	
	return gt
}


var leaf, open, closed = 0, 0, 0
// Potential calculates the potential at x using Plummer potentials for each
// prticle. Potential is returned in units where G Mp^2 is one.
func (gt *Gravitree) Potential(eps float64, x [3]float64) float64 {
	leaf, open, closed = 0, 0, 0
	phi := gt.potential(eps*eps, x, 0)
	return phi
}

// Potential calculates the potential at x using Plummer potentials for each
// particle using only the particles within node i of gt.Tree.
func (gt *Gravitree) potential(eps2 float64, x [3]float64, i int) float64 {
	node := &gt.Tree.Nodes[i]
	pts := gt.Tree.Points
	sum := 0.0

	// If the node is far away, use the monopole approximation.
	r2 := dist2(x, gt.Center[i])
	if r2 > gt.SMD[i] {
		closed++
		n := float64(node.End - node.Start)
		return pointPotential(r2, eps2) * n
	}
	
	// If this is leaf node, just loop over it.
	if node.Left == -1 {
		leaf++
		for j := node.Start; j < node.End; j++ {
			sum += pointPotential(dist2(x, pts[j]), eps2)
		}
		return sum
	}

	// If the node is close, open the node.
	piv := gt.Tree.Points[node.Pivot]
	open++
	return gt.potential(eps2, x, node.Left) +
		gt.potential(eps2, x, node.Right) +
		pointPotential(dist2(x, piv), eps2)
}

// dist2 computes the squared distance between x1 and x2.
func dist2(x1, x2 [3]float64) float64 {
	r2 := 0.0
	for k := 0; k < 3; k++ {
		dx := x1[k] - x2[k]
		r2 += dx*dx
	}
	return r2
}

// pointPotential computes the potential at a squared distance, r2 and force
// softening, eps.
func pointPotential(r2 float64, eps2 float64) float64 {
	return -math.Pow(r2 + eps2, -0.5)
}
