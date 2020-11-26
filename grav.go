package gravitree

import (
	"fmt"
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

			// Start with the first point in the range.
			low, high := tree.Points[node.Start], tree.Points[node.Start]
			for j := node.Start+1; j < node.End; j++ {
				x := tree.Points[j]
				
				// Loops are split up for vectorization.
				for d := 0; d < 3; d++ {
					if x[d] < low[d] { low[d] = x[d] }
				}
				for d := 0; d < 3; d++ {
					if x[d] > high[d] { high[d] = x[d] }
				}
			}
			out[i] = [2][3]float64{ low, high }
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
					if leftLow[d] < low[d] { low[d] = leftLow[d] }
				}
				for d := 0; d < 3; d++ {
					if leftHigh[d] > high[d] { high[d] = leftHigh[d] }
				}
			}
			if right.Start != right.End { // Non-empty node.
				rightLow, rightHigh  := out[node.Right][0], out[node.Right][1]
				for d := 0; d < 3; d++ {
					if rightLow[d] < low[d] { low[d] = rightLow[d] }
				}
				for d := 0; d < 3; d++ {
					if rightHigh[d] > high[d] { high[d] = rightHigh[d] }
				}
			}
			out[i] = [2][3]float64{ low, high }
		}
	}
}
