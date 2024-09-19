package gravitree

import (
	"math"
)

func (t *Tree) SearchSphere(x [3]float64, r float64, buf... []int) []int {
	var out []int
	if len(buf) > 0 {
		out = buf[0]
	} else {
		out = []int{ }
	}
	out = out[:0]

	return t.walkTreeSearchSphere(0, x, r, out)
}

func (t *Tree) walkTreeSearchSphere(
	i int, x [3]float64, r float64, buf []int,
) []int {
	n := &t.Nodes[i]
	// the first three conditionals are base cases, the last one is the
	// recursive case.

	// This square root is painful, but I don't know how to avoid it, tbh
	if dr := math.Sqrt(calcR2(&x, &n.Center)); dr - n.RMax > r {
		// Node and sphere are disjoint
		return buf
	} else if dr + n.RMax < r {
		// Node completely contained within sphere: add everything
		for i := n.Start; i < n.End; i++ {
			buf = append(buf, t.Index[i])
		}
		return buf
	} else if t.Nodes[i].Left == -1 {
		r2 := r*r
		// Leaf node, do a brue force search
		for i := n.Start; i < n.End; i++ {
			dr2 := calcR2(&x, &t.Points[i])
			if r2 > dr2 {
				buf = append(buf, t.Index[i])
			}
		}
		return buf
	} else {
		buf = t.walkTreeSearchSphere(n.Left, x, r, buf)
		buf = t.walkTreeSearchSphere(n.Right, x, r, buf)
		return buf
	}
}

func calcR2(x1, x2 *[3]float64) float64 {
	dx := x1[0] - x2[0]
	dy := x1[1] - x2[1]
	dz := x1[2] - x2[2]
	return dx*dx + dy*dy + dz*dz
}