package gravitree

import (
	"math"
	"testing"
)

func TestPartition(t *testing.T) {
	tests := []struct{
		x []float64
		pivot float64
		left int
	}{
		{[]float64{}, 0,  0},

		{[]float64{0}, 0,  1},
		{[]float64{0}, 1,  1},
		{[]float64{0}, -1,  0},

		{[]float64{0, 2}, -1, 0},
		{[]float64{0, 2}, 1, 1},
		{[]float64{0, 2}, 3, 2},
		{[]float64{2, 0}, -1, 0},
		{[]float64{2, 0}, 1, 1},
		{[]float64{2, 0}, 3, 2},

		{[]float64{0, 2, 4}, -1, 0},
		{[]float64{0, 2, 4}, 1, 1},
		{[]float64{0, 2, 4}, 3, 2},
		{[]float64{0, 2, 4}, 5, 3},

		{[]float64{0, 4, 2}, -1, 0},
		{[]float64{0, 4, 2}, 1, 1},
		{[]float64{0, 4, 2}, 3, 2},
		{[]float64{0, 4, 2}, 5, 3},

		{[]float64{2, 0, 4}, -1, 0},
		{[]float64{2, 0, 4}, 1, 1},
		{[]float64{2, 0, 4}, 3, 2},
		{[]float64{2, 0, 4}, 5, 3},
		
		{[]float64{4, 0, 2}, -1, 0},
		{[]float64{4, 0, 2}, 1, 1},
		{[]float64{4, 0, 2}, 3, 2},
		{[]float64{4, 0, 2}, 5, 3},

		{[]float64{2, 4, 0}, -1, 0},
		{[]float64{2, 4, 0}, 1, 1},
		{[]float64{2, 4, 0}, 3, 2},
		{[]float64{2, 4, 0}, 5, 3},

		{[]float64{4, 2, 0}, -1, 0},
		{[]float64{4, 2, 0}, 1, 1},
		{[]float64{4, 2, 0}, 3, 2},
		{[]float64{4, 2, 0}, 5, 3},
	}
	
	for i := range tests {
		test := tests[i]

		// Create the input arrays.
		vec := make([][3]float64, len(test.x))
		for i := range vec { vec[i][0] = test.x[i] }
		idx := make([]int, len(test.x))
		for i := range idx { idx[i] = i }

		// Perform partition.
		left := partition(vec,  idx, 0, test.pivot)
		x := make([]float64, len(test.x))
		for i := range x { x[i] = vec[i][0] }

		// Run tests.
		if left != test.left {
			t.Errorf("%d) expected left = %d for x = %.0f, pivot = %.0f, " +
				"got %d.", i, test.left, test.x, test.pivot, left)
		} else if !isPartitioned(x, test.pivot, left) {
			t.Errorf("%d) post partition() array, %.0f is not partitioned " +
				"below index %d at %.0f", i, x, left, test.pivot)
		} else if !isIdxPartitioned(test.x, idx, test.pivot, left) {
			t.Errorf("%d) post partition() index array, %d is not " +
				"partitioned below index %d at %.0f", i, idx, left, test.pivot)
		}
	}
}

func isPartitioned(x []float64, pivot float64, left int) bool {
	for i := 0; i < left; i++ {
		if x[i] > pivot { return false }
	}
	for i := left; i < len(x); i++ {
		if x[i] <= pivot { return false }
	}

	return true
}

func isIdxPartitioned(
	xOrig []float64, idx []int, pivot float64, left int,
) bool {
	for i := 0; i < left; i++ {
		if xOrig[idx[i]] > pivot { return false }
	}
	for i := left; i < len(xOrig); i++ {
		if xOrig[idx[i]] <= pivot { return false }
	}

	return true
}

func TestPointSpan(t *testing.T) {
	tests := []struct {
		x [][3]float64
		span [2][3]float64
	}{
		{[][3]float64{}, [2][3]float64{{}, {}}},
		{[][3]float64{{1, 1, 1}}, [2][3]float64{{1, 1, 1}, {1, 1, 1}}},
		{[][3]float64{{1, 1, 1}, {2, 2, 2}},
			[2][3]float64{{1, 1, 1}, {2, 2, 2}}},
		{[][3]float64{{2, 2, 2}, {1, 1, 1}},
			[2][3]float64{{1, 1, 1}, {2, 2, 2}}},
		{[][3]float64{{4, 1, 10}, {2, 3, 11}},
			[2][3]float64{{2, 1, 10}, {4, 3, 11}}},
		{[][3]float64{{4, 1, 10}, {2, 3, 11}, {3, 2, 12}},
			[2][3]float64{{2, 1, 10}, {4, 3, 12}}},
	}

	for i := range tests {
		test := tests[i]
		span := pointSpan(test.x)
		if span != test.span {
			t.Errorf("%d) Expected span %.0f from points %.0f, got %.0f",
				i, test.span, test.x, span) 
		}
	}
}

func TestChooseNodeDimension(t *testing.T) {
	tests := []struct{
		span [2][3]float64
		dim int
	}{
		{[2][3]float64{{0, 0, 0}, {1, 1, 1}}, 1},
		{[2][3]float64{{0, 0, 0}, {1, 0.5, 1}}, 2},
		{[2][3]float64{{0, 0, 0}, {1, 0.5, 0.5}}, 0},
		{[2][3]float64{{0, 0, 0}, {0.5, 0.5, 1}}, 2},
		{[2][3]float64{{0, 0, 0}, {0.5, 1, 0.5}}, 1},
		{[2][3]float64{{-1, -0.5, -0.5}, {0, 0, 0}}, 0},
		{[2][3]float64{{-0.5, -0.5, -1}, {0, 0, 0}}, 2},
		{[2][3]float64{{-0.5, -1, -0.5}, {0, 0, 0}}, 1},
	}

	for i := range tests {
		test := tests[i]
		dim := chooseNodeDimension(test.span)
		if dim != test.dim {
			t.Errorf("%d) Expected dim = %d from chooseNodeDimension on " +
				"span %.1f, got %d", i, test.dim, test.span, dim)
		}
	}
}

func TestROpen2(t *testing.T) {
	pts := [][3]float64{ {0, 0, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}
	span := [2][3]float64{ {0, 0, 0}, {2, 1, 0} }
	tests := []struct {
		criteria OpeningCriteria
		theta float64
		rOpen float64
	}{
		{BarnesHut, 0.5, 4.0},
		{PKDGRAV3, 0.5, 3.3541019662496847},
		{SalmonWarren, 0.5, 2.23606797749979},
	}

	tree := &Tree{ }
	tree.Points = pts
	tree.Nodes = make([]Node, 1)
	node := &tree.Nodes[0]
	node.Start, node.End = 0, 4

	for i := range tests {
		test := tests[i]

		tree.Criteria = test.criteria
		tree.Theta = test.theta

		rOpen2 := tree.ROpen2(0, span)
		rOpen := math.Sqrt(rOpen2)

		if !almostEq(rOpen, test.rOpen, 1e-5) {
			t.Errorf("%d) Expected r_open = %.6f, got %.6f",
				i, test.rOpen, rOpen)
		}
	}
}

func almostEq(x, y, eps float64) bool {
	return x + eps > y && x - eps < y
}
