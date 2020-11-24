package gravitree

import (
	"testing"
	"math/rand"
)

func TestSort3(t *testing.T) {
	x0 := [3]float64{ 0, 0, 3 }
	x1 := [3]float64{ 0, 0, 2 }
	x2 := [3]float64{ 0, 0, 1 }
	tests := []struct {
		a, b, c [3]float64
	}{
		{x0, x1, x2},
		{x1, x0, x2},
		{x0, x2, x1},
		{x2, x0, x1},
		{x1, x2, x0},
		{x2, x1, x0},
	}

	for i := range tests {
		test := tests[i]
		y0, y1, y2 := sort3(test.a, test.b, test.c, 2)
		if x0 != y0 || x1 != y1 || x2 != y2 {
			t.Errorf("(%.1f %.1f %.1f) sorted to (%.1f %.1f %.1f))",
				test.a, test.b, test.c, y0, y1, y2)
		}
	}
}

func TestPartition(t *testing.T) {
	tests := []struct{
		x [][3]float64
		pivot float64
	}{
		{[][3]float64{}, 0},
		{[][3]float64{{0,0,1}}, 1},
		{[][3]float64{{0,0,2}, {0,0,1}}, 2},
		{[][3]float64{{0,0,1}, {0,0,2}}, 2},
		{[][3]float64{{0,0,1}, {0,0,2}, {0,0,3}}, 2},
		{[][3]float64{{0,0,3}, {0,0,2}, {0,0,1}}, 2},
		{[][3]float64{{0,0,2}, {0,0,2}, {0,0,2}}, 2},
		{[][3]float64{{0,0,2}, {0,0,2}, {0,0,2}, {0,0,2}}, 2},
		{[][3]float64{{0,0,2}, {0,0,1}, {0,0,2}, {0,0,2}}, 2},
		{[][3]float64{{0,0,2}, {0,0,3}, {0,0,2}, {0,0,2}}, 2},
		{[][3]float64{{0,0,1}, {0,0,2}, {0,0,3}, {0,0,4}, {0,0,5}}, 3},
		{[][3]float64{{0,0,1}, {0,0,4}, {0,0,3}, {0,0,2}, {0,0,5}}, 3},
		{[][3]float64{{0,0,1}, {0,0,4}, {0,0,3}, {0,0,2}, {0,0,5}}, 3},
		{[][3]float64{{0,0,3}, {0,0,4}, {0,0,5}, {0,0,2}, {0,0,1}}, 3},
	}
	d := 2

	for i := range tests {
		test := tests[i]

		buf := make([][3]float64, len(test.x))
		copy(buf, test.x)

		iPiv := partition(buf, d)
		if len(buf) == 0 { continue }

		if test.pivot != buf[iPiv][d] {
			t.Errorf("%d) Expected pivot = %.1f, got %.1f",
				1+i, test.pivot, buf[iPiv])
		} else if !isPartitioned(buf, iPiv, d) {
			t.Errorf("%d) %.1f is not pivoted at index %d, dimension %d.",
				1+i, buf, iPiv, d)
		}
	}

	rand.Seed(0)
	for i := 0; i < 1000; i++ {
		vec := make([][3]float64, 1000)
		for i := range vec { vec[i][2] = rand.Float64() }

		iPiv := partition(vec, 2)
		if !isPartitioned(vec, iPiv, 2) {
			t.Errorf("Random vector %d icorrectly partitioned.", i)
		}
	}
}

func isPartitioned(x [][3]float64, iPiv, d int) bool {
	pivot := x[iPiv][d]
	for i := 0; i < iPiv; i++ {
		if x[i][d] > pivot { return false }
	}
	for i := iPiv + 1; i < len(x); i++ {
		if x[i][d] < pivot { return false }
	}
	return true
}

func TestKDTreeConsistency(t *testing.T) {
	tubeSpan := [2][3]float64{{-0.5, -0.5, 0}, {0.5, 0.5, 80}}
	cubeSpan := [2][3]float64{{0, 0, 0}, {10, 10, 10}}
	tests := []struct{
		points [][3]float64
		span [2][3]float64
		leafSize int
	}{
		{[][3]float64{{0,0,32}}, tubeSpan, 1},
		{[][3]float64{{0,0,32}, {0,0,1}, {0,0,42}}, tubeSpan, 3},
		{[][3]float64{{0,0,32}, {0,0,1}, {0,0,42}}, tubeSpan, 1},
		{[][3]float64{{0,0,42}, {0,0,1}, {0,0,32}}, tubeSpan, 1},
		{[][3]float64{{0,0,1}, {0,0,32}, {0,0,42}}, tubeSpan, 1},
		{[][3]float64{{0,0,32}, {0,0,42}, {0,0,1}}, tubeSpan, 1},
		{[][3]float64{{0,0,42}, {0,0,32}, {0,0,1}}, tubeSpan, 1},
		{[][3]float64{{0,0,32}, {0,0,32}, {0,0,32}}, tubeSpan, 1},
		{[][3]float64{{0,0,0}, {0,0,10}, {0,0,20}, {0,0,30}, {0,0,40},
			{0,0,50}, {0,0,60}, {0,0,70}, {0,0,80}}, tubeSpan, 1},
		{[][3]float64{{0,0,0}, {0,0,10}, {0,0,20}, {0,0,30}, {0,0,40},
			{0,0,50}, {0,0,60}, {0,0,70}, {0,0,80}}, tubeSpan, 2},
		{[][3]float64{{0,0,80}, {0,0,70}, {0,0,60}, {0,0,50}, {0,0,40},
			{0,0,30}, {0,0,20}, {0,0,10}, {0,0,0}}, tubeSpan, 1},
		{[][3]float64{{0,0,0}, {1,1,2}, {2, 2, 4}}, cubeSpan, 1},
		{[][3]float64{{5,5,5}, {0,0,0}, {10,10,10}, {1,1,2}, {2,3,8},
			{7,4,1}, {3,9,3}, {4,8,4}, {6,2,9}, {8,4,6}, {9,6,7}},
			cubeSpan, 1},
		{[][3]float64{{2,3,9}, {1,4,3}, {3,7,10}, {9,7,7}}, cubeSpan, 1},
	}
	_, _ = cubeSpan, tubeSpan

	for i := range tests {
		test := tests[i]
		tree := NewKDTree(test.points, test.span, test.leafSize)
		
		if !testKDTreeCompleteness(tree) {
			t.Errorf("%d) Tree nodes are not complete.", i+1)
		}
		if !testKDNodeSpanConsistency(tree, 0, test.span) {
			t.Errorf("%d) Tree spans are not consistent.", i+1)
		}
		if !testKDNodePointsInSpan(tree, 0, test.span) {
			t.Errorf("%d) Node has points outside span", i+1)
		}
	}

	rand.Seed(0)
	for i := 0; i < 1000; i++ {
		x := make([][3]float64, 1000)
		for j := range x {
			x[j] = [3]float64{
				rand.Float64()*10, rand.Float64()*10, rand.Float64()*10,
			}
		}

		tree := NewKDTree(x, cubeSpan, 1)
		
		if !testKDTreeCompleteness(tree) {
			t.Errorf("rand %3d) Tree nodes are not complete.", i+1)
		}
		if !testKDNodeSpanConsistency(tree, 0, cubeSpan) {
			t.Errorf("rand %3d) Tree spans are not consistent.", i+1)
		}
		if !testKDNodePointsInSpan(tree, 0, cubeSpan) {
			t.Errorf("rand %3d) Node has points outside span", i+1)
		}
		break
	}
}


func testKDTreeCompleteness(t *KDTree) bool {
	ok := make([]bool, len(t.Nodes))
	if !testKDNodeCompleteness(t, 0, ok) { return false }
	
	for i := range ok {
		if !ok[i] { return false }
	}

	return true
}

func testKDNodeCompleteness(t *KDTree, i int, ok []bool) bool {
	if i == -1 { return true }
	node := &t.Nodes[i]
	if ok[i] { return false }
	ok[i] = true
	if !testKDNodeCompleteness(t, node.Left, ok) { return false }
	if !testKDNodeCompleteness(t, node.Right, ok) { return false }

	return true
}

func testKDNodeSpanConsistency(t *KDTree, i int, span [2][3]float64) bool {
	if i == -1 { return true }
	for d := 0; d < 3; d++ {
		if span[0][d] >= span[1][d] { return false }
	}

	node := &t.Nodes[i]
	if node.Left == -1 && node.Right == -1 { return true }
	
	left, right := t.splitSpan(node, span)
	if left[0] != span[0] || right[1] != span[1] { return false }
	if left[1][node.Dim] != right[0][node.Dim] { return false }

	return testKDNodeSpanConsistency(t, node.Left, left) &&
		testKDNodeSpanConsistency(t, node.Right, right)
}

func testKDNodePointsInSpan(t *KDTree, i int, span [2][3]float64) bool {
	if i == -1 { return true }
	
	node := &t.Nodes[i]

	for i := node.Start; i < node.End; i++ {
		for d := 0; d < 3; d++ {
			if t.Points[i][d] < span[0][d] || t.Points[i][d] > span[1][d] {
				return false
			}
		}
	}

	if node.Left != -1 && node.Right != -1 {
		left, right := t.splitSpan(node, span)
		return testKDNodePointsInSpan(t, node.Left, left) &&
			testKDNodePointsInSpan(t, node.Right, right)
	}
	return true
}
