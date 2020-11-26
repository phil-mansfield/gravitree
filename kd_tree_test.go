package gravitree

import (
	"fmt"
	"testing"
	"math/rand"
)

/////////////
// Testing //
/////////////

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
		{[][3]float64{{0,0,2}, {0,0,1}}, 1},
		{[][3]float64{{0,0,1}, {0,0,2}}, 1},
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
				1+i, test.pivot, buf[iPiv][d])
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
		/*
		{[][3]float64{{0,0,32}}, tubeSpan, 1},
		{[][3]float64{{0,0,32}, {0,0,1}, {0,0,42}}, tubeSpan, 3},
		{[][3]float64{{0,0,32}, {0,0,1}, {0,0,42}}, tubeSpan, 1},
		{[][3]float64{{0,0,42}, {0,0,1}, {0,0,32}}, tubeSpan, 1},
		{[][3]float64{{0,0,1}, {0,0,32}, {0,0,42}}, tubeSpan, 1},
		{[][3]float64{{0,0,32}, {0,0,42}, {0,0,1}}, tubeSpan, 1},
		{[][3]float64{{0,0,42}, {0,0,32}, {0,0,1}}, tubeSpan, 1},
		{[][3]float64{{0,0,32}, {0,0,32}, {0,0,32}}, tubeSpan, 1},
		*/
		{[][3]float64{{0,0,0}, {0,0,10}, {0,0,20}, {0,0,30}, {0,0,40},
			{0,0,50}, {0,0,60}, {0,0,70}, {0,0,80}}, tubeSpan, 1},
		/*
		{[][3]float64{{0,0,0}, {0,0,10}, {0,0,20}, {0,0,30}, {0,0,40},
			{0,0,50}, {0,0,60}, {0,0,70}, {0,0,80}}, tubeSpan, 2},
		{[][3]float64{{0,0,80}, {0,0,70}, {0,0,60}, {0,0,50}, {0,0,40},
			{0,0,30}, {0,0,20}, {0,0,10}, {0,0,0}}, tubeSpan, 1},
		{[][3]float64{{0,0,0}, {1,1,2}, {2, 2, 4}}, cubeSpan, 1},
		{[][3]float64{{5,5,5}, {0,0,0}, {10,10,10}, {1,1,2}, {2,3,8},
			{7,4,1}, {3,9,3}, {4,8,4}, {6,2,9}, {8,4,6}, {9,6,7}},
			cubeSpan, 1},
		{[][3]float64{{5,5,5}, {0,0,0}, {10,10,10}, {1,1,2}, {2,3,8},
			{7,4,1}, {3,9,3}, {4,8,4}, {6,2,9}, {8,4,6}, {9,6,7}},
			cubeSpan, 1},
		{[][3]float64{{2,3,9}, {1,4,3}, {3,7,10}, {9,7,7}}, cubeSpan, 1},
*/
	}
	_, _ = cubeSpan, tubeSpan

	for i := range tests {
		test := tests[i]
		tree := NewKDTree(test.points, test.span, test.leafSize)
		
		if !testKDTreeCompleteness(tree) {
			t.Errorf("%d) Tree nodes are not complete.", i+1)
		}

		missing, repeat := testKDTreePointCompleteness(tree)
		if missing > 0 {
			t.Errorf("%d) Tree missing %d points.", i+1, missing)
		}
		if repeat > 0 {
			t.Errorf("%d) Tree repeating %d points.", i+1, repeat)
		}
		
		if !testKDNodeSpanConsistency(tree, 0, test.span) {
			t.Errorf("%d) Tree spans are not consistent.", i+1)
		}
		if !testKDNodePointsInSpan(tree, 0, test.span) {
			t.Errorf("%d) Node has points outside span", i+1)
		}
	}
	
	for _, leafSize := range []int{ 1, 2, 10, 100 } {
		rand.Seed(0)
		for i := 0; i < 1000; i++ {
			x := make([][3]float64, 1000)
			for j := range x {
				x[j] = [3]float64{
					rand.Float64()*10, rand.Float64()*10, rand.Float64()*10,
				}
			}
			
			tree := NewKDTree(x, cubeSpan, leafSize)
			
			if !testKDTreeCompleteness(tree) {
				t.Errorf("rand %3d leaf %3d) Tree nodes are not complete.",
					i+1, leafSize)
			}
			if !testKDNodeSpanConsistency(tree, 0, cubeSpan) {
				t.Errorf("rand %3d leaf %3d) Tree spans are not consistent.",
					i+1, leafSize)
			}
			if !testKDNodePointsInSpan(tree, 0, cubeSpan) {
				t.Errorf("rand %3d leaf %3d) Node has points outside span",
					i+1, leafSize)
			}
			break
		}
	}
}

func testKDTreePointCompleteness(t *KDTree) (missing, repeat int) {
	count := make([]int, len(t.Points))
	testKDNodePointCompleteness(t, 0, count)
	
	missing, repeat = 0, 0
	for i := range count {
		if count[i] == 0 { missing++ }
		if count[i] > 1 { repeat++ }
	}

	return missing, repeat
}

func testKDNodePointCompleteness(t *KDTree, i int, count []int) {
	node := &t.Nodes[i]
	if node.Left == -1 {
		for i := node.Start; i < node.End; i++ {
			count[i]++
		}
	} else {
		count[node.Pivot]++
		testKDNodePointCompleteness(t, node.Left, count)
		testKDNodePointCompleteness(t, node.Right, count)
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
		if span[0][d] > span[1][d] {
			fmt.Printf("%d) Incorrect span ordering, %.0f\n", i, span)
			return false
		}
	}

	node := &t.Nodes[i]
	if node.Left == -1 && node.Right == -1 { return true }
	
	left, right := t.splitSpan(node, span)
	if left[0] != span[0] || right[1] != span[1] {
		fmt.Printf("%d) left, %.0f, and right, %.0f, sub-spans don't align " +
			"with bas span, %0f.\n", i, left, right, span)
		return false
	}
	if left[1][node.Dim] != right[0][node.Dim] {
		fmt.Printf("%d) left, %.0f, and right, %.0f, sub-spans don't align " +
			"at node dim, %d.\n", i, left, right, node.Dim)
		return false
	}

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

//////////////////
// Benchmarking //
//////////////////

func BenchmarkKDTreeUniform_1e3(b *testing.B) {
	benchmarkKDTreeUniform(b, 4, 1000)
}
func BenchmarkKDTreeUniform_n1e4(b *testing.B) {
	benchmarkKDTreeUniform(b, 4, 10000)
}
func BenchmarkKDTreeUniform_n1e5(b *testing.B) {
	benchmarkKDTreeUniform(b, 4, 100000)
}
func BenchmarkKDTreeUniform_n1e6(b *testing.B) {
	benchmarkKDTreeUniform(b, 4, 1000000)
}
func BenchmarkKDTreeUniform_n1e7(b *testing.B) {
	benchmarkKDTreeUniform(b, 4, 10000000)
}

func BenchmarkKDTreeUniform_leaf1(b *testing.B) {
	benchmarkKDTreeUniform(b, 1, 100000)
}
func BenchmarkKDTreeUniform_leaf2(b *testing.B) {
	benchmarkKDTreeUniform(b, 2, 100000)
}
func BenchmarkKDTreeUniform_leaf4(b *testing.B) {
	benchmarkKDTreeUniform(b, 4, 100000)
}
func BenchmarkKDTreeUniform_leaf8(b *testing.B) {
	benchmarkKDTreeUniform(b, 8, 100000)
}
func BenchmarkKDTreeUniform_leaf16(b *testing.B) {
	benchmarkKDTreeUniform(b, 16, 100000)
}
func BenchmarkKDTreeUniform_leaf32(b *testing.B) {
	benchmarkKDTreeUniform(b, 32, 100000)
}

func BenchmarkKDTreeUniform_alpha0_0(b *testing.B) {
	benchmarkKDTreePowerLaw(b, 4, 100000, 0.0)
}
func BenchmarkKDTreeUniform_alpha0_5(b *testing.B) {
	benchmarkKDTreePowerLaw(b, 4, 100000, -0.5)
}
func BenchmarkKDTreeUniform_alpha1_0(b *testing.B) {
	benchmarkKDTreePowerLaw(b, 4, 100000, -1.0)
}
func BenchmarkKDTreeUniform_alpha1_5(b *testing.B) {
	benchmarkKDTreePowerLaw(b, 4, 100000, -1.5)
}
func BenchmarkKDTreeUniform_alpha2_0(b *testing.B) {
	benchmarkKDTreePowerLaw(b, 4, 100000, -2.0)
}
func BenchmarkKDTreeUniform_alpha2_5(b *testing.B) {
	benchmarkKDTreePowerLaw(b, 4, 100000, -2.5)
}


func benchmarkKDTreeUniform(b *testing.B, leafSize, points int) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	x := make([][][3]float64, nMax)
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range x {
		x[i] = make([][3]float64, points)
		generateUniformPoints(span, x[i])
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		NewKDTree(x[i % len(x)], span, leafSize)
	}
}

func benchmarkKDTreePowerLaw(
	b *testing.B, leafSize, points int, alpha float64,
) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	x := make([][][3]float64, nMax)
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range x {
		x[i] = make([][3]float64, points)
		generatePowerLawPoints(alpha, 1.0, x[i])
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		NewKDTree(x[i % len(x)], span, leafSize)
	}
}
