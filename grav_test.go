package gravitree

import (
	"testing"
)

func TestCenterOfMass(t *testing.T) {
	x := make([][3]float64, 1000)
	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	generatePowerLawPoints(-2, 1, x)

	leafSizes := []int{ 1, 10 }
	for j := range leafSizes {
		tree := NewKDTree(x, span, leafSizes[j])
		centers := make([][3]float64, len(tree.Nodes))
		tree.CenterOfMass(centers)
		
		for i := range tree.Nodes {
			start, end := tree.Nodes[i].Start, tree.Nodes[i].End
			c := centerOfMass(tree.Points[start: end])

			if !vecAlmostEq(c, centers[i], 1e-4) {
				t.Errorf("%d) Got center of mass %.4f, epxected %.4f",
					i, c, centers[i])
			}
		}
	}
}

func TestSpan(t *testing.T) {
	x := make([][3]float64, 2)
	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	generatePowerLawPoints(-2, 1, x)

	leafSizes := []int{ 1, 10 }
	for j := range leafSizes {
		tree := NewKDTree(x, span, leafSizes[j])
		nodeSpan := make([][2][3]float64, len(tree.Nodes))
		tree.Span(nodeSpan)
		
		for i := range tree.Nodes {
			start, end := tree.Nodes[i].Start, tree.Nodes[i].End
			span := naiveSpan(tree.Points[start: end])
			low, high := span[0], span[1]

			if !vecAlmostEq(nodeSpan[i][0], low, 1e-4) {
				t.Errorf("%d) Got low %.4f, epxected %.4f",
					i, nodeSpan[i][0], low)
			}
			if !vecAlmostEq(nodeSpan[i][1], high, 1e-4) {
				t.Errorf("%d) Got low %.4f, epxected %.4f",
					i, nodeSpan[i][1], high)
			}
		}
	}
}

func vecAlmostEq(x1, x2 [3]float64, eps float64) bool {
	for d := 0; d < 3; d++ {
		if x1[d] + eps < x2[d] || x1[d] - eps > x2[d] {
			return false
		}
	}
	return true
}

func centerOfMass(x [][3]float64) [3]float64 {
	c := [3]float64{ }
	for i := range x {
		for d := 0; d < 3; d++ {
			c[d] += x[i][d]
		}
	}
	for d := 0; d < 3; d++ {
		c[d] /= float64(len(x))
	}

	return c
}

func naiveSpan(x [][3]float64) [2][3]float64 {
	if len(x) == 0 { return [2][3]float64{ } }
	
	low, high := x[0], x[0]
	for i := 1; i < len(x); i++ {
		pt := x[i]
		for d := 0; d < 3; d++ {
			if pt[d] < low[d] { low[d] = pt[d] }
			if pt[d] > high[d] { high[d] = pt[d] }
		}
	}
	return [2][3]float64{ low, high }
}

func BenchmarkCenterOfMass_leaf1(b *testing.B) {
	benchmarkCenterOfMass(b, 1, 100000, -2)
}
func BenchmarkCenterOfMass_leaf2(b *testing.B) {
	benchmarkCenterOfMass(b, 2, 100000, -2)
}
func BenchmarkCenterOfMass_leaf4(b *testing.B) {
	benchmarkCenterOfMass(b, 4, 100000, -2)
}
func BenchmarkCenterOfMass_leaf8(b *testing.B) {
	benchmarkCenterOfMass(b, 8, 100000, -2)
}
func BenchmarkCenterOfMass_leaf16(b *testing.B) {
	benchmarkCenterOfMass(b, 16, 100000, -2)
}
func BenchmarkCenterOfMass_leaf32(b *testing.B) {
	benchmarkCenterOfMass(b, 32, 100000, -2)
}

func benchmarkCenterOfMass(
	b *testing.B, leafSize, points int, alpha float64,
) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	t := make([]*KDTree, nMax)
	out := make([][][3]float64, nMax)
	
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range t {
		x := make([][3]float64, points)
		generatePowerLawPoints(alpha, 1.0, x)
		t[i] = NewKDTree(x, span, leafSize)
		out[i] = make([][3]float64, len(t[i].Nodes))
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		tree := t[i % len(t)]
		tree.CenterOfMass(out[i % len(t)])
	}
}

func BenchmarkSpan_leaf1(b *testing.B) {
	benchmarkSpan(b, 1, 100000, -2)
}
func BenchmarkSpan_leaf2(b *testing.B) {
	benchmarkSpan(b, 2, 100000, -2)
}
func BenchmarkSpan_leaf4(b *testing.B) {
	benchmarkSpan(b, 4, 100000, -2)
}
func BenchmarkSpan_leaf8(b *testing.B) {
	benchmarkSpan(b, 8, 100000, -2)
}
func BenchmarkSpan_leaf16(b *testing.B) {
	benchmarkSpan(b, 16, 100000, -2)
}
func BenchmarkSpan_leaf32(b *testing.B) {
	benchmarkSpan(b, 32, 100000, -2)
}

func benchmarkSpan(
	b *testing.B, leafSize, points int, alpha float64,
) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	t := make([]*KDTree, nMax)
	out := make([][][2][3]float64, nMax)
	
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range t {
		x := make([][3]float64, points)
		generatePowerLawPoints(alpha, 1.0, x)
		t[i] = NewKDTree(x, span, leafSize)
		out[i] = make([][2][3]float64, len(t[i].Nodes))
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		tree := t[i % len(t)]
		tree.Span(out[i % len(t)])
	}
}
