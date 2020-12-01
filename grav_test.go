package gravitree

import (
	"math"
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

func TestRStats(t *testing.T) {
	tests := []struct{
		x [][3]float64
		rMax2, sigmaX2 float64
	}{
		{[][3]float64{{0, 0, 0}, {0, 0, 4}, {0, 0, 4}, {0, 0, 4}}, 9, 3},
		{[][3]float64{{0, 0, 0}, {0, 4, 0}, {0, 4, 0}, {0, 4, 0}}, 9, 3},
		{[][3]float64{{0, 0, 0}, {4, 0, 0}, {4, 0, 0}, {4, 0, 0}}, 9, 3},
		{[][3]float64{{1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
			{2, 1, 1}, {1, 2, 1}, {1, 1, 2}}, 1, 6/7.0},
	}

	for i := range tests {
		test := tests[i]
		tree := &KDTree{ Nodes: make([]KDNode, 1), Points: test.x }
		tree.Nodes[0].Start, tree.Nodes[0].End = 0, len(test.x)
		centers := [][3]float64{ centerOfMass(test.x) }
		
		rMax2, sigmaX2 := tree.rStats(centers, 0)

		if !almostEq(rMax2, test.rMax2, 1e-4) {
			t.Errorf("%d) Expected rMax2 = %.4f, got %.4f.",
				i+1, test.rMax2, rMax2)
		}
		if !almostEq(sigmaX2, test.sigmaX2, 1e-4) {
			t.Errorf("%d) Expected sigmaX2 = %.4f, got %.4f.",
				i+1, test.sigmaX2, sigmaX2)
		}
	}
}


func TestSquaredMonopoleDistance(t *testing.T) {
	relErr := 0.01
	tests := []struct{
		x [][3]float64
		monoDist float64
	}{
		{[][3]float64{{0, 0, 0}, {0, 0, 4}, {0, 0, 4}, {0, 0, 4}},
			18.88533865071371},
		{[][3]float64{{0, 0, 0}, {0, 4, 0}, {0, 4, 0}, {0, 4, 0}},
			18.88533865071371},
		{[][3]float64{{0, 0, 0}, {4, 0, 0}, {4, 0, 0}, {4, 0, 0}},
			18.88533865071371},
		{[][3]float64{{1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
			{2, 1, 1}, {1, 2, 1}, {1, 1, 2}}, 9.771692710302995},

	}

	for i := range tests {
		test := tests[i]
		tree := &KDTree{ Nodes: make([]KDNode, 1), Points: test.x }
		tree.Nodes[0].Start, tree.Nodes[0].End = 0, len(test.x)
		centers := [][3]float64{ centerOfMass(test.x) }

		out := []float64{ 0 }

		tree.SquaredMonopoleDistance(relErr, centers, out)

		if !almostEq(out[0], test.monoDist*test.monoDist, 1e-4) {
			t.Errorf("%d) expected mono distance = %.4f, got %.4f",
				i+1, test.monoDist*test.monoDist, out[0])
		}
	}
}

func TestPotential(t *testing.T) {
	x := [][3]float64{ {0, 0, 0} }
	nodes := []KDNode{ {-1, -1, -1, -1, 0, 1} }
	kt := &KDTree{ Points: x, Nodes: nodes }
	gt1 := NewGravitree(kt, 0.01)
	gt1.SMD[0] = 3

	x = [][3]float64{ {0, 0, -1}, {0, 0, 1} }
	nodes = []KDNode{
		{1, 2, -1, 0, 0, 2},
		{-1, -1, -1, -1, 1, 1},
		{-1, -1, -1, -1, 1, 2},
	}
	kt = &KDTree{ Points: x, Nodes: nodes }
	gt2right := NewGravitree(kt, 0.01)
	gt2right.SMD[0] = 3
	
	nodes = []KDNode{
		{1, 2, -1, 1, 0, 2},
		{-1, -1, -1, -1, 0, 1},
		{-1, -1, -1, -1, 2, 2},
	}
	kt = &KDTree{ Points: x, Nodes: nodes }
	gt2left := NewGravitree(kt, 0.01)
	gt2left.SMD[0] = 3
	
	tests := []struct{
		gt *Gravitree
		x [3]float64
		eps, phi float64
	}{
		{gt1, [3]float64{0, 0, 0}, 1.0, -1.0},
		{gt1, [3]float64{0, 0, 2}, 1.0, -1/math.Sqrt(5)},
		{gt1, [3]float64{0, 0, 4}, 1.0, -1/math.Sqrt(17)},
		{gt2right, [3]float64{0, 0, 0}, 1.0, -math.Sqrt(2)},
		{gt2right, [3]float64{0, 0, 4}, 1.0, -2/math.Sqrt(17)},
		{gt2left, [3]float64{0, 0, 0}, 1.0, -math.Sqrt(2)},
		{gt2left, [3]float64{0, 0, 4}, 1.0, -2/math.Sqrt(17)},
	}

	for i := range tests {
		test := tests[i]
		phi := test.gt.Potential(test.eps, test.x)
		if !almostEq(phi, test.phi, 1e-6) {
			t.Errorf("%d) Expected phi = %.6f, got %.6f.", i+1, test.phi, phi)
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

func almostEq(x1, x2, eps float64) bool {
	return !(x1 + eps < x2 || x1 - eps > x2)
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

func BenchmarkRStats_leaf1(b *testing.B) {
	benchmarkRStats(b, 1, 100000, -2)
}
func BenchmarkRStats_leaf2(b *testing.B) {
	benchmarkRStats(b, 2, 100000, -2)
}
func BenchmarkRStats_leaf4(b *testing.B) {
	benchmarkRStats(b, 4, 100000, -2)
}
func BenchmarkRStats_leaf8(b *testing.B) {
	benchmarkRStats(b, 8, 100000, -2)
}
func BenchmarkRStats_leaf16(b *testing.B) {
	benchmarkRStats(b, 16, 100000, -2)
}
func BenchmarkRStats_leaf32(b *testing.B) {
	benchmarkRStats(b, 32, 100000, -2)
}

func benchmarkRStats(
	b *testing.B, leafSize, points int, alpha float64,
) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	t := make([]*KDTree, nMax)
	centers := make([][][3]float64, nMax)
	
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range t {
		x := make([][3]float64, points)
		generatePowerLawPoints(alpha, 1.0, x)
		t[i] = NewKDTree(x, span, leafSize)
		centers[i] = make([][3]float64, len(t[i].Nodes))
		t[i].CenterOfMass(centers[i])
	}

	
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		tree := t[i % len(t)]
		for j := range tree.Nodes {
			tree.rStats(centers[i % len(t)], j)
		}
	}
}

func BenchmarkSquaredMonopoleDistance_leaf1(b *testing.B) {
	benchmarkSquaredMonopoleDistance(b, 1, 100000, -2)
}
func BenchmarkSquaredMonopoleDistance_leaf2(b *testing.B) {
	benchmarkSquaredMonopoleDistance(b, 2, 100000, -2)
}
func BenchmarkSquaredMonopoleDistance_leaf4(b *testing.B) {
	benchmarkSquaredMonopoleDistance(b, 4, 100000, -2)
}
func BenchmarkSquaredMonopoleDistance_leaf8(b *testing.B) {
	benchmarkSquaredMonopoleDistance(b, 8, 100000, -2)
}
func BenchmarkSquaredMonopoleDistance_leaf16(b *testing.B) {
	benchmarkSquaredMonopoleDistance(b, 16, 100000, -2)
}
func BenchmarkSquaredMonopoleDistance_leaf32(b *testing.B) {
	benchmarkSquaredMonopoleDistance(b, 32, 100000, -2)
}

func benchmarkSquaredMonopoleDistance(
	b *testing.B, leafSize, points int, alpha float64,
) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	t := make([]*KDTree, nMax)
	centers := make([][][3]float64, nMax)
	out := make([][]float64, nMax)
	
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range t {
		x := make([][3]float64, points)
		generatePowerLawPoints(alpha, 1.0, x)
		t[i] = NewKDTree(x, span, leafSize)
		centers[i] = make([][3]float64, len(t[i].Nodes))
		t[i].CenterOfMass(centers[i])
		out[i] = make([]float64, len(t[i].Nodes))
	}

	
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		tree := t[i % len(t)]
		tree.SquaredMonopoleDistance(
			0.01, centers[i % len(t)], out[i % len(t)],
		)
	}
}


func BenchmarkPotential_fid(b *testing.B) {
	benchmarkPotential(b, 4, 10000, -1, 0.1)
}

func benchmarkPotential(
	b *testing.B, leafSize, points int, alpha, relAcc float64,
) {	
	nMax := 100
	if nMax > b.N { nMax = b.N }
	gt := make([]*Gravitree, nMax)
	
	b.SetBytes(int64(points * 3 * 8))

	span := [2][3]float64{ {-1, -1, -1}, {1, 1, 1} }
	for i := range gt {
		x := make([][3]float64, points)
		generatePowerLawPoints(alpha, 1.0, x)
		kt := NewKDTree(x, span, leafSize)
		gt[i] = NewGravitree(kt, relAcc)
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		tree := gt[i % len(gt)]
		pts := tree.Tree.Points
		_ = pts
		
		for j := range pts {
			tree.Potential(0.01, pts[j])
		}
	}
}
