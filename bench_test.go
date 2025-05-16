package gravitree

import (	
	"testing"
	"github.com/phil-mansfield/symtable"
)

func readPointFile(filename string) [][3]float64{

	t := symtable.TextFile(filename)
	cols := t.ReadFloat64s([]int{0, 1, 2}) // column indices
	xs := cols[0]
	ys := cols[1]
	zs := cols[2]

	var result [][3]float64

	for i := 0; i < len(xs); i++ {
		var point [3]float64
		point[0] = xs[i]
		point[1] = ys[i]
		point[2] = zs[i]
		result = append(result, point)
	}

	return result
}

func benchmarkPotentialTree(b *testing.B, n int, filename string) {
	x := readPointFile(filename)

	tree := NewTree(x)
	phi := Potential(make([]float64, len(x)))

	b.SetBytes(int64(24 * n))
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		tree.Evaluate(0.0, phi)
	}
}

func benchmarkNewTree(b *testing.B, n int, filename string) {
	x := readPointFile(filename)

	b.SetBytes(int64(24 * n))
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		NewTree(x)
	}
}

func benchmarkBruteForce(b *testing.B, n int, filename string) {
	x := readPointFile(filename)
	phi := make([]float64, len(x))

	b.SetBytes(int64(24 * n))
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		BruteForcePotential(0.0, x, phi)
	}
}

func BenchmarkPotentialTree_1e2(b *testing.B) {
	benchmarkPotentialTree(b, int(1e2), "test_files/einasto_n=2_a=18.dat")
}
func BenchmarkPotentialTree_3e2(b *testing.B) {
	benchmarkPotentialTree(b, int(3e2), "test_files/einasto_n=2.5_a=18.dat")
}
func BenchmarkPotentialTree_1e3(b *testing.B) {
	benchmarkPotentialTree(b, int(1e3), "test_files/einasto_n=3_a=18.dat")
}
func BenchmarkPotentialTree_3e3(b *testing.B) {
	benchmarkPotentialTree(b, int(3e3), "test_files/einasto_n=3.5_a=18.dat")
}
func BenchmarkPotentialTree_1e4(b *testing.B) {
	benchmarkPotentialTree(b, int(1e4), "test_files/einasto_n=4_a=18.dat")
}
func BenchmarkPotentialTree_3e4(b *testing.B) {
	benchmarkPotentialTree(b, int(3e4), "test_files/einasto_n=4.5_a=18.dat")
}
func BenchmarkPotentialTree_1e5(b *testing.B) {
	benchmarkPotentialTree(b, int(1e5), "test_files/einasto_n=5_a=18.dat")
}

func BenchmarkNewTree_1e2(b *testing.B) {
	benchmarkNewTree(b, int(1e2), "test_files/einasto_n=2_a=18.dat")
}
func BenchmarkNewTree_3e2(b *testing.B) {
	benchmarkNewTree(b, int(3e2), "test_files/einasto_n=2.5_a=18.dat")
}
func BenchmarkNewTree_1e3(b *testing.B) {
	benchmarkNewTree(b, int(1e3), "test_files/einasto_n=3_a=18.dat")
}
func BenchmarkNewTree_3e3(b *testing.B) {
	benchmarkNewTree(b, int(3e3), "test_files/einasto_n=3.5_a=18.dat")
}
func BenchmarkNewTree_1e4(b *testing.B) {
	benchmarkNewTree(b, int(1e4), "test_files/einasto_n=4_a=18.dat")
}
func BenchmarkNewTree_3e4(b *testing.B) {
	benchmarkNewTree(b, int(3e4), "test_files/einasto_n=4.5_a=18.dat")
}
func BenchmarkNewTree_1e5(b *testing.B) {
	benchmarkNewTree(b, int(1e5), "test_files/einasto_n=5_a=18.dat")
}

func BenchmarkBruteForce_1e2(b *testing.B) {
	benchmarkBruteForce(b, int(1e2), "test_files/einasto_n=2_a=18.dat")
}
func BenchmarkBruteForce_3e2(b *testing.B) {
	benchmarkBruteForce(b, int(3e2), "test_files/einasto_n=2.5_a=18.dat")
}
func BenchmarkBruteForce_1e3(b *testing.B) {
	benchmarkBruteForce(b, int(1e3), "test_files/einasto_n=3_a=18.dat")
}
func BenchmarkBruteForce_3e3(b *testing.B) {
	benchmarkBruteForce(b, int(3e3), "test_files/einasto_n=3.5_a=18.dat")
}
func BenchmarkBruteForce_1e4(b *testing.B) {
	benchmarkBruteForce(b, int(1e4), "test_files/einasto_n=4_a=18.dat")
}
//func BenchmarkBruteForce_3e4(b *testing.B) {
//	benchmarkBruteForce(b, int(3e4), "test_files/einasto_n=4.5_a=18.dat")
//}
//func BenchmarkBruteForce_1e5(b *testing.B) {
//	benchmarkBruteForce(b, int(1e5), "test_files/einasto_n=5_a=18.dat")
//}
