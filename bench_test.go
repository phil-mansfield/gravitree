package gravitree

import (
	"testing"
	"github.com/phil-mansfield/symtable"
)

// readPointFile reads the points 
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

func benchmarkEinastoTree(b *testing.B, n int, filename string) {
	x := readPointFile(filename)

	tree := NewTree(x)
	phi := Potential(make([]float64, len(x)))

	b.SetBytes(int64(12 * n))
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		tree.Quantity(0.0, phi)
	}
}

func benchmarkEinastoOldPotential(b *testing.B, n int, filename string) {
	x := readPointFile(filename)

	tree := NewTree(x)
	phi := make([]float64, len(x))

	b.SetBytes(int64(12 * n))
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		tree.oldPotential(0.0, phi)
	}
}

func BenchmarkEinastoOldPotential_1e2(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(1e2), "test_files/einasto_n=2_a=18.dat")
}
func BenchmarkEinastoOldPotential_3e2(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(3e2), "test_files/einasto_n=2.5_a=18.dat")
}
func BenchmarkEinastoOldPotential_1e3(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(1e3), "test_files/einasto_n=3_a=18.dat")
}
func BenchmarkEinastoOldPotential_3e3(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(3e3), "test_files/einasto_n=3.5_a=18.dat")
}
func BenchmarkEinastoOldPotential_1e4(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(1e4), "test_files/einasto_n=4_a=18.dat")
}
func BenchmarkEinastoOldPotential_3e4(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(3e4), "test_files/einasto_n=4.5_a=18.dat")
}
func BenchmarkEinastoOldPotential_1e5(b *testing.B) {
	benchmarkEinastoOldPotential(b, int(1e5), "test_files/einasto_n=5_a=18.dat")
}

func BenchmarkEinastoTree_1e2(b *testing.B) {
	benchmarkEinastoTree(b, int(1e2), "test_files/einasto_n=2_a=18.dat")
}
func BenchmarkEinastoTree_3e2(b *testing.B) {
	benchmarkEinastoTree(b, int(3e2), "test_files/einasto_n=2.5_a=18.dat")
}
func BenchmarkEinastoTree_1e3(b *testing.B) {
	benchmarkEinastoTree(b, int(1e3), "test_files/einasto_n=3_a=18.dat")
}
func BenchmarkEinastoTree_3e3(b *testing.B) {
	benchmarkEinastoTree(b, int(3e3), "test_files/einasto_n=3.5_a=18.dat")
}
func BenchmarkEinastoTree_1e4(b *testing.B) {
	benchmarkEinastoTree(b, int(1e4), "test_files/einasto_n=4_a=18.dat")
}
func BenchmarkEinastoTree_3e4(b *testing.B) {
	benchmarkEinastoTree(b, int(3e4), "test_files/einasto_n=4.5_a=18.dat")
}
func BenchmarkEinastoTree_1e5(b *testing.B) {
	benchmarkEinastoTree(b, int(1e5), "test_files/einasto_n=5_a=18.dat")
}