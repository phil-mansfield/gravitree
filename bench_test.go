package gravitree

import (
	"testing"
	"github.com/phil-mansfield/symtable"
)

func readFile(filename string) (points [][3]float64) {

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

func BenchmarkPlummer1e3(b *testing.B) {

	filename := "plummer.txt"
	x := readFile(filename)

	tree := NewTree(x)
	phi := make([]float64, len(x))

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		tree.Potential(0.0, phi)
	}
}