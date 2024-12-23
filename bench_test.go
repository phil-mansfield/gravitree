package gravitree

import (
	"testing"
)

func BenchmarkPlummer1e3(b *testing.B) {

	filename := "plummer.txt"
	x := readFile(filename)

	tree := NewTree(x)
	phi := make(Potential, len(x))

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		tree.Quantity(0.0, phi)
	}
}