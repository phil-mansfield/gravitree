package gravitree

import (
	"testing"
	"math/rand"
	"math"
)

func TestSearchSphere(t *testing.T) {
	rand.Seed(0)

	n := 1000
	x := make([][3]float64, n)
	for i := range x {
		x[i][0] = rand.Float64()
		x[i][1] = rand.Float64()
		x[i][2] = rand.Float64()
	}

	tree := NewTree(x)

	pt := [3]float64{0.5, 0.5, 0.5}
	r := 0.25

	idx := tree.SearchSphere(pt, r)
	ok := make([]bool, n)
	for i := range idx {
		ok[idx[i]] = true
	}

	nOk := 0
	for i := range x {
		dx, dy, dz := x[i][0] - pt[0], x[i][1] - pt[1], x[i][2] - pt[2]
		dr := math.Sqrt(dx*dx + dy*dy + dz*dz)
		inR := dr < r
		if ok[i] == inR {
			nOk++
		}
	}

	if nOk != n {
		t.Errorf("Manual search and Tree.SearchSphere only agreed on %d " + 
			"points, not %d", nOk, n)
	}
}