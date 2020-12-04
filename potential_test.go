package gravitree

import (
	"testing"
)

func TestPotentialPrimitives(t *testing.T) {
	tests := []struct{
		phi []float64
	}{
		{[]float64{-1.5, -2, -1.5, 0, 0, 0}}, // pairwise
		{[]float64{-2.15432, -2.41421, -2.15432, 0, 0, 0}}, // one-sided
		{[]float64{-2.12132, -3, -2.12132, 0, 0, 0}}, // monopole
	}

	tree := &Tree{ }
	tree.Points = [][3]float64{ {0, 0, 0}, {0, 0, 2}, {0, 0, 1},
		{0, 1, 1}, {0, 1, 0}, {0, 1, 2}}
	tree.Index = []int{ 0, 2, 1, 4, 3, 5 }
	phi := make([]float64, 6)
	
	tree.Nodes = []Node{ {Start: 0, End: 3, Center: [3]float64{0, 0, 1}},
		{Start: 3, End: 6, Center: [3]float64{0, 1, 1}} }
	tree.eps2 = 0.0

	for i := range tests {
		test := tests[i]
		
		for j := range phi { phi[j] = 0 }

		switch i {
		case 0: tree.pairwisePotential(0, phi)
		case 1: tree.oneSidedPotential(0, 1, phi)
		case 2: tree.monopolePotential(0, 1, phi)
		}

		if !multArrayAlmostEq(phi, 1.0, test.phi, 1e-3) {
			t.Errorf("%d.0) expected phi = 1*%.4f, but phi = %.4f.",
				i, test.phi, phi)
		}
		
		switch i {
		case 0: tree.pairwisePotential(0, phi)
		case 1: tree.oneSidedPotential(0, 1, phi)
		case 2: tree.monopolePotential(0, 1, phi)
		}

		if !multArrayAlmostEq(phi, 2.0, test.phi, 1e-3) {
			t.Errorf("%d.1) expected phi = 2*%.4f, but phi = %.4f.",
				i, test.phi, phi)
		}
	}
}

func multArrayAlmostEq(
	x []float64, mult float64, y []float64, eps float64,
) bool {
	if len(x) != len(y) { return false }
	
	for i := range x {
		if !almostEq(x[i], mult*y[i], eps) { return false }
	}
	return true
}
