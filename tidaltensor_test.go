package gravitree

import (
	"testing"

	"github.com/phil-mansfield/symtable"
	"gonum.org/v1/gonum/mat"
)

func TestTidalTensorSinglePoint(t *testing.T) {

	tests := []struct {
		tT [][3][3]float64
	}{
		{[][3][3]float64{
			{{1.125, 0, 0},
				{0, 1.125, 0},
				{0, 0, -2.25}},
			{{2, 0, 0},
				{0, 2, 0},
				{0, 0, -4}},
			{{1.125, 0, 0},
				{0, 1.125, 0},
				{0, 0, -2.25}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
		}}, // pairwise
		{[][3][3]float64{
			{{1.44299611, 0, 0},
				{0, -2.14099961, -0.63766135},
				{0, -0.63766135, 0.6980035}},
			{{1.70710678, 0, 0},
				{0, -2.35355339, 0},
				{0, 0, 0.64644661}},
			{{1.44299611, 0, 0},
				{0, -2.14099961, 0.63766135},
				{0, 0.63766135, 0.6980035}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
		}}, // one-sided
		{[][3][3]float64{
			{{1.06066017, 0, 0},
				{0, -0.53033009, -1.59099026},
				{0, -1.59099026, -0.53033009}},
			{{3, 0, 0},
				{0, -6, 0},
				{0, 0, 3}},
			{{1.06066017, 0, 0},
				{0, -0.53033009, 1.59099026},
				{0, 1.59099026, -0.53033009}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}},
			{{0, 0, 0},
				{0, 0, 0},
				{0, 0, 0}}}},
	} // monopole
	tree := &Tree{}
	tree.Points = [][3]float64{{0, 0, 0}, {0, 0, 2}, {0, 0, 1},
		{0, 1, 1}, {0, 1, 0}, {0, 1, 2}}
	tree.Index = []int{0, 2, 1, 4, 3, 5}

	tT := make(TidalTensor, 6)

	tree.Nodes = []Node{{Start: 0, End: 3, Center: [3]float64{0, 0, 1}},
		{Start: 3, End: 6, Center: [3]float64{0, 1, 1}}}
	tree.eps2 = 0.0

	for i := range tests {

		test := tests[i]

		for j := range tT {
			tT[j] = [3][3]float64{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
		}

		switch i {
		case 0:
			tree.Pairwise(0, tT)
		case 1:
			tree.OneSided(0, 1, tT)
		case 2:
			tree.Monopole(0, 1, tT)
		}

		flat_mat := flatten9(tT)
		flat_test := flatten9(test.tT)

		if !multArrayAlmostEq(flat_mat, 1.0, flat_test, 1e-3) {
			t.Errorf("%d.0) expected 1*%.4f, but got %.4f.",
				i, test.tT, tT)
		}
	}
}

// should probably make a generic function so
// we don't have to keep flattening these damn things!
func flatten9(x [][3][3]float64) []float64 {
	flat := make([]float64, len(x)*9)

	for i := range x {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				flat[9*i+j+k] = x[i][j][k]
			}
		}
	}

	return flat
}

func TestTidalTensorPlummer(t *testing.T) {
	t.Skip()
	filename := "plummer.txt"
	x := readFile(filename)

	point_indices := []int{100, 200, 300, 400}

	tests := []float64{
		-0.06791001116763588, -0.08731181478956411,
		-0.049879254759283294, -0.034753012028661934}

	tree := NewTree(x)
	tensor := make(TidalTensor, len(x))
	tree.Quantity(0.0, tensor)

	// Pick a point.
	// Get tensor at that point.
	for k, point_index := range point_indices {

		tensor_at_point := tensor[point_index]

		// calculate eigenvalues of tensor
		// flatten the tensor into a []float64
		// because Eigen takes a 1x9
		// and ravels it.
		flat_mat := make([]float64, 9)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				flat_mat[i*3+j] = tensor_at_point[i][j]
			}
		}

		// this might be slow...
		// benchmark this.
		var eig mat.Eigen
		eig.Factorize(mat.NewDense(3, 3, flat_mat), mat.EigenBoth)
		eigenvalues := eig.Values(nil)

		// fmt.Printf("%v \n", eigenvalues)

		real_part := []float64{0, 0, 0}

		for i := range 3 {
			real_part[i] = float64(real(eigenvalues[i]))
		}

		min_eigenvalue := findMin(real_part)

		// check if it matches
		if !almostEq(min_eigenvalue, tests[k], 1e-3) {
			t.Errorf("(%d.0) min. eigval. expected %.4f, got %.4f",
				k, tests[k], min_eigenvalue)
		}
	}
}

// replace this with symphony_pipeline reader
// at some point.
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

func findMin(arr []float64) float64 {
	min := arr[0]

	for _, val := range arr {
		if val < min {
			min = val
		}
	}

	return min
}
