package utils

import "github.com/phil-mansfield/symtable"

func ReadPointFile(filename string) [][3]float64 {

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

func ReadPhaseSpaceFile(filename string) [][6]float64 {

	t := symtable.TextFile(filename)

	cols := t.ReadFloat64s([]int{0, 1, 2, 3, 4, 5}) // column indices
	xs := cols[0]
	ys := cols[1]
	zs := cols[2]
	vxs := cols[3]
	vys := cols[4]
	vzs := cols[5]

	var result [][6]float64

	for i := 0; i < len(xs); i++ {
		var point [6]float64
		point[0] = xs[i]
		point[1] = ys[i]
		point[2] = zs[i]
		point[3] = vxs[i]
		point[4] = vys[i]
		point[5] = vzs[i]
		result = append(result, point)
	}

	return result
}
