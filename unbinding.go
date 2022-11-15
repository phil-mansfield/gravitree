package gravitree

import (
	"math"
)

func BindingEnergy(x, v [][3]float64, mp, eps float64, E []float64) {
	IterativeBindingEnergy(x, v, mp, eps, 1, E)
}

func IterativeBindingEnergy(
	x, v [][3]float64, mp, eps float64, iters int, E []float64,
) {
	ok := make([]bool, len(x))
	nPrev := len(x)
	inf := math.Inf(+1)
	for i := range ok { ok[i] = true }

	for j := 0; j < iters || nPrev == 0; j++ {
		nOk := 0
		for i := range ok {
			if ok[i] { nOk++ }
		}

		bindingEnergy(x, v, ok, mp, eps, E)

		nCurr := 0
		for i := range ok {
			ok[i] = ok[i] && E[i] < 0
			if ok[i] {
				nCurr++
			} else {
				E[i] = inf
			}
		}
		
		if nCurr == nPrev { break }
		nPrev = nCurr
	}
}

func bindingEnergy(
	x, v [][3]float64, ok []bool, mp, eps float64, E []float64,
) {
	n0 := 0
	for i := range ok {
		if ok[i] { n0++ }
	}
	
	x0 := make([][3]float64, n0)
	v0 := make([][3]float64, n0)
	E0 := make([]float64, n0)

	j := 0
	for i := range x {
		if ok[i] {
			x0[j], v0[j] = x[i], v[i]
			j++
		}
	}

	tree := NewTree(x0)
	tree.Potential(eps, E0)

	j = 0
	nOk := 0
	for i := range ok {
		if ok[i] { nOk++ }
	}

	for i := range ok {
		if ok[i] {
			v2 := 0.0
			for dim := 0; dim < 3; dim++ {
				v2 += v0[j][dim]*v0[j][dim]
			}
			E[i] = E0[j]*mp*4.301e-6 + v2/2
			j++
		}
	}
}

/*
func RadialBindingEnergy(x, v [][3]float64, mp, eps float64, E []float64) {
	ok := make([]bool, len(x))
	for i := range ok { ok[i] = true }
	radialBindingEnergy(x, v, ok, mp, eps, E)
}


func IterativeRadialBindingEnergy(
	x, v [][3]float64, mp, eps float64, iters int, E []float64,
) {
	ok := make([]bool, len(x))
	nPrev := len(x)
	for i := range ok { ok[i] = true }

	for j := 0; j < iters || nPrev == 0; j++ {
		radialBindingEnergy(x, v, ok, mp, eps, E)

		nCurr := 0
		for i := range ok {
			ok[i] = ok[i] && E[i] < 0
			if ok[i] { nCurr++ }
		}

		if nCurr == nPrev { break }
		nPrev = nCurr
	}

	inf := math.Inf(+1)
	for i := range ok {
		if !ok[i] { E[i] = inf }
	}
}


func radialBindingEnergy(
	x, v [][3]float64, ok []bool, mp, eps float64, E []float64,
) {
	n0 := 0
	for i := range ok {
		if ok[i] { n0++ }
	}
	
	x0 := make([][3]float64, len(x))
	v0 := make([][3]float64, len(v))

	j := 0
	for i := range x {
		if ok[i] {
			x0[j], v0[j] = x[i], v[i]
			j++
		}
	}

	rMin2 := eps*eps
	rMax2 := eps*eps*1.01
	bins := 500
	
	r20 := make([]float64, len(x0))
	KE0 := make([]float64, len(x0))
	for i := range r {
		for dim := 0; dim < 3; dim++ {
			r20[i] += x0[i]*x0[i]
			KE0[i] += v0[i]*v0[i]/2
		}
		if rMax2 < r20[i] { rMax2 = r20[i] }
	}

	r, mEnc := enclosedMass(rMin2, rMax2*1.01, mp, bins, r20)
	E0 := massProfileToPhi(r, mEnc)

	j = 0
	for i := range ok {
		if ok[i] {
			E[i] += E0[j] + KE0[j]
			j++
		}
	}
}

func enclosedMass(rMin2, rMax2, mp float64, bins int, r2 []float64) {
	lrMin, lrMax := math.Log10(rMin2)/2, math.Log10(rMax2)
	dlr := (lrMax - lrMin) / float64(bins - 1)

	// Reread this whole function and trace on paper: off by one errors.
	
	m := make([]float64, bins+1)
	for i := range r2 {
		lr := math.Log10(r2)/2
		i := int(lr-lrMin)
		
		if i < 0 {
			i = -1
		} else if i >= bins {
			continue
		}

		m[i+1] += mp
	}

	for i := 1; i < len(m); i++ {
		m[i] += m[i-1]
	}

	r := make([]float64, bins+1)
	// Uhhh, I'm not sure

	return r, m
}

func massProfileToPhi(r, m []float64) []float64 {
	// Tranlate:
	
	//dW = np.zeros(len(m_enc))

    //dr = r_sort[1:] - r_sort[:-1]
    //dr[dr == 0] = np.finfo(dr.dtype).eps
    //r_scaled = (r_sort[1:]*r_sort[:-1]) / dr
    //dW[:-1] = v_circ(m_enc[:-1], r_scaled)**2
                     
    //vesc_lim = v_circ(m_enc[-1], r_sort[-1]) * np.sqrt(2)

    //W = (np.cumsum(dW[::-1])[::-1] + 2*vesc_lim**2)/vmax**2

}
*/
