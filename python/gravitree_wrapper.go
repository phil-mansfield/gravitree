package main

import (
	"github.com/phil-mansfield/gravitree"

	"C"
	"unsafe"
)

/*
//  //export cIterativeBindingEnergy
func cIterativeBindingEnergy(
	np C.int, x *C.double, v *C.double,
	mp C.double, eps C.double, nIter C.int, E *C.double,
) {
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	gv := unsafe.Slice((*[3]float64)(unsafe.Pointer(v)), int(np))	
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(np))	

	gravitree.IterativeBindingEnergy(
		gx, gv, float64(mp), float64(eps), int(nIter), gE)
}

//export cSplitAceleration
func cSplitAcceleration(
	np C.longlong, x *C.double, mp, eps C.double,
	nActive, nTest C.longlong,
	a *C.double,
) {
	np := nActive + nTest
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), int(np))

	t1 := gravitree.NewTree(gx[:nActive])
	t2 := gravitree.NewArrayTree(g[nActive:])
	acc := gravitree.Acceleration(ga[:nActive])
	t1.EvaluateAt(t2, eps, acc)

	for i := range {
	}
}
*/

func paramToOptions(param *C.double) gravitree.TreeOptions {
	p := unsafe.Slice((*float64)(unsafe.Pointer(param)), 4)
	
	return gravitree.TreeOptions{
		LeafSize: int(p[0]),
		Criteria: gravitree.OpeningCriteria(p[1]),
		Theta: p[2],
		Order: gravitree.ApproximationOrder(p[3]),
	}
}

//export cPotential
func cPotential(
	np C.longlong, x *C.double,
	eps C.double, E *C.double,
	param *C.double,
) {
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(np))

	opt := paramToOptions(param)
	tree := gravitree.NewTree(gx, opt)
	tree.Evaluate(float64(eps), gravitree.Potential(gE))
}

//export cPotentialAt
func cPotentialAt(
	n0 C.longlong, x0 *C.double,
	n1 C.longlong, x1 *C.double,
	eps C.double, E *C.double,
	param *C.double,
) {
	gx0 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x0)), int(n0))
	gx1 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x1)), int(n1))
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(n1))

	opt := paramToOptions(param)
	t0 := gravitree.NewTree(gx0, opt)
	t1 := &gravitree.NewArrayTree(gx1).Tree
	t0.EvaluateAt(t1, float64(eps), gravitree.Potential(gE))
}

//export cBruteForcePotential
func cBruteForcePotential(
	n0 C.longlong, x0 *C.double,
	eps C.double, E *C.double,
	param *C.double,
) {
	gx0 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x0)), int(n0))
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(n0))

	gravitree.BruteForcePotential(float64(eps), gx0, gE)
}

//export cBruteForcePotentialAt
func cBruteForcePotentialAt(
	n0 C.longlong, x0 *C.double,
	n1 C.longlong, x1 *C.double,
	eps C.double, E *C.double,
) {
	gx0 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x0)), int(n0))
	gx1 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x1)), int(n1))
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(n1))

	gravitree.BruteForcePotentialAt(float64(eps), gx0, gx1, gE)
}

//export cAcceleration
func cAcceleration(
	np C.longlong, x *C.double,
	eps C.double, a *C.double,
	param *C.double,
) {
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), int(np))

	opt := paramToOptions(param)
	tree := gravitree.NewTree(gx, opt)
	tree.Evaluate(float64(eps), gravitree.Acceleration(ga))
}

//export cAccelerationAt
func cAccelerationAt(
	n0 C.longlong, x0 *C.double,
	n1 C.longlong, x1 *C.double,
	eps C.double, a *C.double,
	param *C.double,
) {
	gx0 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x0)), int(n0))
	gx1 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x1)), int(n1))
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), int(n1))

	opt := paramToOptions(param)
	t0 := gravitree.NewTree(gx0, opt)
	t1 := &gravitree.NewArrayTree(gx1).Tree
	t0.EvaluateAt(t1, float64(eps), gravitree.Acceleration(ga))
}

//export cBruteForceAcceleration
func cBruteForceAcceleration(
	n0 C.longlong, x0 *C.double,
	eps C.double, a *C.double,
	param *C.double,
) {
	gx0 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x0)), int(n0))
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), int(n0))

	gravitree.BruteForceAcceleration(float64(eps), gx0, ga)
}

//export cBruteForceAccelerationAt
func cBruteForceAccelerationAt(
	n0 C.longlong, x0 *C.double,
	n1 C.longlong, x1 *C.double,
	eps C.double, a *C.double,
) {
	gx0 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x0)), int(n0))
	gx1 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x1)), int(n1))
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), int(n1))

	gravitree.BruteForceAccelerationAt(float64(eps), gx0, gx1, ga)
}

func main() { }
