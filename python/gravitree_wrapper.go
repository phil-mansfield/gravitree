package main

//#include <stdint.h>
import "C"
import (
	"github.com/phil-mansfield/gravitree"
	"runtime/cgo"
	"unsafe"
)

func paramToOptions(param *C.double) gravitree.TreeOptions {
	// TODO: turn into a json file: easier to debug and add to
	p := unsafe.Slice((*float64)(unsafe.Pointer(param)), 4)
	
	return gravitree.TreeOptions{
		LeafSize: int(p[0]),
		Criteria: gravitree.OpeningCriteria(p[1]),
		Theta: p[2],
		Order: gravitree.ApproximationOrder(p[3]),
	}
}

//export cNewTree
func cNewTree(
	np C.longlong, x *C.double,
	param *C.double,
) C.uintptr_t {
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	opt := paramToOptions(param)
	tree := gravitree.NewTree(gx, opt)

	return C.uintptr_t(cgo.NewHandle(tree))
}

//export cFreeTree
func cFreeTree(
	ptr C.uintptr_t,
) {
	h := cgo.Handle(ptr)
	h.Delete()
}

//export cPotential
func cPotential(
	ptr C.uintptr_t,
	eps C.double, E *C.double,
) {
	h := cgo.Handle(ptr)
	tree := h.Value().(*gravitree.Tree)
	
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), len(tree.Points))
	tree.Evaluate(float64(eps), gravitree.Potential(gE))
}

//export cPotentialAt
func cPotentialAt(
	ptr C.uintptr_t,
	np C.longlong, x *C.double,
	eps C.double, E *C.double,
) {
	h := cgo.Handle(ptr)
	t0 := h.Value().(*gravitree.Tree)
	
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(np))

	// TODO: figure out a new way to handle this
	t1 := &gravitree.NewArrayTree(gx).Tree
	
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
	ptr C.uintptr_t,
	eps C.double, a *C.double,
) {
	h := cgo.Handle(ptr)
	tree := h.Value().(*gravitree.Tree)
	
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), len(tree.Points))

	tree.Evaluate(float64(eps), gravitree.Acceleration(ga))
}

//export cAccelerationAt
func cAccelerationAt(
	ptr C.uintptr_t,
	n1 C.longlong, x1 *C.double,
	eps C.double, a *C.double,
) {
	h := cgo.Handle(ptr)
	t0 := h.Value().(*gravitree.Tree)
	
	gx1 := unsafe.Slice((*[3]float64)(unsafe.Pointer(x1)), int(n1))
	ga := unsafe.Slice((*[3]float64)(unsafe.Pointer(a)), int(n1))

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
