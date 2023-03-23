package main

import (
	"github.com/phil-mansfield/gravitree"

	"C"
	"unsafe"
)

//export cIterativeBindingEnergy
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

//export cPotentialEnergy
func cPotentialEnergy(
	np C.longlong, x *C.double,
	mp, eps C.double, E *C.double,
) {
	gx := unsafe.Slice((*[3]float64)(unsafe.Pointer(x)), int(np))
	gE := unsafe.Slice((*float64)(unsafe.Pointer(E)), int(np))

	gravitree.PotentialEnergy(gx, float64(mp), float64(eps), gE)
}

func main() {
	
}
