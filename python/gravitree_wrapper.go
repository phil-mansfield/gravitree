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

/*
func ptrToFloat64(ptr unsafe.Pointer, n int) []float64 {
	slice := []float64{ }
	hd := (*reflect.SliceHeader)(unsafe.Pointer(&slice))
	hd.Data, hd.Len, hd.Cap = uintptr(ptr), n, n
	return slice
}

func ptrToVec64(ptr unsafe.Pointer, n int) [][3]float64 {
	slice := [][3]float64{ }
	hd := (*reflect.SliceHeader)(unsafe.Pointer(&slice))
	hd.Data, hd.Len, hd.Cap = uintptr(ptr), n, n
	return slice
}
*/

func main() {
	
}
