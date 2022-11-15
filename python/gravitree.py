import ctypes
import numpy as np
import numpy.ctypeslib as ctypeslib

gravitree_lib = ctypes.cdll.LoadLibrary("/home/users/phil1/code/src/github.com/phil-mansfield/gravitree/python/gravitree_wrapper.so")
_c_iterative_binding_energy = gravitree_lib.cIterativeBindingEnergy
_c_iterative_binding_energy.restype = None
_c_iterative_binding_energy.argtypes = [
    ctypes.c_int,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
]

def binding_energy(x, v, mp, eps, n_iter=1):
    n = len(x)
    E = np.zeros(n, dtype=np.float64)
    x = np.ascontiguousarray(x.reshape((3*n,)), dtype=np.float64)
    v = np.ascontiguousarray(v.reshape((3*n,)), dtype=np.float64)
    _c_iterative_binding_energy(n, x, v, mp, eps, n_iter, E)
    return E

def test():
    x = np.array([
        [0, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
    ])
    v = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 100, 0]
    ])
    eps = 1e-3
    mp = 1e8

    print("test time")
    E1 = binding_energy(x, v, mp, eps)
    E2 = binding_energy(x, v, mp, eps, 2)

    print(E1)
    print(E2)
    
if __name__ == "__main__":
    test()
