import ctypes
import numpy as np
import numpy.ctypeslib as ctypeslib

gravitree_lib = ctypes.cdll.LoadLibrary("./gravitree.so")
_c_iterative_binding_energy = gravitree_lib.cIterativeBindingEnergy
_c_iterative_binding_energy.restype = None
_c_iterative_binding_energy.argtypes = [
    ctypes.c_longlong,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTINUOUS"),
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTINUOUS"),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_longlong,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTINUOUS"),
]
_c_potential_energy = gravitree_lib.cPotentialEnergy
_c_potential_energy.restype = None
_c_potential_energy.argtypes = [
    ctypes.c_longlong,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTINUOUS"),
    ctypes.c_double,
    ctypes.c_double,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTINUOUS"),
]

def binding_energy(x, v, mp, eps, n_iter=1):
    E = np.zeros(len(x), dtype=np.float64)
    _c_iterative_binding_energy(len(x), x, v, mp, eps, n_iter, E)

def potential_energy(x, mp, eps):
    E = np.zeros(len(x), dtype=np.float64)
    _c_potential_energy(len(x), x, mp, eps, E)
    
def test():
    x = np.array([
        [0, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
    ])
    v = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 1e6, 0]
    ])
    eps = 1e-3
    mp = 1e8

    E1 = iterative_binding_energy(x, v, mp, eps)
    E2 = iterative_binding_energy(x, v, mp, eps, 2)

    P = potential_energy(x, mp, eps)
    
    print(E1)
    print(E2)
    print(P)
    
if __name__ == "__main__":
    test()
