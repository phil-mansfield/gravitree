import ctypes
import numpy as np
import numpy.ctypeslib as ctypeslib
import os.path as path
import time

""" gravitree.py is a Python wrapper around the Go gravitree package. It
consists of two features:
1. The Tree class, which allows for tree-code calculations of accelerations and
potentials.
2. unbind(), a function which uses the Tree class to perform iterative
unbinding.

Example usage:

    # Setup

    mp, eps, G = 1e4, 1.0, graivtree.G_COSMO # Msun, kpc, km/s units
    x = # array with entry for each particle

    # Tree initialization

    my_non_default_params = gravitree.TreeParameters(leaf_size=64)
    t = gravitree.Tree(x, mp, eps, G, param=my_non_default_params)

    # Compute quantities

    pe = tree.potenial() # units of (km/s)^2
    acc = tree.acceleration() # units of (km/s)^2 / kpc
    acc *= gravitree.COSMO_ACC_TO_KPC_MYR_ACC # converts to units of kpc/Myr^2

    # Evaluate at another set of points

    x0 = # ...
    pe0 = tree.potential(x0) # evaluated at x0
"""

##################################
# C interoperability boilerplate #
##################################

# First, we load in the many C wrapper functions. These is nothing that a
# normal maintainer needs to worry about here.

# TODO: write a helper function which automates the boilerplate. Something
# like _c_potential = wrap_c_function(gravitree_lib.cPotential, "iDdDD")

file_name = path.abspath(__file__)
lib_name = path.join(path.dirname(file_name), "gravitree_wrapper.so")
gravitree_lib = ctypes.cdll.LoadLibrary(lib_name)

def _wrap_c_func(c_func, arg_string):
    """ _wrap_c_function is an internal helper function which sets up result
    and argument types. c_func is an external function from a LoadLibrary
    call. arg_string is a string with one character for each argument. Those
    characters give the types of each argument:
    i - integer
    I - integer array
    d - double
    D - double array
    You can add more as needed.

    It return the input function to make one-liners easier. (Sue me.)
    """
    # If you are copy-and-pasting this code into a new project, note that
    # sometimes you don't want restype = None. Read up on how this stuff works
    # before modifying it.
    c_func.restype = None
    args = []
    for c in arg_string:
        if c == "i":
            args.append(ctypes.c_int)
        elif c == "d":
            args.append(ctypes.c_double)
        elif c == "I":
            args.append(
                ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")
            )
        elif c == "D":
            args.append(
                ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
            )
        else:
            raise ValueError("Unrecognized type character, '%s'", c)
    c_func.argtypes = args
    return c_func
        

_c_potential = _wrap_c_func(gravitree_lib.cPotential, "iDdDD")
_c_potential_at = _wrap_c_func(gravitree_lib.cPotentialAt, "iDiDdDD")
_c_bf_potential = _wrap_c_func(gravitree_lib.cBruteForcePotential, "iDdD")
_c_bf_potential_at = _wrap_c_func(gravitree_lib.cBruteForcePotentialAt, "iDiDdD")

_c_acceleration = _wrap_c_func(gravitree_lib.cAcceleration, "iDdDD")
_c_acceleration_at = _wrap_c_func(gravitree_lib.cAccelerationAt, "iDiDdDD")
_c_bf_acceleration = _wrap_c_func(gravitree_lib.cBruteForceAcceleration, "iDdD")
_c_bf_acceleration_at = _wrap_c_func(gravitree_lib.cBruteForceAccelerationAt, "iDiDdD")

_c_potential = gravitree_lib.cPotential
_c_potential.restype = None
_c_potential.argtypes = [
    ctypes.c_int,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ctypes.c_double,
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
]

#############
# Constants #
#############

###########################
# Gravitational constants #
###########################

# units: [m^3 / (kg s^2)]
G_MKS = 6.67430e-11
# units: [kpc (km/2)^2 / Msun]
G_COSMO = 4.30092e-6
# Multiply cosmological accelerations [(km/s)^2/kpc] by this to convert
# to kpc/Gyr^2
COSMO_ACC_TO_KPC_GYR_ACC = 1.045

#########################
# Tree opening criteria #
#########################

PKDGRAV_CRITERIA = 0
SALMON_WARREN_HUT_CRITERIA = 1
BARNES_HUT_CRITERIA = 2

#############################
# Field approximation order #
#############################

MONOPOLE_ORDER = 0
QUADRUPOLE_ORDER = 1 # TODO: implement quadrupole order

class TreeParameters(object):
    def __init__(self, leaf_size=16, criteria=PKDGRAV_CRITERIA,
                 theta=0.7, order=MONOPOLE_ORDER):
        """ TreeParameters is a class containing configuration options for
        gravitree. The average user should not mess around with this unless
        they're comfortable doing manual tests on force accuracy.
        
        leaf_size - the maximum number of particles per leaf node
        criteria - the opening criteria for the tree (i.e. how far away do
        you need to be from a node before it gets split up)
        theta - all opening criteria have an additional nuissance parameter
        controlling approximaiton accuracy. For simplicity, all these
        parameters are called "theta", following the Barnes-Hut convention.
        order - the approximation order used when a node is not split up.
        """
        self.leaf_size = leaf_size
        self.criteria = criteria
        self.theta = theta
        self.order = order

    def to_array(self):
        """ to_array converts TreeParameters to a floating point array so
        it can be more easily passed the Go API. All the integer flags used
        here can be exactly represented in floating point, so it's safe to
        store them in that format.
        """
        param = [self.leaf_size, self.criteria, self.theta, self.order]
        return np.array(param, dtype=np.float64)

class Tree(object):
    def __init__(self, x, eps, mp, G, param=None):
        """ Tree is a class wrapping the Go gravitree.Tree methods.

        x - array of position vectors of particles which are generating the
        gravitational field
        eps - Plummer force-softening scale, must be in the same units as x
        mp - particle mass
        G - gravitational constant. This should be chosen to match the units
        of x, eps, and mp and will usually be one the gravitree.G_* constants
        param - a TreeParameters instance if the user wants to do advanced
        configuration shenanigans. If not supplied, the default TreeParameters
        values will be used.

        Performance note - a new tree is generated every time a method like
        potential() or acceleration() is called. This is sub-optimal, but
        I haven't figured out how to get Go to allocate the tree in a way that
        can be stored in Python's memory without errors yet.
        """
        self.x = x
        self.eps = eps
        self.mp = mp
        self.G = G
        
        if param is None:
            self.param = TreeParameters()
        else:
            self.param = param
        
    def potential(self, x=None, brute_force=False):
        """ potential computes the potential at each point

        x - if set to None, potential will be calculated at the the same
        points that are generating the gravitational field. Otherwise the
        potential will be calculated at x
        brute_force - If set to False, a tree will be used to calculate the
        potential. If set to True, the calculation will be done with an O(n^2)
        brute force calculation
        """
        n0 = len(self.x)
        x0 = np.ascontiguousarray(self.x.reshape(3*n0), dtype=np.float64)
        if x is None:
            E = np.zeros(n0, dtype=np.float64)
            if brute_force:
                _c_bf_potential(
                    n0, x0, self.eps, E)
            else:
                _c_potential(
                    n0, x0, self.eps, E, self.param.to_array())
        else:
            n1 = len(x)
            x1 = np.ascontiguousarray(x.reshape(3*n1), dtype=np.float64)
            E = np.zeros(n1, dtype=np.float64)
            if brute_force:
                _c_bf_potential_at(
                    n0, x0, n1, x1, self.eps, E)
            else:
                _c_potential_at(
                    n0, x0, n1, x1, self.eps, E, self.param.to_array())

        E *= self.mp*self.G
        return E

    def acceleration(self, x=None, brute_force=False):
        """ acceleration computes the potential at each point

        x - if set to None, acceleration will be calculated at the the same
        points that are generating the gravitational field. Otherwise the
        acceleration will be calculated at x
        brute_force - If set to False, a tree will be used to calculate the
        acceleration. If set to True, the calculation will be done with a
        O(n^2) brute force calculation

        Note - when using cosmological units (i.e., kpc, km/s, Msun, and
        G_COSMO), the accelrations are returned in units of (km/s)^2/kpc, which
        is almost certianly not the set of units you want. Multiply this by
        gravitree.
        """
        n0 = len(self.x)
        x0 = np.ascontiguousarray(self.x.reshape(3*n0), dtype=np.float64)
        if x is None:
            a = np.zeros(3*n0, dtype=np.float64)
            if brute_force:
                _c_bf_acceleration(
                    n0, x0, self.eps, a)
            else:
                _c_acceleration(
                    n0, x0, self.eps, a, self.param.to_array())
        else:
            n1 = len(x)
            x1 = np.ascontiguousarray(x.reshape(3*n1), dtype=np.float64)
            a = np.zeros(3*n1, dtype=np.float64)
            if brute_force:
                _c_bf_acceleration_at(
                    n0, x0, n1, x1, self.eps, a)
            else:
                _c_acceleration_at(
                    n0, x0, n1, x1, self.eps, a, self.param.to_array())
                
        a = a.reshape(len(a)//3, 3)
        a *= self.mp*self.G
        return a

def unbind(t, v, iters=-1, brute_force=False, return_diagnostics=False,
           method="direct"):
    """ unbind performs an unbinding operation on the particles contained
    in a tree and returns the energy of the particles.

    t - a Tree instance
    v - the velocity of each particle in the same units at 
    iters - the maximum number of iterations to use. If set to -1, iterate
    until convergence
    brute_force - if True, use an O(N^2) exact brute-force calculation
    return_diagnostic - if True, change the return value to (E, d), where E is
    the normal energy return value and d is a list of unbinding information at
    the end of every iter. Each element in d is a tuple (n_tot, n_bound, dt),
    where n_tot is the total number of particles considered that iteration,
    n_bound is the number that were bound at the end of the iteration, and 
    dt is the number of seconds the iteration took.

    Currently supported methods:
    "direct" - The standard naive unbinding where a tree is constructed from
    the particles which were bound in the previous are unbound again.
    "inverse" - A version of HBT+'s umbinding routine where a tree is made for
    the unbound particles, not the bound ones
    """

    eps, mp, G, param, x = t.eps, t.mp, t.G, t.param, t.x
    ke = np.sum(v**2, axis=1)/2
    pe = np.zeros(len(v))
    d = []
    
    if method == "direct":
        ok_prev = np.ones(len(ke), dtype=bool)
        n = 0
        while iters == -1 or n < iters:
            t0 = time.time()
            
            t = Tree(x[ok_prev], eps, mp, G, param=param)
            pe[~ok_prev] = np.inf
            pe[ok_prev] = t.potential(brute_force=brute_force)
            ok = pe + ke < 0
            
            t1 = time.time()

            if return_diagnostics:
                d.append((np.sum(ok_prev), np.sum(ok), t1 - t0))
            
            n += 1
            if np.all(ok_prev == ok): break
            ok_prev = ok
        ok_prev = ok
    elif method == "inverse":
        ok_prev = np.ones(len(ke), dtype=bool)
        n = 0
        while iters == -1 or n < iters:
            t0 = time.time()

            if n == 0 or np.sum(ok_prev)/4 < np.sum(is_changed):
                t = Tree(x[ok_prev], eps, mp, G, param=param)
                pe[~ok_prev] = np.inf
                pe[ok_prev] = t.potential(brute_force=brute_force)
            else:
                t = Tree(x[is_changed], eps, mp, G, param=param)
                pe[~ok_prev] = np.inf
                pe[ok_prev] -= t.potential(x[ok_prev],
                                           brute_force=brute_force)

            ok = pe + ke < 0            
            t1 = time.time()

            if return_diagnostics:
                d.append((np.sum(ok_prev), np.sum(ok), t1 - t0))
            
            n += 1
            if np.all(ok_prev == ok): break
            is_changed = (~ok) & ok_prev
            ok_prev = ok
        ok_prev = ok
    else:
        raise ValueError("Unrecognized binding energy method, '%s'" % method)

    pe[~ok] = np.inf
    
    if return_diagnostics:
        return pe + ke, d
    else:
        return pe + ke
    
def test():
    eps = 0.0
    mp = 2.0
    G = 1e-3
    x = np.array([
        [0,  0, 0],
        [0, +1, 0],
        [0, -1, 0]
    ])
    
    param = TreeParameters(
        leaf_size=16, criteria=PKDGRAV_CRITERIA,
        theta=0.7, order=MONOPOLE_ORDER
    )
    t = Tree(x, eps, mp, G, param=param)

    E = t.potential()
    print("Expected potential:   [-0.004 -0.003 -0.003]")
    print("Calculated potential:", E)
    x = np.array([[0, 0, 1], [0, 0, 2]])
    # Just check these can run without crashing for now
    E = t.potential(x)
    E = t.potential(brute_force=True)
    E = t.potential(x, brute_force=True)

    a = t.acceleration()
    print("""Expected acceleration:
[[ 0.      0.      0.    ]
 [ 0.     -0.0025  0.    ]
 [ 0.      0.0025  0.    ]]""")
    print("Calculated acceleration:\n", a)
    x = np.array([[0, 0, 1], [0, 0, 2]])
    # Just check these can run without crashing for now
    a = t.acceleration(x)
    a = t.acceleration(brute_force=True)
    a = t.acceleration(x, brute_force=True)
    
if __name__ == "__main__":
    test()
