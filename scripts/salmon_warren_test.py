import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc

def r_max(x, x0):
    dx = np.zeros(x.shape)
    for k in range(3): dx[:,k] = x[:,k] - x0[k]
    r2 = np.sum(dx**2, axis=1)
    r2_max = np.max(r2)
    return np.sqrt(r2_max)

def center_of_mass(x, m=None):
    if m is None:
        return np.sum(x, axis=0) / x.shape[0]
    else:
        xm = np.zeros(x.shape)
        for k in range(3): xm[:,k] = x[:,k]*m
        return np.sum(xm, axis=0) / np.sum(m)

def sigma_x2(x, x0, m=None):
    dx = np.zeros(x.shape)
    for k in range(3): dx[:,k] = x[:,k] - x0[k]
    norm2 = np.sum(dx**2, axis=1)
    if m is None:
        return np.sum(norm2) / x.shape[0]
    else:
        return np.sum(norm2*m) / np.sum(m)

def salmon_warren_limit(x, d, m=None):
    x0 = center_of_mass(x, m)
    rm = r_max(x, x0)
    sx2 = sigma_x2(x, x0, m)
    return 1/(1 - rm/d)**2 * (sx2 / d**2)
    
def direct_potential(x, m=None):
    """ This is done in units where G = 1.
    """
    if m is None: m = np.ones(x.shape[0])
    
    N = x.shape[0]
    phi = np.zeros(N)
    for i in range(N):
        dx = np.zeros(x.shape)
        for k in range(3): dx[:,k] = x[:,k] - x[i,k]
        dx = dx[np.arange(N) != i,:]
        inv_r = np.sum(dx**2, axis=1)**-0.5
        phi[i] = -np.sum(inv_r)

    return phi

def direct_potential_split(y, x):
    N = x.shape[0]
    phi = np.zeros(N)
    for i in range(N):
        dy = np.zeros(y.shape)
        for k in range(3): dy[:,k] = x[i,k] - y[:,k]
        inv_r = np.sum(dy**2, axis=1)**-0.5
        phi[i] = -np.sum(inv_r)

    return phi

def test_point_line(unit, d):
    out = np.zeros((len(d), 3))
    for k in range(3):
        out[:,k] = unit[k]*d
    return out

def main():
    palette.configure(False)
    
    x = np.array([[0, 0, 0], [0, 0, 1]], dtype=float)
    x0 = center_of_mass(x)
    for k in range(3): x[:,k] -= x0[k]
    
    d = np.linspace(2, 10, 100)
    mono_phi = -len(x) / d
    lim = salmon_warren_limit(x, d) * mono_phi

    n = len(x)
    
    units = np.array([[0, 0, 1], [0, 1, 0], [0, 1/np.sqrt(2), 1/np.sqrt(2)]])
    colors = [pc("r"), pc("o"), pc("b")]
    
    for i in range(len(units)):
        point_line = test_point_line(units[i], d)
        phi = direct_potential_split(x, point_line)
        
        plt.plot(d, phi, c=colors[i])
        

    plt.fill_between(d, mono_phi+lim, mono_phi-lim,
                     alpha=0.2, color="k")
    plt.plot(d, mono_phi, c="k")
    
    plt.show()

    
if __name__ == "__main__":
    main()
