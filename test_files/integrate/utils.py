import numpy as np
from scipy.special import gammaincc, gamma
from scipy.optimize import fsolve


def load_particles(filename):
    p = np.loadtxt(filename)
    if p.ndim == 1:
        p = p.reshape(1, 6)
    return p


def get_trajectory(fn, idx, snaps, full=False):
    traj = np.zeros((len(snaps), 3))

    if full:
        traj = np.zeros((len(snaps), 6))

    for k, snap in enumerate(snaps):
        file = fn % snap
        pos = load_particles(file)
        if not full:
            traj[k] = pos[idx, :3]
        else:
            traj[k] = pos[idx]
    return traj


def get_quantiles(data, axis=0):
    p50 = np.quantile(data, 0.50, axis)
    p84 = np.quantile(data, 0.84, axis)
    p16 = np.quantile(data, 0.16, axis)
    return p16, p50, p84


# source: https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
def win_avg(x, n):
    return np.convolve(x, np.ones(n) / n, mode="valid")


###### einasto ######


class Einasto:
    def __init__(self, Rs, alpha):
        self.Rs = Rs
        self.alpha = alpha

    def get_d_constant(self):
        n = 1 / self.alpha
        d = 3.0 * n
        d += -(1.0 / 3.0)
        d += 8.0 / (1215.0 * n)
        d += 184.0 / (229635.0 * n * n)
        d += 1084.0 / (31000725 * n * n * n)
        d += -(17557576.0 / (1242974068875.0 * n * n * n * n))
        return d

    def einasto_potential(self, r):
        s = (r / self.Rs) * (self.get_d_constant() ** (1.0 / self.alpha))
        s_alpha = s**self.alpha

        t1 = gammaincc(3.0 / self.alpha, s_alpha)
        t2 = (
            gammaincc(2.0 / self.alpha, s_alpha)
            * gamma(2.0 / self.alpha)
            / gamma(3.0 / self.alpha)
        )
        psi = (1 / r) * (1 - t1 + t2)
        return psi

    def u_eff(self, r, l=0.0, offset=0.0):
        # l is specific angular momentum
        return (0.5 * (l / r) ** 2 - self.einasto_potential(r)) - offset

    def find_circular_velocity(self, r):
        l_circular = np.sqrt(r * self.einasto_potential(r))
        return l_circular / r

    def estimate_apsides(self, pos, vel):
        r0 = np.sqrt(np.sum(pos**2))

        t = np.sum(vel**2) / 2.0
        u = self.einasto_potential(r0)
        e = t - u

        l = np.linalg.norm(np.cross(pos, vel))

        r_apo = fsolve(self.u_eff, r0 * 0.1, args=(l, e))[0]
        r_peri = fsolve(self.u_eff, r0 * 1.5, args=(l, e))[0]

        return r_apo, r_peri

    def enclosed_mass(self, r):
        # Mtot is the total mass
        # since code units are [m] = 1/npts, then
        # Mtot = 1 when all points are enclosed.
        # return lower incomplete gamma function

        a = 3.0 / self.alpha
        x = (2.0 / self.alpha) * (r / self.Rs) ** self.alpha
        res = 1 - gammaincc(a, x)
        return res

    def get_circular_velocity(self, r):
        res = np.sqrt(self.enclosed_mass(r) / r)
        return res
        ###### plummer ######


class Plummer:
    def __init__(self, b):
        self.b = b

    def plummer_potential(self, r):
        return -1 / np.sqrt(r**2 + self.b**2)

    def u_eff(self, r, l=0.0, offset=0.0):
        # l is specific angular momentum
        return (0.5 * (l / r) ** 2 - self.plummer_potential(r)) - offset

    def find_circular_velocity(self, r):
        l_circular = np.sqrt(r * self.plummer_potential(r))
        return l_circular / r

    def estimate_apsides(self, pos, vel):
        r0 = np.sqrt(np.sum(pos**2))

        t = np.sum(vel**2) / 2.0
        u = self.plummer_potential(r0)
        e = t - u

        l = np.linalg.norm(np.cross(pos, vel))

        r_apo = fsolve(self.u_eff, r0 * 0.1, args=(l, e))[0]
        r_peri = fsolve(self.u_eff, r0 * 1.5, args=(l, e))[0]

        return r_apo, r_peri

    def enclosed_mass(self, r):
        # Mtot is the total mass
        # since code units are [m] = 1/npts, then
        # Mtot = 1 when all points are enclosed.
        # return mass enclosed within radius r

        return r**3 / (r**2 + self.b**2) ** (3 / 2)

    def get_circular_velocity(self, r):
        res = np.sqrt(self.enclosed_mass(r) / r)
        return res
