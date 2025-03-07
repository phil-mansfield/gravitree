import numpy as np
from scipy.special import gammaincc, gamma, gammainc
from scipy.optimize import fsolve


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
        h = self.Rs / (2.0 / self.alpha) ** (1.0 / self.alpha)
        x = (r / h) ** self.alpha
        t1 = gammainc(3.0 / self.alpha, x) / x ** (1.0 / self.alpha)
        t2 = (
            gamma(2.0 / self.alpha)
            * gammaincc(2.0 / self.alpha, x)
            / gamma(3.0 / self.alpha)
        )
        psi = -(1.0 / h) * (t1 + t2)
        return psi

    def u_eff(self, r, l=0.0, offset=0.0):
        # l is specific angular momentum
        return (0.5 * (l / r) ** 2 + self.einasto_potential(r)) - offset

    def estimate_apsides(self, pos, vel):
        r0 = np.sqrt(np.sum(pos**2))

        t = np.sum(vel**2) / 2.0
        u = self.einasto_potential(r0)
        e = t + u

        l = np.linalg.norm(np.cross(pos, vel))

        r_apo = fsolve(self.u_eff, r0 * 1e-3, args=(l, e), maxfev=2000)[0]
        r_peri = fsolve(self.u_eff, r0, args=(l, e), maxfev=2000)[0]

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

    def dr(self, r):
        t1 = 4 * np.pi * r**2
        t2 = (4 * np.pi * self.Rs**3) / (3 * self.enclosed_mass(self.Rs))
        t3 = np.exp((-2 / self.alpha) * ((r / self.Rs) ** self.alpha - 1))
        return t1 * t2 * t3

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
        return (0.5 * (l / r) ** 2 + self.plummer_potential(r)) - offset

    def find_circular_velocity(self, r):
        l_circular = np.sqrt(r * self.plummer_potential(r))
        return l_circular / r

    def estimate_apsides(self, pos, vel):
        r0 = np.sqrt(np.sum(pos**2))

        t = np.sum(vel**2) / 2.0
        u = self.plummer_potential(r0)
        e = t + u

        l = np.linalg.norm(np.cross(pos, vel))

        r_peri = fsolve(self.u_eff, r0 * 1e-3, args=(l, e), maxfev=2000)[0]
        r_apo = fsolve(self.u_eff, r0, args=(l, e), maxfev=2000)[0]

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


#
class PointMass:
    def __init__(self, M):
        self.M = M

    def point_mass_potential(self, r):
        return -self.M / r

    def u_eff(self, r, l=0.0, offset=0.0):
        # l is specific angular momentum
        return (0.5 * (l / r) ** 2 + self.point_mass_potential(r)) - offset

    def find_circular_velocity(self, r):
        l_circular = np.sqrt(r * self.point_mass_potential(r))
        return l_circular / r

    def estimate_apsides(self, pos, vel):
        r0 = np.sqrt(np.sum(pos**2))

        t = np.sum(vel**2) / 2.0
        u = self.point_mass_potential(r0)
        e = t + u

        l = np.linalg.norm(np.cross(pos, vel))

        r_apo = fsolve(self.u_eff, r0 * 1e-3, args=(l, e), maxfev=2000)[0]
        r_peri = fsolve(self.u_eff, r0, args=(l, e), maxfev=2000)[0]

        return r_apo, r_peri

    def enclosed_mass(self, r):
        # For a point mass, the enclosed mass is constant and equal to M
        return self.M

    def get_circular_velocity(self, r):
        res = np.sqrt(self.enclosed_mass(r) / r)
        return res
