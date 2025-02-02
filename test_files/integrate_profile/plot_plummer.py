import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from profiles import Plummer
from scipy.optimize import fmin

########

vc = 0.5

########


mpl.rcParams["font.family"] = "monospace"


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


profile = Plummer(0.05)  # Adjust parameters as needed

###

# initial conditions
x0 = np.array([1.0, 0.0, 0.0])
v0 = np.array([0.0, vc * profile.get_circular_velocity(np.linalg.norm(x0)), 0.0])
e0 = ((np.linalg.norm(v0) ** 2.0) / 2.0) + profile.plummer_potential(np.linalg.norm(x0))

print(x0, v0)
l0 = np.linalg.norm(np.cross(x0, v0))
apo, peri = profile.estimate_apsides(x0, v0)

###

integss = ["leapfrog", "rk4"]
dtss = ["1e-2", "1e-3", "1e-4"]

tfn = "time.dat"
efn = "energy.dat"
afn = "acc.dat"

fig, ax = plt.subplots(1, 4, figsize=(15, 5))
colors = {"leapfrog": "cadetblue", "rk4": "firebrick"}
ls = ["-", "-.", "--"]

r_peri = []
r_apo = []

for i, integs in enumerate(integss):
    for j, dts in enumerate(dtss):
        directory = f"snapshots/plummer_int={integs}_dt={dts}"

        e = np.loadtxt(os.path.join(directory, efn))
        a = np.loadtxt(os.path.join(directory, afn))
        t = np.loadtxt(os.path.join(directory, tfn))

        e_err = np.abs(1 - e / e[0])
        ax[0].plot(t, e_err, color=colors[integs], ls=ls[j], alpha=0.2)
        traj_k = get_trajectory(
            os.path.join(directory, "t_%s.dat"), 0, np.arange(0, len(e), 1)
        )

        r_k = np.sqrt(np.sum(traj_k**2, axis=1))

        ax[1].plot(t, r_k, color=colors[integs], ls=ls[j], alpha=0.2)

        ax[2].plot(
            traj_k[:, 0], traj_k[:, 1], color=colors[integs], ls=ls[j], alpha=0.2
        )

        r_apo.append(np.min(r_k))
        r_peri.append(np.max(r_k))


ax[0].plot([], [], c="k", ls="-", label=r"$dt = 10^{-2}$")
ax[0].plot([], [], c="k", ls="-.", label=r"$dt = 10^{-3}$")
ax[0].plot([], [], c="k", ls="--", label=r"$dt = 10^{-4}$")

ax[2].plot([], [], c="cadetblue", ls="-", label=r"leapfrog")
ax[2].plot([], [], c="firebrick", ls="-", label=r"rk4")
ax[2].legend()
ax[0].legend()

apo, peri = profile.estimate_apsides(x0, v0)

# put a fill between for apo and peri on ax[1]
ax[1].axhline(apo, color="green", label=r"$r_\mathrm{apo}$", alpha=0.3)
ax[1].axhline(peri, color="blue", label=r"$r_\mathrm{peri}$", alpha=0.3)
ax[1].legend()

ax[0].set_yscale("log")

ax[0].set_ylabel(r"$|1 - E/E_0|$")
ax[0].set_xlabel(r"$t/t_\mathrm{circ}(r)$")
ax[1].set_ylabel(r"$r/R_\mathrm{vir}$")
ax[1].set_xlabel(r"$t/t_\mathrm{circ}(r)$")
ax[2].set_xlabel(r"$x$")
ax[2].set_ylabel(r"$y$")
ax[2].set_aspect(1.0)

r_sample = np.linspace(0.1, 5, 100)
ax[3].plot(r_sample, profile.u_eff(r_sample, l0, offset=0))
# ax[3].plot(r_sample, profile.u_eff(r_sample, l0, offset=e0))
ax[3].axhline(e0, color="r", label=r"$E_0$")
ax[3].set_yscale("symlog")
ax[3].set_xlim(0.0, 5)
ax[3].legend()
ax[3].axvline(apo, ls="--", c="green", alpha=0.5)
ax[3].axvline(peri, ls="--", c="green", alpha=0.5)
ax[3].axhline(0.0, ls="--", c="k", alpha=0.5)
# ax[3].axhline(e0, ls='--', c='darkorange', alpha=.5)

for i in range(2):
    ax[3].axvline(r_peri[i], color=colors[integss[i]], alpha=1)
    ax[3].axvline(r_apo[i], color=colors[integss[i]], alpha=1)


ax[3].set_xlabel(r"$r/R_\mathrm{vir}$")
ax[3].set_ylabel(r"$U_\mathrm{eff}$")
fig.subplots_adjust(wspace=0.0)

ax[0].set_title("energy conservation")
ax[1].set_title("radius")
ax[2].set_title("trajectory")
ax[3].set_title("effective potential")
fig.suptitle("plummer")
fig.tight_layout()

plt.show()
