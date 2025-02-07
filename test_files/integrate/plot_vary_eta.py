import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams["font.family"] = "monospace"
from utils import get_trajectory, load_particles, get_quantiles

integ = "leapfrog"

dpss = ["0", "0.1", "0.5", "0.9"]
detas = ["0.025", "0.0025", "0.00025", "2.5e-05"]

dps = np.array(dpss, dtype=float)
etas = np.array(detas, dtype=float)

tfn = "_time.dat"
efn = "_energy_bf.dat"

colors = plt.cm.viridis(np.linspace(0, 1, len(dpss)))

fig, ax = plt.subplots(dpi=200)


# integs = ["leapfrog", "rk4"]

for k, _detas in enumerate(detas):
    de_50 = []
    de_16, de_84 = [], []

    for i, _dps in enumerate(dpss):
        directory = (
            f"vary_eta/snapshots_feps=1_int=lfadp_dt=1e-2_th=1e-2_dp={_dps}_eta={_detas}"
        )

        t_i = np.loadtxt(os.path.join(directory, tfn))
        e_i = np.loadtxt(os.path.join(directory, efn))

        e_i /= e_i[0]
        e_i = np.abs(1 - e_i)

        # masks time
        t_m = t_i >= (2 * 2*np.pi)

        p16, p50, p84 = get_quantiles(e_i)
        de_50.append(p50)
        de_16.append(p16)
        de_84.append(p84)

    ax.plot(dps, de_50, color=colors[k], marker="o", markersize=0.5, label=r"$\eta $ = " + _detas)
    ax.fill_between(dps, de_16, de_84, color=colors[k], alpha=0.2)

ax.legend()
ax.set_yscale("log")

# ax.set_xscale("symlog")
# ax.set_xlim(-.1, 2e3)

ax.set_xlabel(r"$\Delta p / R_\mathrm{vir}$")
ax.set_ylabel(r"$|1 - \frac{E}{E_0}|$")

fig.tight_layout()
ax.set_title("energy error after 2 orbits")
plt.show()
