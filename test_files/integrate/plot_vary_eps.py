import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams["font.family"] = "monospace"
from utils import get_trajectory, load_particles, get_quantiles

integ = "leapfrog"
epss = ["0", "0.01", "0.1", "0.5", "0.75", "1.0"]
eps = np.array(epss, dtype=float)
tfn = "_time.dat"
efn = "_energy_bf.dat"

fig, ax = plt.subplots(dpi=200)

colors = ["cadetblue", "darkorange"]
ls = ['-', '-.']
size = [5, 2.5]

integs = ["leapfrog", "rk4"]

for k, integ in enumerate(integs):
    de_50 = []
    de_16, de_84 = [], []

    for i, _eps in enumerate(epss):
        directory = (
            f"vary_feps/snapshots_n=4_feps={epss[i]}_int={integ}_dt=1e-3_th=1e-2"
        )
        t_i = np.loadtxt(os.path.join(directory, tfn))
        e_i = np.loadtxt(os.path.join(directory, efn))

        e_i /= e_i[0]
        e_i = np.abs(1 - e_i)

        # masks time
        t_m = (t_i > 0) & (t_i <= 2.)

        p16, p50, p84 = get_quantiles(e_i[t_m])
        de_50.append(p50)
        de_16.append(p16)
        de_84.append(p84)

    ax.plot(eps, de_50, color=colors[k], marker="s", ls=ls[k], label=integ, markersize=size[k])
    ax.fill_between(eps, de_16, de_84, color=colors[k], alpha=0.2)

ax.legend()
ax.set_yscale("log")
ax.set_xscale("symlog")

# ax.set_xlim(-.1, 2e3)

ax.set_xlabel(r"$\epsilon / R_\mathrm{vir}$")
ax.set_ylabel(r"$|1 - \frac{E}{E_0}|$")

fig.tight_layout()
ax.set_title("energy error after 2-5 orbits")
plt.show()
