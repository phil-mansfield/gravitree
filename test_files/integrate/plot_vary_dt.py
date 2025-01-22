import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams["font.family"] = "monospace"
from utils import get_trajectory, load_particles, get_quantiles

integ = "leapfrog"
dts = ["1e-1", "1e-2", "1e-3", "1e-4"]
_dts = np.array(dts, dtype=float)
tfn = "point_t.dat"
efn = "point_bf_e.dat"

fig, ax = plt.subplots(dpi=200)

colors = ["cadetblue", "darkorange"]

# |delta E / E|

integs = ['leapfrog', 'rk4']

for k, integ in enumerate(integs):
    de_50 = []
    de_16, de_84 = [], []

    for i, dt in enumerate(dts):
        directory = f"vary_dt/snapshots_int={integ}_dt=%s"
        t_i = np.loadtxt(os.path.join(directory % dt, tfn))
        e_i = np.loadtxt(os.path.join(directory % dt, efn))
        e_i /= e_i[0]
        e_i = np.abs(1 - e_i)
        t_m = (t_i > 1.0) & ((t_i < 5.0))
        p16, p50, p84 = get_quantiles(e_i[t_m])
        de_50.append(p50)
        de_16.append(p16)
        de_84.append(p84)

    ax.plot(_dts / (2 * np.pi), de_50, color=colors[k], marker='s', label=integ)
    ax.fill_between(_dts / (2 * np.pi), de_16, de_84, color=colors[k], alpha=0.2)

ax.legend()
ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlabel(r'$\Delta t / t_\mathrm{circ}$')
ax.set_ylabel(r'$|1 - \frac{E}{E_0}|$')

fig.tight_layout()
ax.set_title('energy error after 2-5 orbits')
plt.show()
