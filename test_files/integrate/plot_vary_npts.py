import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams["font.family"] = "monospace"
from utils import get_trajectory, load_particles, get_quantiles

nptss = ["3", "4", "5"]
npts = np.array(nptss, dtype=float)
tfn = "circ_t.dat"
efn = "circ_bf_e.dat"

fig, ax = plt.subplots(dpi=200)

colors = ["cadetblue", "darkorange", "firebrick"]

# |delta E / E|

integs = ['leapfrog', 'rk4']

for k, integ in enumerate(integs):
    de_50 = []
    de_16, de_84 = [], []

    for i, npt in enumerate(nptss):
        directory = f"vary_npts/snapshots_n={npt}_feps=1_int={integ}_dt=1e-3_th=1e-2"
        t_i = np.loadtxt(os.path.join(directory, tfn))
        e_i = np.loadtxt(os.path.join(directory, efn))
        e_i /= e_i[0]
        e_i = np.abs(1 - e_i)
        t_m = (t_i > 1.0) & ((t_i < 5.0))
        p16, p50, p84 = get_quantiles(e_i[t_m])
        de_50.append(p50)
        de_16.append(p16)
        de_84.append(p84)

    ax.plot(npts, de_50, color=colors[k], marker='s', label=integ)
    ax.fill_between(npts, de_16, de_84, color=colors[k], alpha=0.2)

ax.legend()
ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlabel(r'$n_\mathrm{pts}$')
ax.set_ylabel(r'$|1 - \frac{E}{E_0}|$')

fig.tight_layout()
ax.set_title('energy error after 2-5 orbits')
plt.show()
