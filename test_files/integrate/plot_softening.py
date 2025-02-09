import numpy as np
import matplotlib.pyplot as plt
import os

from utils import load_particles, get_trajectory, win_avg

snapshots = np.arange(0, 1000, 1)

fig, ax = plt.subplots(dpi=200)

softening = ["0", "1e-2","1","2","5","10", "50", "1e2"]
# colors = ['cadetblue', 'darkorange', 'firebrick']
# quantiles
_p50, _p16, _p84 = [], [], []

_p50_1o, _p16_1o, _p84_1o = [], [], []

for m, feps in enumerate(softening):
    # load orbits
    w0 = load_particles(f"snapshots_n=4_feps={feps}_int=leapfrog/circ_bf_t_0.dat")

    # scale time axes by average orbital time
    q0 = w0[:, :3]
    p0 = w0[:, 3:]
    r0 = np.sqrt(np.sum(q0**2, axis=1))
    v0 = np.sqrt(np.sum(p0**2, axis=1))
    t_circ = np.zeros((10))

    for k in range(10):
        t_circ[k] = 2 * np.pi * r0[k] / v0[k]
    
    dt_dsnap = 1e-4 * (1e6 / 1000)
    t = snapshots * dt_dsnap
    t_scaled = np.zeros((10, len(t)))

    for k in range(10):
        t_scaled[k, :] = t / t_circ[k]

    # load energy files
    e  = np.loadtxt(os.path.join(f"snapshots_n=4_feps={feps}_int=leapfrog", "circ_bf_e.dat"))

    # normalize by initial energy
    e = np.abs(1-e/e[0, :])

    # compute averages over a certain number
    # of orbits
    ts_min, ts_max = 2., 5.
    ef_avg = []
    ef_std = []

    ef_avg_1o = []
    ef_std_1o = []

    for i in range(10):
        # indices to compute average over
        indices = np.where((t_scaled[i] >= ts_min) & (t_scaled[i] < ts_max))[0]
        ef_avg.append(np.mean(e[indices, i], axis=0))
        ef_std.append(np.std(e[indices, i], axis=0)) 

        indices = np.where((t_scaled[i] >= 0.) & (t_scaled[i] < 1.))[0]
        ef_avg_1o.append(np.mean(e[indices, i], axis=0))
        ef_std_1o.append(np.std(e[indices, i], axis=0)) 
          
    _p50.append(np.percentile(ef_avg, 50))
    _p16.append(np.percentile(ef_avg, 16))
    _p84.append(np.percentile(ef_avg, 84))

    _p50_1o.append(np.percentile(ef_avg_1o, 50))
    _p16_1o.append(np.percentile(ef_avg_1o, 16))
    _p84_1o.append(np.percentile(ef_avg_1o, 84))



print(_p16, _p84)
ax.fill_between(np.array(softening, dtype=float) * .004, _p16, _p84, alpha=.2)
ax.plot(np.array(softening, dtype=float) * .004, _p50, marker='s', label=r'$t_\mathrm{min} = 2; t_\mathrm{max} = 5$')

ax.fill_between(np.array(softening, dtype=float) * .004, _p16_1o, _p84_1o, alpha=.2, color='red')
ax.plot(np.array(softening, dtype=float) * .004, _p50_1o, marker='s', color='red', label=r'$t_\mathrm{min} = 0; t_\mathrm{max} = 1$')

ax.legend()

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'$\langle |1 - E/E_0|\rangle|_{t_\mathrm{min}\leq t_\mathrm{orb} <t_\mathrm{max}}$')
ax.set_xlabel(r'$\epsilon/R_\mathrm{vir}$')

ax.set_title('energy error averaged between 2-5 orbits')

plt.show()