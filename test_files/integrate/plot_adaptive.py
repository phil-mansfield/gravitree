import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from glob import glob

# usage: plots outputs of ./run_vary_eta.sh
# > python3 plot_adaptive.py

mpl.rcParams["font.family"] = "monospace"
from utils import get_trajectory, load_particles, get_quantiles

integ = "leapfrog"

detas = ["0.025", "0.0025", "0.00025", "2.5e-05"]

etas = np.array(detas, dtype=float)

dtfn = "_dt_bf.dat"
tfn = "_time.dat"
efn = "_energy_bf.dat"

colors = plt.cm.viridis(np.linspace(0, 1, len(detas)))

fig, ax = plt.subplots(dpi=200)

for _eta in detas:
    directory = (
                f"vary_eta/snapshots_feps=1_int=lfadp_dt=1e-2_th=1e-2_dp=0.5_eta={_eta}"
            )
    
    files = glob(os.path.join(directory, "snapshot*"))
    snaps = np.arange(0, len(files), 1)
    
    traj = get_trajectory(os.path.join(directory, "snapshot_bf_%i.dat"), 0, snaps, full=False)
    r = np.sqrt(np.sum(traj**2, axis=1))

    r_sort = np.argsort(r)

    dt = np.loadtxt(os.path.join(directory, dtfn))

    r_sorted = r[r_sort]
    t_sorted = dt[r_sort]

    # ax.scatter(r_sorted, t_sorted, label=f'eta={_eta}', color=colors[detas.index(_eta)])
    bin_width = 0.1
    bins = np.arange(0, np.max(r_sorted) + bin_width, bin_width)
    bin_indices = np.digitize(r_sorted, bins)

    bin_medians = []
    bin_perc_16 = []
    bin_perc_84 = []
    bin_centers = []

    for i in range(1, len(bins)):
        bin_data = t_sorted[bin_indices == i]
        if len(bin_data) > 0:
            bin_medians.append(np.median(bin_data))
            bin_perc_16.append(np.percentile(bin_data, 16))
            bin_perc_84.append(np.percentile(bin_data, 84))
            bin_centers.append((bins[i] + bins[i - 1]) / 2)

    ax.plot(bin_centers, bin_medians, label=fr'$\eta={_eta}$', color=colors[detas.index(_eta)])
    ax.fill_between(bin_centers, bin_perc_16, bin_perc_84, color=colors[detas.index(_eta)], alpha=0.3)

ax.legend()
ax.set_yscale('log')
ax.set_ylabel(r'$\Delta t_\mathrm{grav}$')
ax.set_xlabel(r'$r/R_\mathrm{vir}$')
# integs = ["leapfrog", "rk4"]

fig, ax = plt.subplots(dpi=200)

for _eta in detas:
    directory = (
                f"vary_eta/snapshots_feps=1_int=lfadp_dt=1e-2_th=1e-2_dp=0.5_eta={_eta}"
            )
    
    e = np.loadtxt(os.path.join(directory, efn))
    dt = np.loadtxt(os.path.join(directory, dtfn))

    e0 = e[0]
    e_diff = np.abs(1 - e / e0)

    ax.plot(dt, e_diff, label=fr'$\eta={_eta}$', color=colors[detas.index(_eta)])

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\Delta t_\mathrm{grav}$')
ax.set_ylabel(r'$|1 - E/E_0|$')

fig, ax = plt.subplots(dpi=200)

for _eta in detas:

    directory = (
                f"vary_eta/snapshots_feps=1_int=lfadp_dt=1e-2_th=1e-2_dp=0.5_eta={_eta}"
            )
    
    files = glob(os.path.join(directory, "snapshot*"))
    snaps = np.arange(0, len(files), 1)
    
    traj = get_trajectory(os.path.join(directory, "snapshot_bf_%i.dat"), 0, snaps, full=True)
    
    v = np.sqrt(np.sum(traj[:, 3:]**2, axis=1))  # Calculate velocity
    e = np.loadtxt(os.path.join(directory, efn))
    e0 = e[0]
    e_diff = np.abs(1 - e / e0)

    ax.plot(v, e_diff, label=fr'$\eta={_eta}$', color=colors[detas.index(_eta)])

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$v$')
ax.set_ylabel(r'$|1 - E/E_0|$')
fig, ax = plt.subplots(dpi=200)

for _eta in detas:

    directory = (
                f"vary_eta/snapshots_feps=1_int=lfadp_dt=1e-2_th=1e-2_dp=0.5_eta={_eta}"
            )
    
    files = glob(os.path.join(directory, "snapshot*"))
    snaps = np.arange(0, len(files), 1)
    
    traj = get_trajectory(os.path.join(directory, "snapshot_bf_%i.dat"), 0, snaps, full=True)
    
    v = np.sqrt(np.sum(traj[:, 3:]**2, axis=1))  # Calculate velocity
    dt = np.loadtxt(os.path.join(directory, dtfn))
    e = np.loadtxt(os.path.join(directory, efn))
    e0 = e[0]
    e_diff = np.abs(1 - e / e0)

    sc = ax.scatter(v, dt, c=np.log10(e_diff), cmap='viridis', label=fr'$\eta={_eta}$', alpha=0.7)

cbar = plt.colorbar(sc, ax=ax)
cbar.set_label(r'$\log_{10}(|1 - E/E_0|)$')

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$v$')
ax.set_ylabel(r'$\Delta t_\mathrm{grav}$')



fig, ax = plt.subplots(dpi=200)

for _eta in detas:
    directory = (
                f"vary_eta/snapshots_feps=1_int=lfadp_dt=1e-2_th=1e-2_dp=0.5_eta={_eta}"
            )
    
    dt = np.loadtxt(os.path.join(directory, dtfn))
    acc = np.loadtxt(os.path.join(directory, "_acc.dat"))

    ax.scatter(acc, dt, label=fr'$\eta={_eta}$', color=colors[detas.index(_eta)])

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'$\Delta t_\mathrm{grav}$')
ax.set_xlabel(r'Acceleration')
plt.show()