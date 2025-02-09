import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from scipy.signal import argrelextrema
import argparse
import os

# usage: plot dipole simulated 1 unit apart, no softening, and with eta = 0.025
# > python3 plot_dipole.py 0.5 0 0.025
# usage: plot with 1 * eps softening
# > python3 plot_dipole.py 0.5 1 0.025

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
        file = fn.format(snap)
        pos = np.loadtxt(file)
        pos = pos.reshape(-1, 6)
        if not full:
            traj[k] = pos[idx, :3]
        else:
            traj[k] = pos[idx]
    return traj


def get_quantity(fn, idx, snaps):
    q = np.zeros(len(snaps))

    for k, snap in enumerate(snaps):
        file = fn.format(snap)
        dat = np.loadtxt(file)
        if dat.ndim == 0:
            dat = dat.reshape(
                1,
            )
        q[k] = dat[idx]

    return q


def estimateOrbitalPeriod(v_orb, snaps, time_units=False, dt=0.1):
    # arg rel extrema
    v = np.sqrt(np.sum(v_orb**2, axis=1))
    extr_idx = argrelextrema(v, np.less)
    estimated_period = snaps[extr_idx[0][0]]
    if time_units:
        estimated_period *= dt
    return estimated_period


def main(args):
    lim = 2
    # plot of a tracer over time
    tracer_indices = [0]
    import matplotlib.gridspec as gridspec

    fig = plt.figure(dpi=200)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1])

    ax = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_energy = fig.add_subplot(gs[1, :])

    axis_pairs = [[0, 1], [0, 2], [1, 2]]
    axis_labels = ['x', 'y', 'z']

    # integs = ["leapfrog", "rk4"]
    integs = ["leapfrog"]
    int_color = ["cadetblue", "firebrick"]
    int_ls = ["-", "-."]


    for a in ax:
        a.scatter([-float(args.dp), float(args.dp)], [0, 0], c='k', marker='o', s=1)

    for j, integ in enumerate(integs):

        sim_dir = f"vary_dipole/snapshots_feps={args.feps}_int={integ}_dt=1e-5_th=1e-2_dp={args.dp}_eta={args.eta}"
        files = glob(os.path.join(sim_dir, "snapshot*"))
        snaps = np.arange(0, len(files), 1)
        time = np.loadtxt(os.path.join(sim_dir, "_time.dat"))

        for i, idx in enumerate(tracer_indices):
            for k, a in enumerate(ax):
                traj_bf_i = get_trajectory(
                    os.path.join(sim_dir, "snapshot_bf_{0}.dat"), idx, snaps
                )

                a.plot(
                    traj_bf_i[:, axis_pairs[k][0]],
                    traj_bf_i[:, axis_pairs[k][1]],
                    color=int_color[j],
                    lw=0.5,
                    marker="s",
                    markersize=0.1,
                    alpha=0.2,
                    linestyle=int_ls[j],
                    label=f"{integ}",
                )

                a.set_xlabel(axis_labels[axis_pairs[k][0]])
                a.set_ylabel(axis_labels[axis_pairs[k][1]])

        for i, idx in enumerate(tracer_indices):
            e_exact = np.loadtxt(os.path.join(sim_dir, "_energy_bf.dat"))
            e_exact = e_exact.reshape(-1, 1)
            ax_energy.plot(time, np.abs(1 - e_exact[:, i] / e_exact[0, i]), color=int_color[j], alpha=0.2, linestyle=int_ls[j], label=f"{integ}")

    fig.suptitle(fr"$dp=: {float(args.dp):.3f}$")

    ax_energy.set_yscale("log")
    ax_energy.axhline(0, ls="--", alpha=0.1, color="k")
    ax_energy.set_ylabel(r"$|1 - E / E_0|$")
    ax_energy.set_xlabel("orbits")
    ax_energy.legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dp", type=str, help="Einasto profile to use")
    parser.add_argument("feps", type=str, default=0, help="softening scale")
    parser.add_argument("eta", type=float, default=0.025, help="adaptive eta")
    args = parser.parse_args()
    main(args)

