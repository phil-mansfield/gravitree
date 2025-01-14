import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from scipy.signal import argrelextrema


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
        if not full:
            traj[k] = pos[idx, :3]
        else:
            traj[k] = pos[idx]
    return traj

def estimateOrbitalPeriod(v_orb, snaps, time_units=False, dt=.1):
    #arg rel extrema
    v = np.sqrt(np.sum(v_orb**2, axis=1))
    extr_idx = argrelextrema(v, np.less)
    estimated_period = snaps[extr_idx[0][0]]
    if time_units:
        estimated_period *= dt
    return estimated_period

def main():
    halo = load_particles(f"../einasto_n=3_a=18.dat")
    lim = 55
    # plot of a tracer over time
    tracer_indices = np.arange(0, 10)
    fig, ax = plt.subplots(1, 3, dpi=200, figsize=(12, 4))
    axis_pairs = [[0, 1], [0, 2], [1, 2]]

    for k, a in enumerate(ax):
        a.set_xlim(-lim, lim)
        a.set_ylim(-lim, lim)
        a.scatter(
            halo[:, axis_pairs[k][0]],
            halo[:, axis_pairs[k][1]],
            c="k",
            alpha=0.2,
            s=0.1,
        )

    snaps = np.arange(0, 2000, 10)

    for i, idx in enumerate(tracer_indices):
        for k, a in enumerate(ax):
            traj_i = get_trajectory("snapshots/circ_t_{0}.dat", i, snaps)

            a.plot(
                traj_i[:, axis_pairs[k][0]],
                traj_i[:, axis_pairs[k][1]],
                color="b",
                lw=0.5,
                marker="o",
                markersize=0.1,
                alpha=0.2,
                label="gt",
            )

            traj_bf_i = get_trajectory("snapshots/circ_bf_t_{0}.dat", i, snaps)

            a.plot(
                traj_bf_i[:, axis_pairs[k][0]],
                traj_i[:, axis_pairs[k][1]],
                color="r",
                lw=0.5,
                marker="s",
                markersize=0.1,
                alpha=0.2,
                label="bf",
            )

            a.set_xlabel(axis_pairs[k][0])
            a.set_ylabel(axis_pairs[k][0])

    fig.suptitle('circular orbit')
    fig.tight_layout()
    plt.savefig("circular_orbit.pdf")
    plt.close()


if __name__ == "__main__":
    main()
