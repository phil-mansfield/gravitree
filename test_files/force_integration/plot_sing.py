import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from scipy.signal import argrelextrema
import argparse



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
        pos = load_particles(file)
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
            dat = dat.reshape(1,)
        q[k] = dat[idx]
    
    return q

def estimateOrbitalPeriod(v_orb, snaps, time_units=False, dt=.1):
    #arg rel extrema
    v = np.sqrt(np.sum(v_orb**2, axis=1))
    extr_idx = argrelextrema(v, np.less)
    estimated_period = snaps[extr_idx[0][0]]
    if time_units:
        estimated_period *= dt
    return estimated_period

def main():
    lim = 1
    # plot of a tracer over time
    tracer_indices = np.arange(0, 1)
    fig, ax = plt.subplots(1, 3, dpi=200, figsize=(12, 4))
    axis_pairs = [[0, 1], [0, 2], [1, 2]]

    for k, a in enumerate(ax):
        # a.set_xlim(.08, .12)
        a.set_xlim(-lim, lim)
        a.set_ylim(-lim, lim)

    snaps = np.arange(0, 100, 1)

    for i, idx in enumerate(tracer_indices):
        for k, a in enumerate(ax):
            traj_i = get_trajectory("snapshots/single_t_{0}.dat", i, snaps)

            a.plot(
                traj_i[:, axis_pairs[k][0]],
                traj_i[:, axis_pairs[k][1]],
                color="b",
                lw=0.5,
                marker="o",
                markersize=2,
                alpha=0.2,
                label="gt",
            )

            traj_bf_i = get_trajectory("snapshots/single_bf_t_{0}.dat", i, snaps)

            a.plot(
                traj_bf_i[:, axis_pairs[k][0]],
                traj_bf_i[:, axis_pairs[k][1]],
                color="r",
                lw=0.5,
                marker="s",
                markersize=2,
                alpha=0.2,
                label="bf",
            )

            a.set_xlabel(axis_pairs[k][0])
            a.set_ylabel(axis_pairs[k][0])

    fig.suptitle('single particle orbit')
    fig.tight_layout()
    plt.savefig("sing_orbit.pdf")
    plt.close()

    # plot energy conservation

    fig, ax = plt.subplots(dpi=200)
    for i, idx in enumerate(tracer_indices):

        e_approx = get_quantity("snapshots/single_e_{0}.dat", idx, snaps)
        e_exact  = get_quantity("snapshots/single_bf_e_{0}.dat",idx, snaps)
        ax.plot(snaps,
                np.abs(e_approx - e_exact) / np.abs(e_exact),
                color='r',
                alpha=.2
        )

    ax.set_yscale('symlog')

    print(np.abs(e_approx - e_exact) / np.abs(e_exact))
    ax.axhline(0, ls='--', alpha=.1, color='k')
    plt.savefig('sing_energy.pdf')

if __name__ == "__main__":
    main()
