import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from scipy.signal import argrelextrema
import argparse
import os


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

def main(args):
    halo = load_particles(f"../einasto_n={args.n}_a=18.dat")
    lim = 2
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

    r_halo = np.sqrt(np.sum(halo**2, axis=1))
    print(f"n_part(x<1){np.sum(r_halo <= 1.0)}")

    snaps = np.arange(0, 1000, 1)

    for i, idx in enumerate(tracer_indices):
        for k, a in enumerate(ax):
            traj_i = get_trajectory(os.path.join(f"snapshots_n={args.n}", "circ_t_{0}.dat"), i, snaps)

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

            traj_bf_i = get_trajectory(os.path.join(f"snapshots_n={args.n}", "circ_bf_t_{0}.dat"), idx, snaps)

            a.plot(
                traj_bf_i[:, axis_pairs[k][0]],
                traj_bf_i[:, axis_pairs[k][1]],
                color="r",
                lw=0.5,
                marker="s",
                markersize=0.1,
                alpha=0.2,
                label="bf",
            )

            a.set_xlabel(axis_pairs[k][0])
            a.set_ylabel(axis_pairs[k][1])

    fig.suptitle('circular orbit')
    fig.tight_layout()
    plt.savefig(f"circular_orbit_n={args.n}.pdf")
    plt.close()

   # plot energy conservation

    fig, ax = plt.subplots(dpi=200)
    for i, idx in enumerate(tracer_indices):

        e_exact  = np.loadtxt(os.path.join(f"snapshots_n={args.n}", "circ_bf_e.dat"))
        e_approx  = np.loadtxt(os.path.join(f"snapshots_n={args.n}", "circ_e.dat"))
        ax.plot(snaps,
                e_exact[:, i] / e_exact[0, i],
                color='r',
                alpha=.2
        )

        ax.plot(snaps,
                e_approx / e_approx[0],
                color='r',
                alpha=.2
        )

    ax.set_yscale('log')
    
    ax.axhline(0, ls='--', alpha=.1, color='k')
    fig.tight_layout()
    ax.set_ylabel(r'$|E| / |E_\mathrm{init}|$')
    ax.set_xlabel('snapshot')
    plt.savefig(f'circ_energy_n={args.n}.pdf')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int, help='Einasto profile to use')
    args = parser.parse_args()
    main(args)
