import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from profiles import Einasto, Plummer, PointMass


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
        file = fn % snap
        pos = load_particles(file)
        if not full:
            traj[k] = pos[idx, :3]
        else:
            traj[k] = pos[idx]
    return traj



def main():

    mpl.rcParams["font.family"] = "monospace"

    parser = argparse.ArgumentParser()
    parser.add_argument('profile', type=str, default='einasto', help='`plummer`, `einasto`, or `point`')
    parser.add_argument('vc', type=float, default=1.0, help='fraction of circular speed')
    args = parser.parse_args()
    
    profile_type = args.profile 
    vc = args.vc

    x0 = np.array([1.0, 0.0, 0.0])
    v0 = None
    e0 = None
    profile = None

    if vc not in [0.5, 0.05, 0.25, 1.0]:
        raise ValueError(f"Invalid vc value: {vc}. Allowed values are 0.5, 0.05, 0.25, or 1.0")

    if vc == 1.0:
        sim_dir = "snapshots"
    else:
        sim_dir = f"snapshots_{vc}"

    # NOTE: The parameters set in this are hard-coded in `simulate.go`
    if profile_type == 'plummer':
        profile = Plummer(0.05)
        v0 = np.array([0.0, vc * profile.get_circular_velocity(np.linalg.norm(x0)), 0.0])
        e0 = ((np.linalg.norm(v0) ** 2.0) / 2.0) + profile.plummer_potential(np.linalg.norm(x0))
    elif profile_type == 'einasto':
        profile = Einasto(0.1, 0.18)
        v0 = np.array([0.0, vc * profile.get_circular_velocity(np.linalg.norm(x0)), 0.0])
        e0 = ((np.linalg.norm(v0) ** 2.0) / 2.0) + profile.einasto_potential(np.linalg.norm(x0))
    elif profile_type == 'point':
        profile = PointMass(1.)
        v0 = np.array([0.0, vc * profile.get_circular_velocity(np.linalg.norm(x0)), 0.0])
        e0 = ((np.linalg.norm(v0) ** 2.0) / 2.0) + profile.point_mass_potential(np.linalg.norm(x0))
    else:
        raise ValueError(f"Unknown profile type: {profile_type}")

    l0 = np.linalg.norm(np.cross(x0, v0))
    apo, peri = profile.estimate_apsides(x0, v0)

    # plot    
    integss = ["leapfrog", "rk4", "lfadp"]
    dtss = ["1e-2", "1e-3", "1e-4"]

    tfn = "time.dat"
    efn = "energy.dat"
    afn = "acc.dat"
    fig, ax = plt.subplots(1, 4, figsize=(15, 5))
    colors = {"leapfrog": "cadetblue", "rk4": "firebrick", "lfadp": "purple"}
    ls = ["-", "-.", "--"]
    markers = {"leapfrog": "o", "rk4": "o", "lfadp": "^"}
    sizes = {"leapfrog": 0.1, "rk4": 0.1, "lfadp": 1.0}

    r_peri = []
    r_apo = []

    for i, integs in enumerate(integss):
        for j, dts in enumerate(dtss):
            directory = os.path.join(sim_dir, f"{profile_type}_int={integs}_dt={dts}")

            e = np.loadtxt(os.path.join(directory, efn))
            a = np.loadtxt(os.path.join(directory, afn))
            t = np.loadtxt(os.path.join(directory, tfn))

            e_err = np.abs(1 - e / e[0])
            ax[0].plot(t, e_err, color=colors[integs], ls=ls[j], alpha=0.2, marker=markers[integs], markersize=sizes[integs])
            traj_k = get_trajectory(
                os.path.join(directory, "t_%s.dat"), 0, np.arange(0, len(e), 1)
            )

            r_k = np.sqrt(np.sum(traj_k**2, axis=1))

            ax[1].plot(t, r_k, color=colors[integs], ls=ls[j], alpha=0.2, marker=markers[integs], markersize=sizes[integs])

            ax[2].plot(
                traj_k[:, 0], traj_k[:, 1], color=colors[integs], ls=ls[j], alpha=0.2, marker=markers[integs], markersize=sizes[integs]
            )

            r_apo.append(np.min(r_k))
            r_peri.append(np.max(r_k))


    ax[0].plot([], [], c="k", ls="-", label=r"$dt = 10^{-2}$")
    ax[0].plot([], [], c="k", ls="-.", label=r"$dt = 10^{-3}$")
    ax[0].plot([], [], c="k", ls="--", label=r"$dt = 10^{-4}$")

    ax[2].plot([], [], c="cadetblue", ls="-", marker="o", markersize=0.1, label=r"leapfrog")
    ax[2].plot([], [], c="firebrick", ls="-", marker="o", markersize=0.1, label=r"rk4")
    ax[2].plot([], [], c="purple", ls="-", marker="^", markersize=1.0, label=r"adaptive")
    ax[2].legend()
    ax[0].legend()

    apo, peri = profile.estimate_apsides(x0, v0)

    ax[1].axhline(apo, color="green", label=r"$r_\mathrm{apo}$", alpha=0.3)
    ax[1].axhline(peri, color="blue", label=r"$r_\mathrm{peri}$", alpha=0.3)
    ax[1].legend()

    ax[0].set_yscale("log")

    ax[0].set_ylabel(r"$|1 - E/E_0|$")
    ax[0].set_xlabel(r"$t/t_\mathrm{circ}(r)$")
    ax[1].set_ylabel(r"$r/R_\mathrm{vir}$")
    ax[1].set_xlabel(r"$t/t_\mathrm{circ}(r)$")
    ax[2].set_xlabel(r"$x$")
    ax[2].set_ylabel(r"$y$")
    ax[2].set_aspect(1.0)

    r_sample = np.linspace(0.1, 5, 100)
    ax[3].plot(r_sample, profile.u_eff(r_sample, l0, offset=0))
    ax[3].axhline(e0, color="r", label=r"$E_0$")
    ax[3].set_yscale("symlog")
    ax[3].set_xlim(0.0, 5)
    ax[3].legend()
    ax[3].axvline(apo, ls="--", c="green", alpha=0.5)
    ax[3].axvline(peri, ls="--", c="green", alpha=0.5)
    ax[3].axhline(0.0, ls="--", c="k", alpha=0.5)

    for i in range(2):
        ax[3].axvline(r_peri[i], color=colors[integss[i]], alpha=1)
        ax[3].axvline(r_apo[i], color=colors[integss[i]], alpha=1)


    ax[3].set_xlabel(r"$r/R_\mathrm{vir}$")
    ax[3].set_ylabel(r"$U_\mathrm{eff}$")
    fig.subplots_adjust(wspace=0.0)

    ax[0].set_title("energy conservation")
    ax[1].set_title("radius")
    ax[2].set_title("trajectory")
    ax[3].set_title("effective potential")
    fig.suptitle(fr"{profile_type}; $v_0 = {vc} \times v_\mathrm{{circ}}(r=1.0)$")
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()