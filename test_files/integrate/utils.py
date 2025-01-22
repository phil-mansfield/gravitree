import numpy as np

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

def get_quantiles(data, axis=0):
    p50 = np.quantile(data, .50, axis)
    p84 = np.quantile(data, .84, axis)
    p16 = np.quantile(data, .16, axis)
    return p16, p50, p84

# source: https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
def win_avg(x, n):
    return np.convolve(x, np.ones(n)/n, mode='valid')