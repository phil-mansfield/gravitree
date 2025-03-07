import numpy as np
import matplotlib.pyplot as plt

n_str = "4"
pot_fn_fmt = "pot_table_at_n=%s_ic=%d_it=%d.dat"
acc_fn_fmt = "force_table_at_n=%s_ic=%d_it=%d.dat"
x_fmt = "einasto_n=%s_a=18.dat"

criteria_names = [
    "Barnes-Hut criterion",
    "PKDGRAV3 criterion",
    "Salmon-Warren criterion",
]
thetas = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0]
theta_names = [None, "0.01", "0.025", "0.05", "0.1", "0.25", "0.5", "1.0"]


def main():
    base_pot_fn = pot_fn_fmt % (n_str, 0, 0)
    base_acc_fn = acc_fn_fmt % (n_str, 0, 0)

    x = np.loadtxt(x_fmt % n_str)
    r = np.sqrt(np.sum(x**2, axis=1))

    pot_0 = np.loadtxt(base_pot_fn)
    acc_0 = np.loadtxt(base_acc_fn)
    acc_norm_0 = np.sqrt(np.sum(acc_0**2, axis=1))

    quantiles = np.array([0.5-0.95/2, 0.5-0.68/2, 0.5, 0.5+0.68/2, 0.5+0.95/2, 0.99])
    pot_f_err_table = np.zeros((len(criteria_names), len(thetas), len(quantiles)))
    acc_f_err_table = np.zeros((len(criteria_names), len(thetas), len(quantiles)))

    colors = [None, "r", "o", "navy", "g", "b", "purple", "k"]

    for ic in range(3):
        for it in range(1, 8):
            pot_fn = pot_fn_fmt % (n_str, ic, it)
            pot = np.loadtxt(pot_fn)
            pot_err = np.sqrt(np.sum((pot - pot_0) ** 2))
            pot_f_err = pot_err / np.abs(pot_0)

            acc_fn = acc_fn_fmt % (n_str, ic, it)
            acc = np.loadtxt(acc_fn)
            acc_norm = np.sqrt(np.sum(acc**2, axis=1))
            acc_err = np.sqrt(np.sum((acc_norm - acc_norm_0) ** 2))
            acc_f_err = acc_err / acc_norm_0
            for iq in range(len(quantiles)):
                pot_f_err_table[ic, it - 1, iq] = np.quantile(pot_f_err, quantiles[iq])
                acc_f_err_table[ic, it - 1, iq] = np.quantile(acc_f_err, quantiles[iq])

    
    colors = ["r", "orange", "b"]
    fig, ax = plt.subplots(1, 2, dpi=200, figsize=(10, 5))
    for ic in range(len(criteria_names)):
        
        ax[0].plot(
            thetas, pot_f_err_table[ic, :, 2], c=colors[ic], label=criteria_names[ic]
        )
        ax[0].plot(
            thetas, pot_f_err_table[ic, :, 5], c=colors[ic], ls='--'
        )
        ax[0].fill_between(
            thetas,
            pot_f_err_table[ic, :, 1],
            pot_f_err_table[ic, :, 3],
            color=colors[ic],
            alpha=0.2,
        )

        ax[1].plot(
            thetas, acc_f_err_table[ic, :, 2], c=colors[ic], label=criteria_names[ic]
        )
        ax[1].plot(
            thetas, acc_f_err_table[ic, :, 5], c=colors[ic], ls='--'
        )
        ax[1].fill_between(
            thetas,
            acc_f_err_table[ic, :, 1],
            acc_f_err_table[ic, :, 3],
            color=colors[ic],
            alpha=0.2,
        )
    for a in ax:
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel(r"$\theta$")

    ax[0].set_ylabel(r"$|\Delta \Phi|/|\Phi|$")
    ax[1].set_ylabel(r"$|\Delta \vec{a}|/|\vec{a}|$")
    ax[0].legend(frameon=True)
    fig.tight_layout()
    plt.savefig('eval_at_errors.png')
    plt.close()


if __name__ == "__main__":
    main()
