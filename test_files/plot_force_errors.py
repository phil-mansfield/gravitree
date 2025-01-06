import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc

n_str = "4"
filename_fmt = "force_table_n=%s_ic=%d_it=%d.dat"
x_fmt = "einasto_n=%s_a=18.dat"

criteria_names = ["Barnes-Hut criterion", "PKDGRAV3 criterion",
	"Salmon-Warren criterion"]
thetas =[0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0]
theta_names = [None, "0.01", "0.025", "0.05", "0.1", "0.25", "0.5", "1.0"]


def main():
	palette.configure(False)

	base_filename = filename_fmt % (n_str, 0, 0)
	x = np.loadtxt(x_fmt % n_str)
	acc_0 = np.loadtxt(base_filename)
	acc_norm_0 = np.sqrt(np.sum(acc_0**2, axis=1))
	r = np.sqrt(np.sum(x**2, axis=1))

	quantiles = [0.5-0.95/2, 0.5-0.68/2, 0.5, 0.5+0.68/2, 0.5+0.95/2, 0.99]
	f_err_table = np.zeros(
		(len(criteria_names), len(thetas), len(quantiles))
	)

	colors = [None, pc("r"), pc("o"), pc("n"), pc("g"), pc("b"), pc("p"), pc("k")]

	for ic in range(3):
		plt.figure()
		for it in range(1, 8):
			filename = filename_fmt % (n_str, ic, it)
			acc = np.loadtxt(filename)
			acc_err = np.sqrt(np.sum((acc - acc_0)**2, axis=1))
			f_acc_err = acc_err / acc_norm_0
			
			plt.title(criteria_names[ic])
			plt.plot(r, f_acc_err, ".", alpha=0.2, c=colors[it])
			plt.plot([], [], ".", c=colors[it],
				label=r"$\theta = %s$" % theta_names[it])

			for iq in range(len(quantiles)):
				f_err_table[ic, it-1, iq] = np.quantile(f_acc_err, quantiles[iq])

		plt.title(criteria_names[ic])
		plt.xscale("log")
		plt.yscale("log")
		plt.xlabel(r"$r/R_{\rm vir}$")
		plt.ylabel(r"$|\Delta \vec{a}|/|\vec{a}|$")
		plt.legend(loc="lower left", frameon=True)
		plt.ylim(1e-12, 1)


	plt.figure()
	colors = [pc("r"), pc("o"), pc("b")]
	for ic in range(len(criteria_names)):
		plt.plot(thetas, f_err_table[ic,:,2], c=colors[ic], label=criteria_names[ic])
		plt.plot(thetas, f_err_table[ic,:,5], "--", lw=2, c=colors[ic])
		plt.fill_between(thetas, f_err_table[ic,:,1], f_err_table[ic,:,3],
			color=colors[ic], alpha=0.2)
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel(r"$\theta$")
	plt.ylabel(r"$|\Delta \vec{a}|/|\vec{a}|$")
	plt.legend(loc="lower right", frameon=True)
	plt.show()

if __name__ == "__main__": main()