import numpy as np
import numpy.random as random
import sys
import abc
import scipy.special as special
import scipy.interpolate as interpolate

from gen_points import random_angles

def print_usage():
	print("""gen_points usage intructions:

python3 gen_shell.py <r> <n_points>

Parameters:
    r - the radius of the shell.
    n_points - the number of points to generate. If a floating point number is
        given, it will be truncated to an integer.
""")


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print_usage()
		exit(1)
	
	try:
		r = float(sys.argv[1])
		n_points = int(float(sys.argv[2]))
	except:
		print_usage()
		exit(1)

	phi, theta = random_angles(n_points)
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)

	for i in range(len(x)):
		print("%.5f %.5f %.5f" % (x[i], y[i], z[i]))
