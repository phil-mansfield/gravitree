import numpy as np
import numpy.random as random
import sys
import abc
import scipy.special as special
import scipy.interpolate as interpolate

"""
To add a new profile model, do the following:

1) Make a new class that inherits from AbstractProfileShape. Just follow the
   pattern in the EinastoProfile instance and you'll be fine.
2) Add an entry to the models dictionary to parse_profile_name().
3) Update the comment in print_usage to include your new model.
"""

def print_usage():
	print("""gen_points usage intructions:

python3 gen_points.py <profile_name> <r_s> <n_points>

Parameters:
    profile_name - a comma-separated list containing the name of the profile
        shape and any secondary parameters (e.g. einasto,alpha=0.18).
        Unspecified parameters will be set to a default value.
    r_s - the radius where the logarithmic slope of the profile is -2.
    n_points - the number of points to generate. If a floating point number is
        given, it will be truncated to an integer.


Supported profile shapes:
    einasto - parameters: alpha (default=0.18)
""")


class AbstractProfileShape(abc.ABC):
	def variables(self):
		""" variables returns a list of profile secondary variables. These
		variables will be passed to the profile's constructor as keywords.
		"""

	def mass_fraction_enclosed(self, r_rs):
		""" mass_enclosed returns the fraction of the profile's mass enclosed
		within a given r/r_rs value, r_rs.
		"""

class EinastoProfile(AbstractProfileShape):
	def __init__(self, alpha=0.18):
		self.alpha=alpha

	def variables(self):
		return ["alpha"]

	def mass_fraction_enclosed(self, r_rs):
		return special.gammainc(3/self.alpha, 2*r_rs**self.alpha/self.alpha)


def parse_profile_name(arg):
	models = {
		"einasto": EinastoProfile
	}

	tokens = arg.split(",")
	model_name, variables = tokens[0], tokens[1:]
	kwargs = { }
	for v in variables:
		assert(v.count("=") == 1)
		name, value  = v.split("=")
		kwargs[name] = float(value)

	return models[model_name](**kwargs)
	
def random_angles(n):
	phi = 2*np.pi*random.random(n)
	theta = np.arccos(2*random.random(n) - 1)
	return phi, theta


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print_usage()
		exit(1)

	try:
		model = parse_profile_name(sys.argv[1])
		r_s = float(sys.argv[2])
		n_points = int(float(sys.argv[3]))
	except:
		print_usage()
		exit(1)

	r_rs = 10**np.linspace(-3, 3, 600)
	r_rs[0] = 0

	f = interpolate.interp1d(model.mass_fraction_enclosed(r_rs), r_rs)
	m_max = model.mass_fraction_enclosed(r_rs[-1])

	m = np.minimum(random.random(n_points), m_max)
	r = f(m) * r_s
	phi, theta = random_angles(n_points)
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)

	for i in range(len(x)):
		print("%.5f %.5f %.5f" % (x[i], y[i], z[i]))
