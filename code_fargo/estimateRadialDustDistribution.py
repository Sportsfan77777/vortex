"""
Uses Lyra+Lin 2013 (http://iopscience.iop.org/article/10.1088/0004-637X/775/1/17/pdf)
to calculate dust distribution
"""

import numpy as np
from matplotlib import pyplot as plot

def calculate_xi(aspect_ratio):
	return 1 + aspect_ratio**(-2)

def calculate_vorticity(aspect_ratio):
	return 1.5 / (aspect_ratio - 1)

def scale_function_sq(aspect_ratio):
	xi = calculate_xi(aspect_ratio)
	vorticity = calculate_vorticity(aspect_ratio)

	first_term = 2.0 * vorticity * aspect_ratio
	second_term = xi**(-1) * (2 * (vorticity)**2 + 3)

	return first_term - second_term

def get_dust(x, aspect_ratio, max_density, S = 3):
	f_sq = scale_function_sq(aspect_ratio)
	print f_sq

	coeff = max_density * (S + 1)**(1.5)

	argument = x**2 * f_sq
	exp = np.exp(-(S + 1) * argument / 2.0)

	return coeff * exp

ten_density = 2.8
ten_aspect_ratio = 9.89

thousand_density = 1.85
thousand_aspect_ratio = 26.17

ten_dust = lambda x : get_dust(x, ten_aspect_ratio, ten_density)
thousand_dust = lambda x : get_dust(x, thousand_aspect_ratio, thousand_density)

xs = np.linspace(0, 5, 500)

ten_ys = [ten_dust(x) for x in xs]
thousand_ys = [thousand_dust(x) for x in xs]

### PLOTTING ###

# Parameters
linewidth = 4
fontsize = 16

my_dpi = 100

def make_plot(use_offset = False):

	# Curves
	if use_offset:
		offset = 0.56
	else:
		offset = 1.0
	plot.plot(xs * offset, ten_ys, linewidth = linewidth, label = r"$T_{growth} = 10$ $\rm{orbits}$")

	if use_offset:
		offset = 0.33
	else:
		offset = 1.0
	plot.plot(xs * offset, thousand_ys, linewidth = linewidth, label = r"$T_{growth} = 1000$ $\rm{orbits}$")

	# Axes
	if use_offset:
		plot.xlim(0, 2)
	else:
		plot.xlim(0, 3)
	plot.ylim(0, 30)

	# Annotate
	plot.xlabel(r"Distance from Vortex Center: $|r - r_c| / H$", fontsize = fontsize)
	plot.ylabel("Normalized Dust Density", fontsize = fontsize)
	plot.title("Comparing Dust Distributions for\n" + r"$m_p = 1$ $M_J$, $\nu_{disk} = 10^{-7}$", fontsize = fontsize + 1)
	plot.legend(loc = "upper right")

	# Save + Show
	plot.savefig("radial_dust_comparison.png", bbox_inches = "tight", dpi = my_dpi)
	plot.show()

make_plot()
