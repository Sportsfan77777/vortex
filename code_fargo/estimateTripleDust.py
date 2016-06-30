"""
Uses Lyra+Lin 2013 (http://iopscience.iop.org/article/10.1088/0004-637X/775/1/17/pdf)
to calculate dust distribution
"""

import sys

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import patches

# Constant
default_S = 3

if len(sys.argv) > 1:
	default_S = int(sys.argv[1])

def max_y(s = default_S):
	if s == 1:
		return 10
	elif s == 3:
		return 30
	elif s == 5:
		return 50
	elif s == 10:
		return 120

# Helper Functions

def semi_minor(angle, aspect_ratio, radius):
	return radius * (angle * (np.pi / 180.0)) * 0.06 * 2 #aspect_ratio

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

def get_dust(x, aspect_ratio, max_density, S = default_S):
	f_sq = scale_function_sq(aspect_ratio)
	#print f_sq

	coeff = max_density * (S + 1)**(1.5)

	argument = x**2 * f_sq
	exp = np.exp(-(S + 1) * argument / 2.0)

	return coeff * exp

ten_density = 2.8
ten_radius_over_dr = 1.7 / 0.3
ten_extent = 120.0 * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

thousand_density = 1.85
thousand_radius_over_dr = 1.5 / 0.25
thousand_extent = 240.0 * (np.pi / 180)
thousand_aspect_ratio = thousand_radius_over_dr * thousand_extent

print thousand_aspect_ratio

ten_dust = lambda x : get_dust(x, ten_aspect_ratio, ten_density)
thousand_dust = lambda x : get_dust(x, thousand_aspect_ratio, thousand_density)

num_xs = 500

ten_xs = np.linspace(-180, 180, num_xs)
ten_ys = [ten_dust(semi_minor(x, ten_aspect_ratio / 2.0, ten_radius_over_dr)) for x in ten_xs]

thousand_xs = np.linspace(-180, 180, num_xs)
thousand_ys = [thousand_dust(semi_minor(x, thousand_aspect_ratio / 2.0, thousand_radius_over_dr)) for x in thousand_xs]

### PLOTTING ###

# Parameters
linewidth = 4
fontsize = 16

alpha = 0.3

my_dpi = 100

def add_to_plot(S):
	pass

def make_plot(aspect_ratios):
	small, large = aspect_ratios
	S_range = 3 * np.array([10**(-3), 10**(0), 10**(1)]) # um, mm, cm

	# Figure
	fig = plot.figure()
	ax = fig.add_subplot(111)

	# Identify vortex range
	extent1 = 60
	start1 = np.searchsorted(ten_xs, -extent1)
	end1 = np.searchsorted(ten_xs, extent1)

	extent2 = 120
	start2 = np.searchsorted(thousand_xs, -extent2)
	end2 = np.searchsorted(thousand_xs, extent2)

	# Curves
	plot.plot(ten_xs, ten_ys, c = "b", alpha = alpha, linewidth = linewidth, linestyle = "--")
	plot.plot(thousand_xs, thousand_ys, c = "g", alpha = alpha, linewidth = linewidth, linestyle = "--")

	plot.plot(ten_xs[start1 : end1], ten_ys[start1 : end1], c = "b", linewidth = linewidth, label = r"$T_{growth} = 10$ $\rm{orbits}$")
	plot.plot(thousand_xs[start2 : end2], thousand_ys[start2 : end2], c = "g", linewidth = linewidth, label = r"$T_{growth} = 1000$ $\rm{orbits}$")

	# Axes
	y_max = max_y(default_S)
	plot.xlim(-180, 180)
	plot.ylim(0, y_max)

	angles = np.linspace(-180, 180, 7)
	plot.xticks(angles)

	# Hatched Regions
	plot.text(0, 7, "Vortex Extents", style = "italic", fontsize = fontsize - 2, horizontalalignment = "center")

	ellipse1 = patches.Ellipse(xy = (0, 5), width = 2 * extent1, height = 1.5, color = "b", alpha = alpha, hatch = "+")
	ax.add_artist(ellipse1)

	ellipse2 = patches.Ellipse(xy = (0, 2.5), width = 2 * extent2, height = 1.5, color = "g", alpha = alpha, hatch = "X")
	ax.add_artist(ellipse2)


	#plot.fill_between([-60, 60, 60, -60], [3, 3, 4.5, 4.5], hatch = "+", color = "b", edgecolor = "b", alpha = alpha, linewidth = 0)
	#plot.fill_between([-120, 120, 120, -120], [0.5, 0.5, 2, 2], hatch = "X", color = "g", edgecolor = "g", alpha = alpha, linewidth = 0)

	# Annotate
	plot.xlabel(r"Angular Distance from Vortex Center: $\phi - \phi_c$", fontsize = fontsize)
	plot.ylabel("Normalized Dust Density", fontsize = fontsize)
	plot.title("Steady-state (mm-size) Dust Distributions for\n" + r"$m_p = 1$ $M_J$, $\nu_{disk} = 10^{-7}$", fontsize = fontsize + 1)
	plot.legend(loc = "upper right")

	# Save + Show
	plot.savefig("triple_dust_comparison.png", bbox_inches = "tight", dpi = my_dpi)
	plot.show()

# Set up figure

### Plot Sample ###
# Full Range
aspect_ratios = []

scale_height = 0.06

### First Case: 1 MJ, v = 10^-7 ###
ten_radius = 1.7
ten_scale_height = 0.06 * ten_radius
ten_density = 2.7 / ten_scale_height
ten_radius_over_dr = ten_radius / 0.3
ten_extent = 120.0 * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

thousand_radius = 1.5
thousand_scale_height = scale_height * thousand_radius
thousand_density = 1.68 / thousand_scale_height
thousand_radius_over_dr = thousand_radius / 0.25
thousand_extent = 240.0 * (np.pi / 180)
thousand_aspect_ratio = thousand_radius_over_dr * thousand_extent

print thousand_aspect_ratio

aspect_ratios.append((ten_aspect_ratio, thousand_aspect_ratio))

### Second Case: 5 MJ, v = 10^-7 ###

ten_radius = 1.8
ten_scale_height = scale_height * ten_radius
ten_density = 3.6 / ten_scale_height
ten_radius_over_dr = ten_radius / 0.36
ten_extent = 120.0 * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

four_thousand_radius = 1.7
four_thousand_scale_height = scale_height * four_thousand_radius
four_thousand_density = 1.94 / four_thousand_scale_height
four_thousand_radius_over_dr = four_thousand_radius / 0.24
four_thousand_extent = 180.0 * (np.pi / 180)
four_thousand_aspect_ratio = four_thousand_radius_over_dr * four_thousand_extent

print thousand_aspect_ratio

aspect_ratios.append((ten_aspect_ratio, four_thousand_aspect_ratio))

### Third Case: 5 MJ, v = 10^-6 ###

ten_radius = 1.7
ten_scale_height = scale_height * ten_radius
ten_density = 2.82 / ten_scale_height
ten_radius_over_dr = ten_radius / 0.3
ten_extent = 90.0 * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

# taken at 750
thousand_radius = 1.6
thousand_scale_height = scale_height * thousand_radius
thousand_density = 1.84 / thousand_scale_height
thousand_radius_over_dr = thousand_radius / 0.18
thousand_extent = 180.0 * (np.pi / 180)
thousand_aspect_ratio = thousand_radius_over_dr * thousand_extent

print thousand_aspect_ratio

aspect_ratios.append((ten_aspect_ratio, thousand_aspect_ratio))


fig = plot.figure()
for S in S_range:
    add_to_plot(S)

make_plot(show = True)
plot.close(fig)

