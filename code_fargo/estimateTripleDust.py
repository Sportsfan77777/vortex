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

### PLOTTING ###
# Options
growth = [1000, 4000, 1000]
dashes = [[8, 8], [10**5, 1], [10, 4]]

# Parameters
linewidth = 4
fontsize = 16

alpha = 0.3

my_dpi = 100

def add_to_plot(number, aspect_ratios, densities, r_ratios, str_r_ratios, angles):
    # Unwrap
    aspect1, aspect2 = aspect_ratios
    density1, density2 = densities
    r_over_dr1, r_over_dr2 = r_ratios
    str_r_over_dr1, str_r_over_dr2 = str_r_ratios
    extent1, extent2 = angles

    # Use these particle sizes
    S_range = 3 * np.array([10**(-3), 10**(0), 10**(1)]) # um, mm, cm

    # Subplot
    for S in S_range:
        # Distributions
        dust1 = lambda x : get_dust(x, aspect1, density1, S = S)
        dust2 = lambda x : get_dust(x, aspect2, density2, S = S)

        # Data
        num_xs = 500

        xs1 = np.linspace(-180, 180, num_xs)
        ys1 = [dust1(semi_minor(x, aspect1 / 2.0, r_over_dr1)) for x in xs1]

        xs2 = np.linspace(-180, 180, num_xs)
        ys2 = [dust1(semi_minor(x, aspect2 / 2.0, r_over_dr2)) for x in xs2]

        # Identify vortex range
        start1 = np.searchsorted(ten_xs, -extent1)
        end1 = np.searchsorted(ten_xs, extent1)

        start2 = np.searchsorted(thousand_xs, -extent2)
        end2 = np.searchsorted(thousand_xs, extent2)

        # Labels
        label1 = r"$S = " + str(S)+ r"$T_{growth} = 10$ $\rm{orbits}$"
        label2 = r"$S = " + str(S)+ r"$, $T_{growth} = " + str(growth[number]) + "$ $\rm{orbits}$"

        # Select Subplot
        ax = fig.add_subplot(1, 3, number + 1)

        # Curves
        plot.plot(xs1, ys1, c = "b", alpha = alpha, linewidth = linewidth, dashes = dashes[number])
        plot.plot(xs2, ys2, c = "g", alpha = alpha, linewidth = linewidth, dashes = dashes[number])

        plot.plot(xs1[start1 : end1], ys1[start1 : end1], c = "b", linewidth = linewidth, label = label1)
        plot.plot(xs2[start2 : end2], ys2[start2 : end2], c = "g", linewidth = linewidth, label = label2)    
    
    # Annotate
    plot.xlabel(r"$\phi - \phi_{center}$", fontsize = fontsize)
    plot.ylabel("Normalized Dust Density", fontsize = fontsize)
    plot.title("Dust Distributions for\n" + r"$m_p = 1$ $M_J$, $\nu_{disk} = 10^{-7}$", fontsize = fontsize + 1)
    plot.legend(loc = "lower right")

def make_plot():
    # Axes
    y_max = max_y(default_S)
    plot.xlim(-180, 180)
    #plot.ylim(0, 1000)

    angles = np.linspace(-180, 180, 7)
    plot.xticks(angles)

    # Save + Show
    plot.savefig("triple_dust_comparison.png", bbox_inches = "tight", dpi = my_dpi)
    plot.show()


##### Parameters #####

# Collections
aspect_ratios = []
densities = []
r_ratios = []
str_r_ratios = []
angles = []

# Constant
scale_height = 0.06

### First Case: 1 MJ, v = 10^-7 ###
ten_radius = 1.7
ten_dr = 0.3
ten_angle = 120

ten_scale_height = 0.06 * ten_radius
ten_density = 2.7 / ten_scale_height
ten_radius_over_dr = ten_radius / ten_dr
ten_extent = ten_angle * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

thousand_radius = 1.5
thousand_dr = 0.25
thousand_angle = 240

thousand_scale_height = scale_height * thousand_radius
thousand_density = 1.68 / thousand_scale_height
thousand_radius_over_dr = thousand_radius / thousand_dr
thousand_extent = thousand_angle * (np.pi / 180)
thousand_aspect_ratio = thousand_radius_over_dr * thousand_extent

print thousand_aspect_ratio

# Collect 1 #

aspect_ratios.append((ten_aspect_ratio, thousand_aspect_ratio))
densties.append((ten_density, thousand_density))
r_ratios.append((ten_radius_over_dr, thousand_radius_over_dr))
str_r_ratios.append(("%.1f / %.02f" % (ten_radius, ten_dr), ("%.1f / %.02f" % (thousand_radius, thousand_dr))))
angles.append((ten_angle, thousand_angle))

### Second Case: 5 MJ, v = 10^-7 ###
ten_radius = 1.8
ten_dr = 0.36
ten_angle = 120

ten_scale_height = scale_height * ten_radius
ten_density = 3.6 / ten_scale_height
ten_radius_over_dr = ten_radius / ten_dr
ten_extent = ten_angle * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

four_thousand_radius = 1.7
four_thousand_dr = 0.24
four_thousand_angle = 180

four_thousand_scale_height = scale_height * four_thousand_radius
four_thousand_density = 1.94 / four_thousand_scale_height
four_thousand_radius_over_dr = four_thousand_radius / four_thousand_dr
four_thousand_extent = four_thousand_angle * (np.pi / 180)
four_thousand_aspect_ratio = four_thousand_radius_over_dr * four_thousand_extent

print thousand_aspect_ratio

# Collect 2 #

aspect_ratios.append((ten_aspect_ratio, four_thousand_aspect_ratio))
densties.append((ten_density, four_thousand_density))
r_ratios.append((ten_radius_over_dr, four_thousand_radius_over_dr))
str_r_ratios.append(("%.1f / %.02f" % (ten_radius, ten_dr), ("%.1f / %.02f" % (four_thousand_radius, four_thousand_dr))))
angles.append((ten_angle, four_thousand_angle))

### Third Case: 5 MJ, v = 10^-6 ###
ten_radius = 1.7
ten_dr = 0.3
ten_angle = 90

ten_scale_height = scale_height * ten_radius
ten_density = 2.82 / ten_scale_height
ten_radius_over_dr = ten_radius / ten_dr
ten_extent = ten_angle * (np.pi / 180)
ten_aspect_ratio = ten_radius_over_dr * ten_extent

print ten_aspect_ratio

# taken at 750
thousand_radius = 1.6
thousand_dr = 0.18
thousand_angle = 180

thousand_scale_height = scale_height * thousand_radius
thousand_density = 1.84 / thousand_scale_height
thousand_radius_over_dr = thousand_radius / thousand_dr
thousand_extent = thousand_angle * (np.pi / 180)
thousand_aspect_ratio = thousand_radius_over_dr * thousand_extent

print thousand_aspect_ratio

# Collect 3 #

aspect_ratios.append((ten_aspect_ratio, thousand_aspect_ratio))
densties.append((ten_density, thousand_density))
r_ratios.append((ten_radius_over_dr, thousand_radius_over_dr))
str_r_ratios.append(("%.1f / %.02f" % (ten_radius, ten_dr), ("%.1f / %.02f" % (thousand_radius, thousand_dr))))
angles.append((ten_angle, thousand_angle))

##### Set up plot #####

fig = plot.figure()

for i, (aspect_ratios_i, densities_i, r_ratios_i, str_r_ratios_i, angles_i) in enumerate(zip(aspect_ratios, densities, r_ratios, str_r_ratios, angles)):
    add_to_plot(i, aspect_ratios_i, densities_i, r_ratios_i, str_r_ratios_i, angles_i)

make_plot()
plot.close(fig)

