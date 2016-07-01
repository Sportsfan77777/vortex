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

# Helper Functions

def resort(arr):
    print arr
    print len(arr)
    return [arr[0], arr[2], arr[4], arr[1], arr[3], arr[5]]

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
dashes = [[8, 6, 4, 6], [10**5, 1], [3, 3]]

mass = [1, 5, 5]
viscosity = [-7, -7, -6]

sizes = [r"$\rm{cm}$", r"$\rm{mm}$", r"$\rm{\mu m}$"]

# Parameters
linewidths = [4, 5, 4]
fontsize = 16

alpha = 0.3

my_dpi = 100

def add_to_plot(number, aspect_ratios, densities, r_ratios, rs, drs, angles):
    # Unwrap
    aspect1, aspect2 = aspect_ratios
    density1, density2 = densities
    r_over_dr1, r_over_dr2 = r_ratios
    r1, r2 = rs
    dr1, dr2 = drs
    extent1, extent2 = angles

    print density1, density2

    # Use these particle sizes
    S_range = 3 * np.array([10**(1), 10**(0), 10**(-3)]) # cm, mm, um

    # Subplot
    for i, S in enumerate(S_range):
        # Distributions
        dust1 = lambda x : get_dust(x, aspect1, density1, S = S)
        dust2 = lambda x : get_dust(x, aspect2, density2, S = S)

        # Data
        num_xs = 500

        xs1 = np.linspace(0, 180, num_xs)
        ys1 = [dust1(semi_minor(x, aspect1 / 2.0, r_over_dr1)) for x in xs1]

        xs2 = np.linspace(0, 180, num_xs)
        ys2 = [dust2(semi_minor(x, aspect2 / 2.0, r_over_dr2)) for x in xs2]

        # Identify vortex range
        start1 = np.searchsorted(xs1, -extent1 / 2)
        end1 = np.searchsorted(xs1, extent1 / 2)

        start2 = np.searchsorted(xs2, -extent2 / 2)
        end2 = np.searchsorted(xs2, extent2 / 2)

        # Labels
        label1 = sizes[i] #+ r"$T_{growth} = 10$ $\rm{orbits}$"
        label2 = sizes[i] #+ r"$, $T_{growth} = " + str(growth[i]) + r"$ $\rm{orbits}$"

        # Select Subplot
        ax = fig.add_subplot(1, 3, number + 1)

        # Curves
        #plot.plot(xs1, ys1, c = "b", alpha = alpha, linewidth = linewidths[i], dashes = dashes[i])
        #plot.plot(xs2, ys2, c = "g", alpha = alpha, linewidth = linewidths[i], dashes = dashes[i])

        plot.plot(xs1[start1 : end1], ys1[start1 : end1], c = "b", linewidth = linewidths[i], dashes = dashes[i], label = label1, zorder = 5)
        plot.plot(xs2[start2 : end2], ys2[start2 : end2], c = "g", linewidth = linewidths[i], dashes = dashes[i], label = label2, zorder = 1)

    # Axes
    plot.xlim(0, 180)
    plot.ylim(10**1, 10**4)
    plot.yscale("log")

    if number == 2:
        angles = np.linspace(0, 120, 5)
    else:
        angles = np.linspace(0, 120, 5)
    plot.xticks(angles)

    if number != 0:
        # Remove unless 1st frame
        ax.set_yticklabels([])

    # Annotate
    margin = 128
    diff = 5

    top = 3.9
    base1 = 2.9
    split = 2.73
    base2 = 1.7
    bottom = 1.55

    # Vertical
    plot.plot([margin - diff, margin - diff], [10**(bottom), 10**(top)], c = "k")

    plot.plot([margin - diff, 180], [10**(top), 10**(top)], c = "k") # Horizontal

    plot.text(margin, 10**(base1 + 0.8), r"$T_{growth} = $" + "%d" % 10)
    plot.text(margin, 10**(base1 + 0.6), r"$\chi = $ " + "%.1f" % aspect1)
    plot.text(margin, 10**(base1 + 0.4), r"$r = $" + "%.1f" % (r1))
    plot.text(margin, 10**(base1 + 0.2), r"$dr = $" + "%.02f" % (dr1))
    plot.text(margin, 10**(base1), r"$\rho_{peak} = $" + "%.1f" % density1)

    plot.plot([margin - diff, 180], [10**(split), 10**(split)], c = "k") # Horizontal

    plot.text(margin, 10**(base2 + 0.8), r"$T_{growth} = $" + "%d" % growth[number])
    plot.text(margin, 10**(base2 + 0.6), r"$\chi = $ " + "%.1f" % aspect2)
    plot.text(margin, 10**(base2 + 0.4), r"$r = $" + "%.1f" % (r2))
    plot.text(margin, 10**(base2 + 0.2), r"$dr = $" + "%.02f" % (dr2))
    plot.text(margin, 10**(base2), r"$\rho_{peak} = $" + "%.1f" % density2)

    plot.plot([margin - diff, 180], [10**(bottom), 10**(bottom)], c = "k") # Horizontal

    plot.xlabel(r"$|\phi - \phi_{center}|$ $\rm{(in}$ $\rm{degrees)}$", fontsize = fontsize)
    if number == 0:
        plot.ylabel("Normalized Dust Density", fontsize = fontsize)
    plot.title(r"$m_p = " + str(mass[number]) + r" $ $M_J$, $\nu_{disk} = 10^{" + str(viscosity[number]) + r"}$", fontsize = fontsize + 1, y = 1.01)

    if number == 2:
        plot.legend(loc = "upper right", bbox_to_anchor = [1.35, 1])

def make_plot():
    # Save + Show
    plot.savefig("triple_dust_comparison.png", bbox_inches = "tight", dpi = my_dpi)
    plot.show()


##### Parameters #####

# Collections
aspect_ratios = []
densities = []
r_ratios = []
rs = []
drs = []
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
densities.append((ten_density, thousand_density))
r_ratios.append((ten_radius_over_dr, thousand_radius_over_dr))
rs.append((ten_radius, thousand_radius))
drs.append((ten_dr, thousand_dr))
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
densities.append((ten_density, four_thousand_density))
r_ratios.append((ten_radius_over_dr, four_thousand_radius_over_dr))
rs.append((ten_radius, four_thousand_radius))
drs.append((ten_dr, four_thousand_dr))
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
densities.append((ten_density, thousand_density))
r_ratios.append((ten_radius_over_dr, thousand_radius_over_dr))
rs.append((ten_radius, thousand_radius))
drs.append((ten_dr, thousand_dr))
angles.append((ten_angle, thousand_angle))

##### Set up plot #####

fig = plot.figure(figsize = (15, 5))
# Remove hspace between plots
plot.subplots_adjust(wspace = 0, hspace = 0)

for i, (aspect_ratios_i, densities_i, r_ratios_i, rs_i, drs_i, angles_i) in enumerate(zip(aspect_ratios, densities, r_ratios, rs, drs, angles)):
    add_to_plot(i, aspect_ratios_i, densities_i, r_ratios_i, rs_i, drs_i, angles_i)

make_plot()
plot.close(fig)

