"""
plots sin waves with e^() bounds
"""


import numpy as np 
from matplotlib import pyplot as plot




def exp(r, x):
    return np.exp(r * x)

rate = 0.1

decay = lambda x : exp(-rate, x)
flat = lambda x : exp(0, x)
grow = lambda x : exp(rate, x)

xs = np.linspace(0, 10, 100)

# Plotting
fontsize = 14
linewidth = 3

def make_plot(choice):
    """ Three options: growing, decaying, flat """

    choice_f = flat


    if choice is "growing":
        choice_f = grow
    elif choice is "decaying":
        choice_f = decay

    ys = [choice_f(x) * np.sin(np.pi * x) for x in xs]
    upper_bounds = np.array([choice_f(x) for x in xs])
    lower_bounds = upper_bounds * -1

    plot.plot(xs, ys, 'b', linewidth = linewidth + 1)
    plot.plot(xs, upper_bounds, 'r', linewidth = linewidth)
    plot.plot(xs, lower_bounds, 'r', linewidth = linewidth)

    # Limit
    plot.ylim(-3, 3)

    # Annotate
    plot.xlabel("Time", fontsize = fontsize)
    plot.ylabel(r"$\delta \Sigma$", fontsize = fontsize + 5)
    plot.title("", fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s_mode.png" % choice)
    plot.show()
    plot.cla()


make_plot("growing")
make_plot("decaying")
make_plot("flat")


