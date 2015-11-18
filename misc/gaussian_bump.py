"""
plot gaussian bumps

arg = width of gaussian
"""

import sys
import numpy as np 
from matplotlib import pyplot as plot


sigma = int(sys.argv[1])

def gaussian(x, sigma = sigma):
	#const = (2 * np.pi * sigma**2)**(-0.5)
	const = 1
	return const * np.exp(-x**2 / sigma**2)


xs = np.linspace(-10, 10, 1000)
ys = np.array([gaussian(x) for x in xs])

plot.xlim(-10, 10)
plot.ylim(0, 1)

plot.plot(xs, ys, 'b', linewidth = 4)
plot.savefig("bump%d.png" % sigma)
plot.show()