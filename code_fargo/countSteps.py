"""
count number of timesteps in overrun fargo simulations
(should be Ninterm, but it's not)
"""

import numpy as np


count_zero = 0
count_one = 0
count_two = 0

data = np.loadtxt("bigplanet0.dat")
timesteps = data[:,0]

for t in timesteps:
	if t == 0:
		count_zero += 1
	elif t == 1:
		count_one += 1
	elif t == 2:
		count_two += 1
	else:
		break

print count_zero, count_one, count_two


