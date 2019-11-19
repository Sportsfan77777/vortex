import sys
import numpy as np

data = np.loadtxt("bigplanet0.dat")
data = data[-2:-1] # second-to-last line, in case it is in the middle of writing

frame = data[0]
print frame