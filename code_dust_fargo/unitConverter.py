import numpy as np
import matplotlib.pyplot as plot

# some constants
hh = 6.262e-27
kk = 1.381e-16
cc = 2.9979e10

# a frequency "grid"
nu = 1e9 * np.array([30, 100])


# choose a beam size (this is the FWHM in milliarcseconds)
beam = 20.

# specify a temperature
T = 60.   # e.g., appropriate for r = 10 AU

# specify an optical depth at each frequency
tau = 0.1 * (nu/100e9)

# calculate the beam area in steradians
dOmega = (np.pi * (beam*0.001)**2 / (4.*np.log(2.))) * (np.pi/(180.*3600.))**2

# crude estimate of the surface brightness (in microJy/beam)
Bnu = (2*hh*nu**3/cc**2)/(np.exp(hh*nu/(kk*T))-1.)
Snu = 1e29 * Bnu * (1.-np.exp(-tau)) * dOmega