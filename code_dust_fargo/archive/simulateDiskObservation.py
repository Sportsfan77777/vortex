"""
Simulates Disk Observations

Taken from:
https://casaguides.nrao.edu/index.php/Protoplanetary_Disk_Simulation_(CASA_4.4)
"""

# In CASA
default("simobserve")

# This reports image header parameters in the Log Messages window
imhead("ppdisk672_GHz_50pc.fits")

# In CASA
ia.open("ppdisk672_GHz_50pc.fits")

# In CASA
axesLength = ia.shape()
# Divide the first two elements of axesLength by 2.
center_pixel = [ x / 2.0 for x in axesLength[:2] ]
# Feed center_pixel to ia.toworld and and save the RA and Dec to ra_radians and dec_radians
(ra_radians, dec_radians) = ia.toworld( center_pixel )['numeric'][:2]
ia.close()

ra_hms  = qa.formxxx(str(ra_radians)+"rad",format='hms',prec=5)
dec_dms = qa.formxxx(str(dec_radians)+"rad",format='dms',prec=5)

project = "psim2"
skymodel = "ppdisk672_GHz_50pc.fits"

setpointings       =  True
direction          =  "J2000 18h00m00.031s -22d59m59.6s"
mapsize            =  "0.76arcsec"

obsmode            =  "int"
totaltime          =  "1200s"

antennalist        =  "alma.out20.cfg"

thermalnoise = ''

simobserve()

default ("simanalyze")
project = "psim2"
image = True

modelimage = "ppdisk672_GHz_50pc.fits"
vis = project + ".alma.out20.ms"
imsize = [192, 192]

niter = 10000
threshold = "1e-7Jy"
weighting = "natural"

analyze = True  
showuv = False
showresidual = True  
showconvolved = True

graphics = "both"
verbose = True
overwrite = True

simanalyze()