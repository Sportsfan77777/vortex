"""
Usage: trackVortexCenters.py

Plots vortex centers over time.
"""

directories = ["cm", "hcm", "mm", "hmm", "hum", "um"]
sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001])

### Task Functions ###

def retrieve_density(frame, sizes):
    """ Step 0: Retrieve density """
    density = np.zeros((num_rad, num_theta, len(sizes)))

    for i, size_i in enumerate(sizes):
        fn_i = "../%s-size/gasddens%d.dat" % (size_i, frame)
        density[:, :, i] = fromfile(fn_i).reshape(num_rad, num_theta)

    return density

def method1():
	""" argmax """
	pass

def method2():
	""" center of threshold """
	pass


###############################################################################

##### PLOTTING #####

def make_plot(frame, show = False):

	### Line Plot ###

	plot.legend()