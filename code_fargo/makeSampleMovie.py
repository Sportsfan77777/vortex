"""
makes a movie out of the sample selection in a particular directory

It is assumed that the sample is small

python makeSampleMovie directory_name
"""

import sys
import os
import pickle
import subprocess
import shutil

# Get Sample

directory = sys.argv[1]
sample = pickle.load("%s/sample.p" % directory)

##### Make Movie #####

# cd directory
os.chdir(directory)

# Get all files
files = glob.glob("*.png" % directory)

# Create a set of all of the prefixes in prefix%d.png
prefixes = set()
for filename in files:
	png_name = (filename.split("."))[0]
	# Get rid of number part
	prefix = "".join(char for char in png_name if not char.isdigit())
	prefixes.add(prefix)

# For each prefix, make a movie
for prefix in prefixes:
	# Create tmp files numbered 0 to num_samples
	for i, frame in enumerate(sample):
		old = "%s%d.png" % (prefix, frame)
		new = "tmp4movie_%s%d.png" % (prefix, i) # re-numbered from frame to i

	# Movie Command
	fps = 1
	path = ("tmp4movie_%s" % prefix) + "%d.png"
	output = "sample_%s.mp4" % prefix

	command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)


# Delete all tmp files
tmp_files = glob.glob("tmp4movie*")
for tmp_file in tmp_files:
	os.remove(tmp_file)

