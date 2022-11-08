'''
Full sequence of commands for synthetic images
'''

import sys
import os
import subprocess

id_number = 0
cores = 16
frames = "500 3000 50"

step1 = "python generateSingleGrainSyntheticInput.py %s -c %d --id %d" % (frames, cores, id_number)
os.system(step1)

step2 = "python makeSyntheticImages.py %s -c %d --id %d" % (frames, cores, id_number)
os.system(step2)

step3 = "python convolveIntensityMaps.py %s -c %d -b 0.2 --id %d" % (frames, cores, id_number)
os.chdir("lambda0870")
os.system(step3)

step4 = "python plotCartesianConvolvedIntensityMaps.py %s -c %d --id %d" % (frames, cores, id_number)
os.chdir("beam004")
os.system(step4)