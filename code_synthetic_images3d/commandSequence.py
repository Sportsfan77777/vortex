'''
Full sequence of commands for synthetic images
'''

import sys
import os
import subprocess

id_number = 0
frames = "500"

step1 = "python generateSingleGrainSyntheticInput.py %s --id %d -o" % (frames, id_number)
os.system(step1)

step2 = "python makeSyntheticImages.py %s --id %d" % (frames, id_number)
os.system(step2)

step3 = "python convolveIntensityMaps.py %s -b 0.2 --id %d" % (frames, id_number)
os.chdir("lambda0870")
os.system(step3)

step4 = "python plotCartesianConvolvedIntensityMaps.py %s --id %d" % (frames, id_number)
os.chdir("beam004")
os.system(step4)