"""
adds entry (or changes entry) in a pickle file

requires pickle stores a dictionary!

Usage:
python addToPickle.py storage.p Thing 55 # Adds Thing to storage.p with a value of 55
"""

import pickle
import sys

fn = sys.argv[1]
entry = sys.argv[2]
value = sys.argv[3]

dictionary = pickle.load(open(fn, 'rb'))
dictionary[entry] = value
pickle.dump(dictionary, open(fn, 'wb'))