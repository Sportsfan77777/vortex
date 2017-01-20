"""
writes a title to a directory (used in plots)

Usage:
python writeTitle.py "title"

title is stored in 'title.p'
"""

import sys
import pickle

title = sys.argv[1]

fn = "title.p"
pickle.dump(title, open(fn, 'wb'))