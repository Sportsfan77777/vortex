"""
turn matplotlib.use('Agg') on

Usage: 
python switchAgg.py on fn
python switchAgg.py off fn
"""

import sys

off = "#matplotlib.use('Agg')"
on = "matplotlib.use('Agg')"

def find_replace(old, new, fn):
    # Read
    with open(fn, 'r') as f:
       txt = f.read()

    # Replace
    txt = txt.replace(old, new)

    # Write
    with open(fn, 'w') as f:
       f.write(txt)

######## Find and Replace! ########
switch = sys.argv[1]
fn = sys.argv[2] # Argument is fn

if switch is "on":
    # Turn on
    find_replace(off, on, fn)
elif switch is "off":
    # Turn off
    find_replace(on, off, fn)

