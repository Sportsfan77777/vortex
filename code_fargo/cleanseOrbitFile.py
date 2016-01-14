"""
gets rid of any line in the orbit file that is not 6 entries

http://stackoverflow.com/questions/4710067/deleting-a-specific-line-in-a-file-python
"""

fn = "orbit0.dat"

f = open(fn, "r+")
lines = f.readlines()
f.seek(0)
for line in lines:
    line_sp = line.split()
    if len(line_sp) == 6:
        # There should be six entries in each line (http://fargo.in2p3.fr/Output)
        f.write(line)
f.truncate()
f.close()