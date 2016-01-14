"""
gets rid of any line in the torque file that is not 10 entries

http://stackoverflow.com/questions/4710067/deleting-a-specific-line-in-a-file-python
"""

fn = "tqwk0.dat"

f = open(fn, "r+")
lines = f.readlines()
f.seek(0)
for line in lines:
	line_sp = line.split()
    if len(line_sp) == 10:
    	# There should be ten entries in each line (http://fargo.in2p3.fr/Output)
        f.write(line)
f.truncate()
f.close()