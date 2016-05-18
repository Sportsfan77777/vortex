"""
choose when the vortex dies

python chooseVortexEnd.py
"""

import pickle as p
import numpy as np

end_fn = "end.p"

end_candidates_fn = "end_candidates.p"
end_candidates = p.load(open(end_candidates_fn, "rb"))

if len(end_candidates) == 1:
	end = (end_candidates[0])[-1]
	p.dump(end, open(end_fn, "wb"))
else:
	print end_candidates
	end = int(input("Choose the end! "))
	p.dump(end, open(end_fn, "wb"))