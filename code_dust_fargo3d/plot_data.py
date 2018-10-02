from pylab import *
import sys
sys.path.insert(0, '../../utils/python')
from reader import Fields
from advanced import Parameters

n = 1
outputdir = "../../outputs/fargo_dusty" 
fluids = ["gas","dust1","dust2","dust3"]
field  = "dens"
fig = figure(figsize=(20,5))


p = Parameters(outputdir)
for i,fluid in enumerate(fluids):
    data = Fields("../../outputs/fargo_dusty", fluid, n).get_field(field).reshape(p.ny,p.nx)
    ax = fig.add_subplot(1,len(fluids),i+1)
    ax.imshow(log10(data),origin='lower',aspect='auto')

show()
