
#include "fargo3d.h"

void rescale () {
DT *= sqrt(R0*R0*R0/G/MSTAR);
MASSTAPER *= sqrt(R0*R0*R0/G/MSTAR);
SIGMA0 *= MSTAR/(R0*R0);
OMEGAFRAME *= sqrt(G*MSTAR/(R0*R0*R0));
}
