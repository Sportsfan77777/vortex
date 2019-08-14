#include "mp.h"

void init_planetarysys_withSG (sys)
     PlanetarySystem *sys;
{
  extern boolean SGZeroMode;
  int k, ipl;
  real r, dist, ri, rip1, dr, sgacc;
  if ( !SGZeroMode )
    mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
  else
    GLOBAL_AxiSGAccr = SG_Accr;
  /* Planetary system initialization in self-gravity cases: 
     planets are put in a fixed circular orbit, we need to know 
     radial sg acceleration felt by planets. */
  for ( k = 0; k < sys->nb; k++ ) {
    r = sys->x[k];
    /* dist denotes the planet's semi-major axis */
    dist = (real)(r / (1. + ECCENTRICITY));
    ipl = 0;
    while (GlobalRmed[ipl] <= dist) ipl++;
    ri = GlobalRmed[ipl];
    rip1 = GlobalRmed[ipl+1];
    dr = rip1 - ri;
    sgacc = (dist - ri)*GLOBAL_AxiSGAccr[ipl+1] + (rip1 - dist)*GLOBAL_AxiSGAccr[ipl];
    sgacc /= dr;
    /* sgacc is the radial sg acc. at the planet's semi-major axis */
    sys->vy[k] *= (real)sqrt (1. - dist*dist*sgacc);
  }
}
