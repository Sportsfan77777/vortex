#include "mp.h"

void compute_selfgravity (Rho, RadialVelocity, AzimutalVelocity, DeltaT, SGUpdate)
     PolarGrid *Rho;
     PolarGrid *RadialVelocity;
     PolarGrid *AzimutalVelocity;
     real DeltaT;
     boolean SGUpdate;
{
  extern boolean SGZeroMode;
  if ( SG_initcounter == 0 ) {
    init_sg ();
    if ( !SGZeroMode )
      compute_fftkernel ();
  }
  /* Only mode m=0 of disk self-gravity is taken into account */
  if ( SGZeroMode )
    compute_SGZeroMode (Rho);
  /* Complete disk self-gravity treatment in this case */
  if ( !SGZeroMode ) {
    compute_fftdensity (Rho);
    /* Here we compute radial and azimutal components of sg acceleration
       as a convolution product of reduced density and kernel arrays */
    compute_sgacc (Rho);
  }
  if ( SGUpdate ) {
    /* Computes polar components of acceleration and
       updates values of vrad, vtheta at each step */
    update_sgvelocity (RadialVelocity, AzimutalVelocity, DeltaT);
  }
  SG_initcounter++;
}
