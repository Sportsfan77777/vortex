#include "mp.h"

void update_sgvelocity (VRad, VTheta, Dt)
     PolarGrid *VRad, *VTheta;
     real Dt;
{
  extern boolean SGZeroMode;
  int i, j, l, nr;
  int jm1, lm1;
  real *vrad, *vtheta;
  vrad = VRad->Field;
  vtheta = VTheta->Field;
  nr = VTheta->Nrad;
  
  /* Here we update velocity fields to take into account
     self-gravity */
  for ( i = 0 ; i < nr; i++ ) {
    for (j = 0; j < NSEC; j++) {
      l = i*NSEC + j;
      /* We compute VRAD - half-centered in azimuth - from
	 centered-in-cell radial sg acceleration. */
      if ( i > 0 ) {
	if ( !SGZeroMode )
	  vrad[l] += Dt*( (Rinf[i] - Rmed[i-1]) * SG_Accr[l] +		\
			  (Rmed[i] - Rinf[i]) * SG_Accr[l-NSEC] ) * InvDiffRmed[i];
	if ( SGZeroMode )
	  vrad[l] += Dt*( (Rinf[i] - Rmed[i-1]) * SG_Accr[i+IMIN] +	\
			  (Rmed[i] - Rinf[i]) * SG_Accr[i+IMIN-1] ) * InvDiffRmed[i];
      }
      /* We compute VTHETA - half-centered in radius - from
	 centered-in-cell azimutal sg acceleration. */
      if ( !SGZeroMode ) {
	if (j == 0) 
	  jm1 = NSEC-1;
	else
	  jm1 = j-1;
	lm1 = i*NSEC + jm1;
	vtheta[l] += 0.5 * Dt * ( SG_Acct[l] + SG_Acct[lm1] );
      }
    }
  }
}
