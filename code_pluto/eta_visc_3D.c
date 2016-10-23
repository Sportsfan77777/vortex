#include "pluto.h"
#include <stdio.h>

/*---------------------------------------------------------------------------*/
/*---- Specification of explicit first and second viscosity coefficients ----*/
/*---------------------------------------------------------------------------*/

void eta_visc_func(real *v, real x1, real x2, real x3, 
                   real *eta1_visc, real *eta2_visc )
{/*assume sphericla polars*/
  *eta1_visc = v[RHO]*kinematic_visc(x1, x2); /*visc*fz;*/
//  *eta1_visc = kinematic_visc(x1, x2); //kinematic visc prop to 1/rho

//rough model to account for viscosity reduction when density increases
/*assume spherical polars*/
/*  double z, fz, R, zeta;
  double visc; 

  R = x1*sin(x2);
  z = x1*cos(x2);
  zeta = (z/bigH(R))/sqrt(2.0);

//normal, smooth viscosity
  visc = 0.0;
  if(g_inputParam[nu] > 0.0)       visc  = g_inputParam[nu]*g_unitLength*g_unitLength*omega_k(g_unitLength);
  if(g_inputParam[nu_alpha] > 0.0) visc  = g_inputParam[nu_alpha]*sqrt(csq(R))*bigH(R);

  fz = v[RHO]*sqrt(2.0*CONST_PI)*bigH(R)*exp(zeta*zeta)*erfc(zeta);
  fz/= surface_density(g_unitLength)*erfc(g_inputParam[visc_jump_height]/sqrt(2.0)); 
  fz = exp(-fz);

  *eta1_visc = v[RHO]*visc*(1.0 + (g_inputParam[visc_jump_amp] - 1.0)*fz); 
*/
  *eta2_visc = 0.0;
}
