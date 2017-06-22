/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet interaction problem.

  Simulate the interaction of a planet embedded in a disk as described
  in section 3.4 of Mignone et al., A&A (2012) 545, A152.
  This test is a nice benchmark for the FARGO module and the
  \c ROTATING_FRAME switch.
  For testing-purposes  no viscosity is used here.
  The initial condition consists of a locally isothermal configuration
  with temperature profile \f$\propto T^{-1}\f$  yielding a disk vertical
  height to radius of \c 0.05.
  The gravitational potential due to the presence of the star and the
  planet is defined in BodyForcePotential() function.

  The conventions used throught the implementation are the following:
 
  - \c r  = spherical radius
  - \c R  = cylindrical radius
  - \c z  = cylindrical height
  - \c th = meridional angle
  
  The test can be carried out in polar (2D or 3D) or spherical (3D)
  coordinates and the following parameters determine the initial configuration:
   
  -# <tt>g_inputParam[Mstar]</tt>: controls the star mass (in solar masses)
  -# <tt>g_inputParam[Mdisk]</tt>: controls the disk mass (in solar masses)
  -# <tt>g_inputParam[Mplanet]</tt>: sets the planet mass (in earth masses)
  -# <tt>g_inputParam[Viscosity]</tt>: sets the amount of viscosity


  Computation can be carried in the rotating or in the observer's frame
  of reference (\c ROTATING_FRAME to \c YES or \c NO, respectively).
  In particular:

  - Configurations #01 and #02 are in 2D polar coordinates without and with
    FARGO, in the rotating frame.
  - Configurations #03, #04 and #05 are in spherical 3D coordinates with
    and without FARGO in the rotating frame
  - Configurations #06 and #07 are in 2D polar coordinates but in the
    observer's frame
  - Configuration #08 employs static AMR  (grid levels are spatial
    dependent but not dependent on time) in the rotating frame.
    

  \image html hd_disk_planet.08.png "Density map for configuration #08 using AMR" width=1cm

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012

  \b References:
     - "A Conservative orbital advection scheme for simulations
        of magnetized shear flows with the PLUTO Code"
        Mignone et al, A&A (2012)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MIN_DENSITY 1e-8

static void NormalizeDensity (const Data *d, Grid *g);
#if ROTATING_FRAME == NO
 #define g_OmegaZ  0.0
#endif

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, cs;

  #if EOS == IDEAL
   g_gamma = 1.01;
  #endif
  
  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
  #elif GEOMETRY == SPHERICAL
   r  = x1;
   th = x2;
   R  = r*sin(th);
   z  = r*cos(th);
  #endif
  
  // Initialize Field Variables
  us[RHO] = density3D(R, z);

  us[VX1] = radialVelocity_rComponent(R, th, z); // r 
  us[VX2] = radialVelocity_thetaComponent(R, th, z); // theta (out of the plane)
  us[VX3] = azimuthalVelocity3D(R, z); // phi (azimuthal)

  #if EOS == IDEAL
   us[PRS] = pressure(R, z);
  #elif EOS == ISOTHERMAL
//   g_isoSoundSpeed = cs;
     g_isoSoundSpeed = CONST_PI*0.1; // keep in mind, this is for Omega = 2 * pi
  #endif

#if DUST == YES
  us[RHO_D] = 1.e-4;
  us[VX1_D] = us[VX2_D] = us[VX3_D] = 0.0;
  us[VX2_D] = us[iVPHI];
#endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3, R, z, OmegaK, v[256];
  static int do_once = 0;
  
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  #if DIMENSIONS == 3
  if (side == 0){
    if (do_once){
      NormalizeDensity(d, grid);
      do_once = 0;
    }
  }
  #endif

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2;
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[RHO][k][j][i]   = density3D(R, z);
      d->Vc[VX1][k][j][i]   = radialVelocity_rComponent(R, th, z);
      #if GEOMETRY == CYLINDRICAL
          d->Vc[VX2][k][j][i]   = 0.0; // vz or vtheta
      #elif GEOMETRY == SPHERICAL
          d->Vc[VX2][k][j][i]   = radialVelocity_thetaComponent(R, th, z); // vtheta
      #elif GEOMETRY == POLAR && DIMENSIONS == 3
          d->Vc[VX3][k][j][i]   = 0.0; // vz
      #endif
      d->Vc[iVPHI][k][j][i] = azimuthalVelocity3D(R, z);
      #if EOS == IDEAL
          d->Vc[PRS][k][j][i]  = pressure(R, z);
      #endif
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2;
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[RHO][k][j][i]   = density3D(R, z);
      d->Vc[VX1][k][j][i]   = radialVelocity_rComponent(R, th, z);
      #if GEOMETRY == CYLINDRICAL
          d->Vc[VX2][k][j][i]   = 0.0; // vphi
      #elif GEOMETRY == SPHERICAL
          d->Vc[VX2][k][j][i]   = radialVelocity_thetaComponent(R, th, z); // vtheta
      #elif GEOMETRY == POLAR && DIMENSIONS == 3
          d->Vc[VX3][k][j][i]   = 0.0; // vtheta
      #endif
      d->Vc[iVPHI][k][j][i] = azimuthalVelocity3D(R, z);
      #if EOS == IDEAL
          d->Vc[PRS][k][j][i]  = pressure(R, z);
      #endif
    }
  }

 if (side == X2_BEG){
    X2_BEG_LOOP(k,j,i){
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2;
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[RHO][k][j][i]   = density3D(R, z);
      d->Vc[VX1][k][j][i]   = radialVelocity_rComponent(R, th, z);
      #if GEOMETRY == CYLINDRICAL
          d->Vc[VX2][k][j][i]   = 0.0; // vz or vtheta
      #elif GEOMETRY == SPHERICAL
          d->Vc[VX2][k][j][i]   = radialVelocity_thetaComponent(R, th, z); // vtheta
      #elif GEOMETRY == POLAR && DIMENSIONS == 3
          d->Vc[VX3][k][j][i]   = 0.0; // vz
      #endif
      d->Vc[iVPHI][k][j][i] = azimuthalVelocity3D(R, z);
      #if EOS == IDEAL
          d->Vc[PRS][k][j][i]  = pressure(R, z);
      #endif
    }
  }

/* original below ************************

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] *= -1.0;
      #if GEOMETRY == POLAR
       R = x1[i];
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
#if DUST == YES      
//      NDUST_LOOP(nv) d->Vc[nv][k][j][i] = 0.0;
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
#endif      
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      #if GEOMETRY == POLAR
       R = x1[i];
//       d->Vc[iVR][k][j][i] = 0.0;
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
       d->Vc[iVR][k][j][i]  = 0.0;
       d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
#if DUST == YES      
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
#endif      
    }
  }
  ********************** */
}

/* ************************************************************** */
void NormalizeDensity (const Data *d, Grid *grid)
/*
 *
 * Normalize density and pressure as   rho -> K*rho, where
 *
 *   K = M/(\sum rho*dV)
 *
 **************************************************************** */
{
  int   i, j, k;
  double mc;
        
  //mc  = 0.5*g_inputParam[Mdisk]*CONST_Msun;
  mc = 0.1;
  mc /= UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
  DOM_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] *= mc;
    #if EOS == IDEAL
     d->Vc[PRS][k][j][i] *= mc;
    #endif
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double R, r, z, th, x, y;
  double angle_cell, angle_planet, angular_separation;
  double xp, yp, t, phi;

  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
   x  = R*cos(x2);
   y  = R*sin(x2);
  #elif (GEOMETRY == SPHERICAL)
   r  = x1;
   th = x2;
   R = r*sin(th);
   z = r*cos(th);
   x = r*sin(th)*cos(x3);
   y = r*sin(th)*sin(x3);
  #endif

/* ---------------------------------------------
             planet position
   --------------------------------------------- */

  #if ROTATING_FRAME == NO
   double OmegaZ;
   t = g_time;
   if (g_stepNumber == 2) t += g_dt;
   OmegaZ  = sqrt(1.0 + g_inputParam[P_Mplanet]/g_inputParam[P_Mstar]*CONST_Mearth/CONST_Msun);
   //OmegaZ *= 2.0*CONST_PI;

   xp = cos(OmegaZ*t);
   yp = sin(OmegaZ*t);
  #else
   xp = 1.0/sqrt(2.0);  /* initial planet position */
   yp = 1.0/sqrt(2.0); 
  #endif

  angle_cell = x3; // phi (angle in the plane)
  angle_planet = atan2(yp, xp); 
  angular_separation = angle_cell - angle_planet; // between cell and planet

  phi = stellarPotential(r) + planetPotential(r, R, angular_separation);

  return phi;
}
#endif
