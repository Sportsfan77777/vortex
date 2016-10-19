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

#define MIN_DENSITY 1e-8 // probably neglect for 2-d

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
  double scrh; // 2-d or 3-d related thing

  //printf("");
  //fflush(stdout);

  // Initialize Disk

  #if EOS == IDEAL
   g_gamma = 1.01;
  #endif

  #if ROTATING_FRAME == YES
   g_OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
   g_OmegaZ *= 2.0*CONST_PI;
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

  /// above this, don't touch anything

  /// set up disk parameters here                                                                                                                                       -----------------------------------------------------
  // to use own parameters: g_inputParam
  
  H      = 0.06*R; // R is a vector, so this works (this is consistent with eos.c IFF h/r = const.)
  OmegaK = 2.0*CONST_PI/(R*sqrt(R));
  cs     = H*OmegaK;
  
  scrh   = (0.5*CONST_PI - th)*r/H; // th = theta = z-direction angle (not phi)
  us[RHO] = 1.0/(R)*exp(-0.5*scrh*scrh); /// change radial density profile here
  us[VX1] = us[VX2] = us[VX3] = 0.0;

  us[iVPHI] = R*(OmegaK - g_OmegaZ);
  #if EOS == IDEAL
   us[PRS] = us[RHO]*cs*cs;
  #elif EOS == ISOTHERMAL
//   g_isoSoundSpeed = cs;
   g_isoSoundSpeed = CONST_PI*0.1;
  #endif

}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 * Add torque here (parallelize it)
 * set analysis rate in "pluto.ini"
 *
 * To parallelize it, DOM_LOOP (domain loop) takes care of individual torque cell calculations
 * To sum them up, (each computer has its own torque_cell variable. Call MPI_Allreduce with MPI_SUM to sum them up)
 *  ---- Parallel data reduction ---- 
 * #ifdef PARALLEL
 * MPI_Allreduce (&Ekin, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 * Ekin = scrh;
 * MPI_Allreduce (&Eth_max, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
 * Eth_max = scrh;
 * MPI_Barrier (MPI_COMM_WORLD);
 * #endif
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
  double *x1, *x2, *x3, R, OmegaK, v[256];
  static int do_once = 1;
  
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
  //////**  **//////
  // Inner [Reflective]
  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      // k = z, j = phi, i = r (i iterates from 0 to # ghost radial cells)
      // iterate through ghost zone (2 cells w/ linear interpolation)
      for (nv = 0; nv < NVAR; nv++){
        // copy from active to ghost (reflection)
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - (i + 1)]; // IBEG = zeroth index of active zone
      }
      d->Vc[VX1][k][j][i] *= -1.0; // reflect v_r (VX_1 = index of v_r in this grid = iVR should also work!!!)
      #if GEOMETRY == POLAR
       R = x1[i];
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ); // re-set v_phi to Keplerian
    }
  }
  // Outer [Open (aka Outflow, but can be inflow unless you change it)]
  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      // k = z, j = phi, i = r
      // iterate through ghost zone (2 cells w/ linear interpolation)
      for (nv = 0; nv < NVAR; nv++){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND]; // copy from end of active zone to ghost (open)
      }
      #if GEOMETRY == POLAR
       R = x1[i];
//       d->Vc[iVR][k][j][i] = 0.0;
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
       d->Vc[iVR][k][j][i]  = 0.0;
       d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ); // re-set v_phi to Keplerian
    }
  }
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
        
  mc  = 0.5*g_inputParam[Mdisk]*CONST_Msun;
  mc /= UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;

  // Domain_Loop iterates only through active zone (the real 2-d grid)
  DOM_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] *= mc;
    #if EOS == IDEAL
     d->Vc[PRS][k][j][i] *= mc; // PRS = pressure
    #endif
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 * This is for extra forces.
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
 * Implement the presence of a planet here (by defining its potential)
 * 
 * Usage:
 * This function is called at each time step to find out where the planet is
 * (through its potential)
 *
 * To implement migration,
 * write separate function to calculate torque (in parallel) and da/dt (in parallel)
 * Apply net (da/dt * dt) here.
 * 
 * In practice, it would be better to copy FARGO's RK4 integrator here 
 * (keep variable time steps in mind --- use Saul's book)
 *
 *************************************************************************** */
{
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi;
  double taper_orbits, taper_time, taper;

  #if GEOMETRY == POLAR
   R  = x1; // radial distance from grid cell to star
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
             planet position (initial) (migration not implemented)
   --------------------------------------------- */

  #if ROTATING_FRAME == NO
   double OmegaZ;
   t = g_time; // You can do Mass Tapering here --- use the time variable and the Mplanet input
   if (g_stepNumber == 2) t += g_dt;
   OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
   OmegaZ *= 2.0*CONST_PI;

   xp = cos(OmegaZ*t);
   yp = sin(OmegaZ*t);
  #else
   xp = -1.0;  /* set initial planet position here */
   yp = 0.0;  // Change this so that planet is at phi = pi, not 45 degrees
  #endif

  // The planet potential part is here
  t = g_time;

  taper_orbits = 100;
  taper_time = (2.0 * CONST_PI) * taper_orbits;

  taper = 1.0;
  if (t < taper_time) taper = t / taper_time;

  d = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z); // distance from cell to planet
  rsm = 0.03*R; // smoothing radius --- change this so that's it's a function of planet mass

  if (d > rsm) phiplanet = taper * g_inputParam[Mplanet]/d; // use regular potential
  else phiplanet = taper * g_inputParam[Mplanet]/d*(pow(d/rsm,4.)-2.*pow(d/rsm,3.)+2.*d/rsm); // soften it somehow
  
  // The star potential part is here (along with the planet potential)
  // phi = potential, not angle
  phi  = - 4.0*CONST_PI*CONST_PI/g_inputParam[Mstar];
  phi *= (g_inputParam[Mstar]/r + phiplanet*CONST_Mearth/CONST_Msun);

  return phi; // the total potential
}
#endif
