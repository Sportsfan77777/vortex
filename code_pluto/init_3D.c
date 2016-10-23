/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet problem.

  Disk-Planet setup.
 
  Reference paper:
   "A Conservative orbital advection scheme for simulations
    of magnetized shear flows with the PLUTO Code"
    Mignone et al, A&A (2012)
 
  -------------------------------------------------------------
   Independent of the coordinate system (polar/spherical), we
   adopt the following conventions:
 
    r = spherical radius
    R = cylindrical radius
    z = cylindrical height
   th = meridional angle

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "float.h"
#define MIN_DENSITY 0.0
#define g_OmegaZ  0.0
#define bigG 1.0
#define mstar 1.0
#define MSEED 6

static void NormalizeDensity (const Data *d, Grid *g);
static double ran2(long int *idum);

/*MKL: accretion disk functions*/
static double grav_pot3_cylin(const double bigR, const double z);
static double bump(const double bigR);
static double dlogbump(const double bigR);
static double d2logbump(const double bigR);
static double dlogomega_iso(const double bigR);
static double dlog_rho0(const double bigR);
static double azivel(const double bigR, const double z);
static double radvel(const double x1, const double x2);
static double polvel(const double x1, const double x2);
static double inner_hole(const double bigR);

int    pert; 
double sig0, bwidth, dampin, dampout, planeton, switchon, planeton_switchon;
double H0;
/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int    m;
  double r, th, R, z, H, OmegaK, cs, totmodes, delta;
  double scrh;

  #if EOS == IDEAL
   g_gamma = g_inputParam[gmma];
  #endif

  g_unitLength   = g_inputParam[r0];
  g_unitDensity  = density3D(g_unitLength, 0.0);
  g_unitVelocity = g_unitLength*omega_k(g_unitLength);

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
  
  sig0      = sqrt(csq(g_inputParam[rout]))*omega_k(g_inputParam[rout]);
  sig0     /= g_inputParam[qout]*CONST_PI*bigG*pow(g_inputParam[rout]/g_unitLength, -g_inputParam[smallp]);//*inner_hole(g_inputParam[rout]);
  sig0     /= bump(g_inputParam[rout]);

  H0        = bigH(g_unitLength);

  bwidth = g_inputParam[bump_width]*g_unitLength;
  dampin = g_inputParam[damp_in]*g_inputParam[rmin];
  dampout= g_inputParam[damp_out]*g_inputParam[rout];
  planeton= g_inputParam[planet_on]*2.0*CONST_PI;
  switchon= g_inputParam[switch_on]*2.0*CONST_PI;
  planeton_switchon = planeton + switchon;  

  us[RHO] = density3D(R, z);
  us[VX1] = us[VX2] = us[VX3] = 0.0;

  pert = 0;  

//background vel field (include viscous accretion)
  us[VX1]   = radvel(x1, x2);
  us[iVPHI] = azivel(R, z);
  us[VX2]   = polvel(x1, x2);

  #if EOS == IDEAL
   us[PRS] = us[RHO]*csq(R);
  #elif EOS == ISOTHERMAL
   g_isoSoundSpeed = sqrt(csq(g_unitLength));
  #endif

}

static double grav_pot3_cylin(const double bigR, const double z)
{
  double d;
  double star_pot;
  
  d =sqrt(bigR*bigR + z*z);

  star_pot  = -bigG*mstar/d;
  return star_pot;
}

double omega_k(const double bigR)
{
  double omega_sq;

  omega_sq = bigG*mstar/pow(bigR, 3.0);
  return sqrt(omega_sq);
}

double csq(const double bigR)
{
  double soundspeed_sq;

  soundspeed_sq = pow(g_inputParam[smallh], 2.0);
  soundspeed_sq*= bigG*mstar/g_unitLength; 
  soundspeed_sq*= pow(bigR/g_unitLength, -g_inputParam[smallq]); 

  return soundspeed_sq;
}

double bigH(const double bigR)
{
  return sqrt(csq(bigR))/omega_k(bigR);
}

double bump(const double bigR)
{  
  return 1.0 + (g_inputParam[amp] - 1.0)*exp(-0.5*pow((bigR - g_unitLength)/bwidth, 2.0)); 
}

static double dlogbump(const double bigR)
{
  double dbump;

  dbump   = (bump(bigR) - 1.0);
  dbump  *=-(bigR - g_unitLength)/pow(bwidth, 2.0);
  dbump  /= bump(bigR);

  return dbump;
}

static double d2logbump(const double bigR)
{
  double DR2, d2bump, B;

  DR2 = pow(bwidth, 2.0);
  B   = bump(bigR);

  d2bump = dlogbump(bigR)/B;
  d2bump*= -(bigR - g_unitLength)/DR2;
  d2bump += -(1.0 - 1.0/B)/DR2;

  return d2bump;
}

static double inner_hole(const double bigR)
{  
  double fac = sin(CONST_PI/2.0 - 4.0*g_inputParam[smallh]);
  double Hin;

  Hin = 4.0*bigH(g_inputParam[rmin]);
  return 1.0 - sqrt(g_inputParam[rmin]*fac/(bigR + Hin));
}


double surface_density(const double bigR)
{
  double sig_s;

  sig_s = sig0*pow(bigR/g_unitLength, -g_inputParam[smallp]);
  return sig_s*bump(bigR);//*inner_hole(bigR);
}

double density3D(const double bigR, const double z)
{
  double vertical, sigma;

  vertical = -(grav_pot3_cylin(bigR, z) - grav_pot3_cylin(bigR,0.0));
  vertical/= csq(bigR); 
 
  sigma = surface_density(bigR);
  sigma/= sqrt(2.0*CONST_PI)*bigH(bigR);
  
  return sigma*exp(vertical);
}

static double dlog_rho0(const double bigR)
{
  double fac = sin(CONST_PI/2.0 - 4.0*g_inputParam[smallh]);
  double dsigma_s, dbump, dbigH, dtaper, dtaper_in;
  double dcs, domega_k, Hin;

  Hin = 4.0*bigH(g_inputParam[rmin]);

  dsigma_s = -g_inputParam[smallp]/bigR;
//  dsigma_s+= 0.5*sqrt(g_inputParam[rmin]*fac)*pow(bigR+Hin, -1.5)/inner_hole(bigR);

  dbump   = (bump(bigR) - 1.0);
  dbump  *=-(bigR - g_unitLength)/pow(bwidth, 2.0);
  dbump  /= bump(bigR);

  dcs     =  -0.5*g_inputParam[smallq]/bigR;
  domega_k = -1.5/bigR;

  dbigH   = dcs - domega_k;
  
  return dsigma_s + dbump - dbigH;
}

static double azivel(const double bigR, const double z)
{
  double dcs2, cs2, dphi0, vsq; 
  
  cs2     = csq(bigR);
  dcs2    =-g_inputParam[smallq]/bigR;
  dphi0   = bigR*pow(omega_k(bigR),2.0);
  
  vsq     = cs2*dlog_rho0(bigR) +
    dcs2*(grav_pot3_cylin(bigR,z)-grav_pot3_cylin(bigR,0.0)+cs2) +
    dphi0; 
  vsq    *= bigR;
  return sqrt(vsq); 
}


static double dlogomega_iso(const double bigR)
{//dln(omega)/dln R for strictly isothermal case, or just return the kep value of -1.5 if not strictly isothermal

  double cs2, omega, d2rho, d2phi;
  double result;

  if(g_inputParam[smallq] > 0.0){
  return -1.5; //keplerian value of shear
  }else{
  cs2   = csq(g_unitLength);
  omega = azivel(bigR, 0.0)/bigR;

  d2rho   = (1.5 + g_inputParam[smallp])/pow(bigR, 2.0) + d2logbump(bigR); //artificially surpress bump
  d2phi   = -2.0*bigG*mstar/pow(bigR, 3.0);

  result = cs2*d2rho + d2phi;
  result/= bigR*omega*omega;
  result-= 1.0/bigR;
  result/= 2.0;

  return result*bigR; //this is dlog(omega)/dlog(R) for strictly isothermal
  }
}

double pep_eos(const double r, const double theta, const double phi)
{
 double const hp=0.5, eosn=3.5;
 double Hp, dp, x, y, z, bigR, t, omega_p, rsm, xp, yp;
 double csmod;

 if(g_inputParam[mplanet] > 0.0) {
 bigR = r*sin(theta);
 x = bigR*cos(phi);
 y = bigR*sin(phi);
 z = r*cos(theta);

 t = g_time;

  if(g_intStage == 2){
#if TIME_STEPPING == RK2
    t += 0.5*g_dt;
#elif TIME_STEPPING == RK3
    t += 0.25*g_dt;
#endif
  }
  if(g_intStage == 3) t += 1.0/3.0*g_dt;

   omega_p = omega_k(g_unitLength);
   xp =  g_unitLength*cos(omega_p*t + CONST_PI);
   yp =  g_unitLength*sin(omega_p*t + CONST_PI);

   rsm = g_unitLength*pow(g_inputParam[mplanet]/mstar/3.0, 1.0/3.0);
   rsm *= g_inputParam[soft];
   dp = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z + rsm*rsm);

   Hp = hp*dp;

   csmod = Hp*sqrt( 1.0 + (g_inputParam[mplanet]/mstar)*pow(bigR/dp, 3.0) );
   csmod/=pow( pow(bigH(bigR), eosn) + pow(Hp, eosn), 1.0/eosn);
   csmod-= 1.0;
   csmod*= turn_on(t);
   csmod+= 1.0;
  } else {
  csmod = 1.0;
  }

   return csmod;
}

static double radvel(const double x1, const double x2)
{//assume spherical polars
  double vr, visc, R;

  R = x1*sin(x2);
  visc = kinematic_visc(x1, x2);
  vr = dlogomega_iso(R)*visc/x1; 
  
  return vr;
}

static double polvel(const double x1, const double x2)
{//assume spherical polars
  double vtheta, visc, R;

  R = x1*sin(x2);
  visc = kinematic_visc(x1, x2);
  vtheta = dlogomega_iso(R)*visc;
  vtheta/= x1*tan(x2);

  return vtheta;
}

double kinematic_visc(const double x1, const double x2)
{
  /*assume spherical polars*/
  double z, fz, R, transition_height, transition_radius;
  double visc, vert_coord, gauss;

  R = x1*sin(x2);
  z = x1*cos(x2);

  /*normal, smooth viscosity*/
  visc = 0.0;
  if(g_inputParam[nu] > 0.0)       visc  = g_inputParam[nu]*g_unitLength*g_unitLength*omega_k(g_unitLength);
  if(g_inputParam[nu_alpha] > 0.0) visc  = g_inputParam[nu_alpha]*sqrt(csq(R))*bigH(R);

//note: steady state disks require sigma*nu = const in 2D, or rho*nu = f(z) in 3D

  if(g_inputParam[visc_jump_H]>0.0){
// vertical jump as a function of normalized height at r0, Z = z/(scaleheight_0)
// set rho*nu = const. here 
//

//    this viscosity is a step:  
//     visc *= (density3D(g_unitLength, z)*dlogomega_iso(g_unitLength))/(density3D(R, z)*dlogomega_iso(R)); 

//     same as above but remove the bump factor 
       visc *= ( density3D(g_unitLength, z)/bump(g_unitLength) )/( density3D(R, z)/bump(R) );

//     visc *= (density3D(g_unitLength, z)*dlogomega_iso(g_unitLength))/(dlogomega_iso(R));//the the above mult by rho
     vert_coord = z/H0;
  } else {
// vertical jump as function of angle from midplane, psi = pi/2 - theta 
// set sigma*nu = const. here 
//    visc *= (surface_density(g_unitLength)*dlogomega_iso(g_unitLength))/(surface_density(R)*dlogomega_iso(R));


//    visc *=  surface_density(g_unitLength)/surface_density(R);

    vert_coord = CONST_PI/2.0 - x2;
  }

//append vertical dependence here 
  transition_height = g_inputParam[visc_jump_height];

  fz = 0.5*(1.0 + tanh((vert_coord - transition_height)/g_inputParam[visc_jump_width]));
  fz+= 0.5*(1.0 - tanh((vert_coord + transition_height)/g_inputParam[visc_jump_width]));

  fz*= g_inputParam[visc_jump_amp] - 1.0;
  fz+= 1.0;

  if(g_inputParam[visc_rtrans] > 0.0){
    double fr;
    fr = 0.5*(1.0 + tanh((R - g_unitLength)/(g_inputParam[bump_width]*bigH(g_unitLength))));
    fz *= fr;
    fz += (1.0 - fr)*g_inputParam[visc_jump_amp];
  }


  return visc*fz;
}

double turn_on(const double time)
{
  double temp;
  
  if(time < planeton){

     return 0.0;
   } else if( (time >= planeton) && (time < planeton_switchon)){
     temp = sin(0.5*CONST_PI*(time - planeton)/switchon);
     temp*= temp;

     return temp;
   } else {

     return 1.0;
   }
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
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv, gjbeg;
  double *x1, *x2, *x3, R, OmegaK, v[256];
  static int do_once = 1;
  double r, z, damp_time, vphi_zero, d_zero, p_zero, vrad_zero, vtheta_zero;
  double taper, theta_min, dt; 
  double *x2r, *x2_glob;  
 
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x; x2r = grid[JDIR].xr; x2_glob = grid[JDIR].xl_glob; 
  x3 = grid[KDIR].x;

  theta_min = x2_glob[grid[JDIR].gbeg];


  if(g_intStage == 1) dt = g_dt;
  if(g_intStage == 2){
#if TIME_STEPPING == RK2 
    dt = 0.5*g_dt;
#elif TIME_STEPPING == RK3
    dt = 0.25*g_dt;
#endif
  }
  if(g_intStage == 3) dt = 2.0/3.0*g_dt;


  if (side == X1_BEG){
    if(g_inputParam[rad_obc] > 0.0){
    X1_BEG_LOOP(k,j,i){
      for (nv = 0; nv < NVAR; nv++){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
      }
      if(d->Vc[VX1][k][j][IBEG] > 0.0 ) d->Vc[VX1][k][j][i] = 0.0;
  
       r = x1[i];
       R = r*sin(x2[j]);
//       z = r*cos(x2[j]);
//       d->Vc[VX3][k][j][i] = azivel(R,z);
    
     d->Vc[VX3][k][j][i] = sqrt(mstar*bigG/R - g_inputParam[smallq]*csq(R) );

    }
   } else {
    X1_BEG_LOOP(k,j,i){
       r = x1[i]; 
       R = r*sin(x2[j]); 
       z = r*cos(x2[j]); 
      
       d_zero = density3D(R, z); 
      
       d->Vc[VX1][k][j][i] = radvel(x1[i], x2[j]); 
       d->Vc[VX2][k][j][i] = polvel(x1[i], x2[j]); 
      
       d->Vc[RHO][k][j][i] = d_zero; 
       d->Vc[VX3][k][j][i] = azivel(R,z); 
      
 #if EOS == IDEAL 
       d->Vc[PRS][k][j][i] = d_zero*csq(R); 
 #endif	 
    }
  }
  }

  if (side == X1_END){

    if(g_inputParam[rad_obc] > 0.0){

    X1_END_LOOP(k,j,i){

       for (nv = 0; nv < NVAR; nv++){  
           d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];  
         }  
      
       if(d->Vc[VX1][k][j][IEND] < 0.0) d->Vc[VX1][k][j][i] = 0.0;  

       r = x1[i];
       R = r*sin(x2[j]);
//       z = r*cos(x2[j]);
//       d->Vc[VX3][k][j][i] = azivel(R,z);
      
      d->Vc[VX3][k][j][i] = sqrt(mstar*bigG/R - g_inputParam[smallq]*csq(R) );  

   }
   } else {
    X1_END_LOOP(k,j,i){
      r = x1[i];
      R = r*sin(x2[j]);
      z = r*cos(x2[j]);
      
      d_zero = density3D(R, z);
      
      d->Vc[VX1][k][j][i] = radvel(x1[i], x2[j]);
      d->Vc[VX2][k][j][i] = polvel(x1[i], x2[j]);
      
      d->Vc[RHO][k][j][i] = d_zero;
      d->Vc[VX3][k][j][i] = azivel(R,z);
      
#if EOS == IDEAL
      d->Vc[PRS][k][j][i] = d_zero*csq(R);
#endif
    }
    }
  }

  if (side == X2_BEG){
    double v_zero, d_zero;
    /*initial state*/
  
//    KTOT_LOOP(k){
//    for(j=0; j<=JBEG; j++){
//    ITOT_LOOP(i){
 
    X2_BEG_LOOP(k,j,i){
      r = x1[i];
      R = r*sin(x2[j]);
      z = r*cos(x2[j]);
      
      d_zero = density3D(R, z); 
      v_zero = azivel(R,z); 
      
      d->Vc[VX1][k][j][i] = radvel(x1[i], x2[j]);
      d->Vc[VX2][k][j][i] = polvel(x1[i], x2[j]);
      
      d->Vc[RHO][k][j][i] = d_zero; 
      d->Vc[VX3][k][j][i] = v_zero; 
      
#if EOS == IDEAL
      d->Vc[PRS][k][j][i] = d_zero*csq(R); 
#endif
//     }
//    }
    }
  }

  if (side == 0){
  double psi, psi_max, damp_rate;  
  double sint, cost, vR, dvrad;
  long int iseed, iseed_m;

  psi_max = CONST_PI/2.0 - theta_min; 

  //time for random pert?
//  print1 ("time is %f \n",g_time);
  if((planeton_switchon >= g_time) && (planeton_switchon <= g_time + g_dt) && (pert == 0)){
    pert = 1;
    iseed  =  x1[IEND/2]*x2[JEND/2]*x3[KEND/2]*1.0e4;
    iseed_m= -iseed;
    ran2(&iseed_m);

//random pert to radial velocity
  print1 ("do random pert... \n");
  TOT_LOOP(k,j,i){
    R = x1[i]*sin(x2[j]);

    dvrad = g_inputParam[pert_amp]*(2.0*ran2(&iseed) - 1.0)*sqrt(csq(R));
    dvrad*= exp(-0.5*pow((R - g_unitLength)/bigH(g_unitLength), 2.0));
    d->Vc[VX1][k][j][i] += dvrad;
  }
  }

   TOT_LOOP(k,j,i){

    sint = sin(x2[j]); cost = cos(x2[j]);

    r = x1[i];
    R = r*sint;
    z = r*cost;
    psi = CONST_PI/2.0 - x2[j];

    d_zero = density3D(R, z);
    p_zero = d_zero*csq(R);
    vphi_zero = azivel(R,z);
    vrad_zero = radvel(x1[i], x2[j]);
    vtheta_zero = polvel(x1[i], x2[j]);

    damp_rate = g_inputParam[tdamp]*vphi_zero/R;

    if((psi >= g_inputParam[damp_up])
       && (r > dampin)
       && (r < dampout)){
      
      taper = pow((psi -g_inputParam[damp_up] )/(psi_max - g_inputParam[damp_up]), 2.0);
      taper*= pow( (r - dampin)/(g_unitLength - dampin), 2.0);
      taper*= pow( (r - dampout)/(g_unitLength - dampout), 2.0);
      taper*= damp_rate;
      
      d->Vc[VX1][k][j][i] -= (d->Vc[VX1][k][j][i] - vrad_zero)*taper*dt;
      d->Vc[VX2][k][j][i] -= (d->Vc[VX2][k][j][i] - vtheta_zero)*taper*dt;
      d->Vc[VX3][k][j][i] -= (d->Vc[VX3][k][j][i] - vphi_zero)*taper*dt;
    
    }
    
    if(r <= dampin){
      
      taper = pow( (dampin - r)/(dampin - g_inputParam[rmin]), 2.0);
      taper*= damp_rate;
      
      d->Vc[VX1][k][j][i] -= (d->Vc[VX1][k][j][i] - vrad_zero)*taper*dt;
      d->Vc[VX2][k][j][i] -= (d->Vc[VX2][k][j][i] - vtheta_zero)*taper*dt;
//      d->Vc[VX3][k][j][i] -= (d->Vc[VX3][k][j][i] - vphi_zero)*taper*dt;
      
    } else if(r >= dampout){
      
      taper = pow( (r - dampout)/(g_inputParam[rout] - dampout), 2.0);
      taper*= damp_rate;
      
      d->Vc[VX1][k][j][i] -= (d->Vc[VX1][k][j][i] - vrad_zero)*taper*dt;
      d->Vc[VX2][k][j][i] -= (d->Vc[VX2][k][j][i] - vtheta_zero)*taper*dt;
//      d->Vc[VX3][k][j][i] -= (d->Vc[VX3][k][j][i] - vphi_zero)*taper*dt;
    } 
    
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
  double *dV1, *dV2, *dV3;
  double dV, mass, gmass, mc;

/*        
  dV1 = grid[IDIR].dV; 
  dV2 = grid[JDIR].dV; 
  dV3 = grid[KDIR].dV;

  mass = 0.0;
  DOM_LOOP(k,j,i){
    dV    = dV1[i]*dV2[j]*dV3[k];
    mass += dV*d->Vc[RHO][k][j][i];
  }
                        
#ifdef PARALLEL
  gmass = 0.;
  MPI_Allreduce (&mass, &gmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  mass = gmass;
#endif
        
  mc  = 0.5*g_inputParam[Mdisk]*CONST_Msun;
  mc /= g_unitDensity*g_unitLength*g_unitLength*g_unitLength*mass;
  DOM_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] *= mc;
    #if EOS == IDEAL
     d->Vc[PRS][k][j][i] *= mc;
    #endif

  }
*/
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
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi, mp, temp;

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

   double omega_p;
   t = g_time;
   /* if (t > 0.0) t += g_dt; */

//  if(g_intStage == 1) dt = g_dt;
  if(g_intStage == 2){
#if TIME_STEPPING == RK2 
    t += 0.5*g_dt;
#elif TIME_STEPPING == RK3
    t += 0.25*g_dt;
#endif
  }
  if(g_intStage == 3) t += 1.0/3.0*g_dt;

   omega_p = omega_k(g_unitLength);
   xp =  g_unitLength*cos(omega_p*t + CONST_PI);
   yp =  g_unitLength*sin(omega_p*t + CONST_PI);

   rsm = g_unitLength*pow(g_inputParam[mplanet]/mstar/3.0, 1.0/3.0);
   rsm *= g_inputParam[soft]; 
   d = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z + rsm*rsm);

   if(g_time < planeton){
      mp = 0.0;
   } else if( (g_time >= planeton) && (g_time < planeton_switchon)){
     temp = sin(0.5*CONST_PI*(g_time - planeton)/switchon);
     temp*= temp;
     mp   = g_inputParam[mplanet]*temp;
   } else {
     mp   = g_inputParam[mplanet];
   }

   phiplanet = -bigG*mp/d;

//planet indirect potential. rp=g_unitLength. no z-contrib (planet at midplane)
   phiplanet += bigG*mp*(x*xp + y*yp)/pow(g_unitLength, 3.0);

   phi = grav_pot3_cylin(R, z) + phiplanet;
  return phi;
}
#endif

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  *  \brief  Extracted from the Numerical Recipes in C (version 2) code.  
 *   *   Modified to use doubles instead of floats. - T. A. Gardiner - Aug. 12, 2003
 *    *   
 *     *
 *      * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 *       * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 *        * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 *         * values).  Call with idum = a negative integer to initialize;
 *          * thereafter, do not alter idum between successive deviates in a
 *           * sequence.  RNMX should appriximate the largest floating point value
 *            * that is less than 1. 
 *             */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


