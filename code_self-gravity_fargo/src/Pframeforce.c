#include "mp.h"

extern boolean AllowAccretion, Corotating, Indirect_Term, Cooling;
extern Pair DiskOnPrimaryAcceleration;
static Pair IndirectTerm;
static real q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static real vt_int[MAX1D], vt_cent[MAX1D];

void ComputeIndirectTerm () {
  IndirectTerm.x = -DiskOnPrimaryAcceleration.x;
  IndirectTerm.y = -DiskOnPrimaryAcceleration.y; 
  if (Indirect_Term == NO) {
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
  }
}
/* Below : work in non-rotating frame */
/* centered on the primary */
void FillForcesArrays (sys, Rho, Energy)
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns, k, NbPlanets;
  real x, y, angle, distance, distancesmooth, iplanet;
  real xplanet, yplanet, RRoche,smooth, mplanet;
  real PlanetDistance, *Pot, pot, smoothing;
  real InvPlanetDistance3, InvDistance;
  real xbin, ybin, mbin, distbin, Invdistbin3;
  Pot= Potential->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  NbPlanets = sys->nb;

  /* Indirect term star on gas here */
  ComputeIndirectTerm ();

#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) Pot[i] = 0.0;
  /* -- Gravitational potential from planet on gas -- */
  for (k = 0; k < NbPlanets; k++) {
    xplanet = sys->x[k];
    yplanet = sys->y[k];
    mplanet = sys->mass[k]*MassTaper;
    PlanetDistance = sqrt(xplanet*xplanet+yplanet*yplanet);
    InvPlanetDistance3 =  1.0/PlanetDistance/PlanetDistance/PlanetDistance;
    RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
    if (RocheSmoothing)
      smoothing = RRoche*ROCHESMOOTHING;
    else
      smoothing = compute_smoothing (PlanetDistance);
    smooth = smoothing*smoothing;
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,distance,distancesmooth,pot)
    for (i = 0; i < nr; i++) {
      InvDistance = 1.0/Rmed[i];
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	angle = (real)j/(real)ns*2.0*PI;
	x = Rmed[i]*cos(angle);
	y = Rmed[i]*sin(angle);
	distance = (x-xplanet)*(x-xplanet)+(y-yplanet)*(y-yplanet);
	distancesmooth = sqrt(distance+smooth);
	pot = -G*mplanet/distancesmooth; /* Direct term from planet */
	if (Indirect_Term == YES)
	  pot += G*mplanet*InvPlanetDistance3*(x*xplanet+y*yplanet); /* Indirect term from planet  */
	Pot[l] += pot;
      }
    }
  }
  /* -- Gravitational potential from star on gas -- */
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,pot)
  for (i = 0; i < nr; i++) {
    InvDistance = 1.0/Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      angle = (real)j/(real)ns*2.0*PI;
      x = Rmed[i]*cos(angle);
      y = Rmed[i]*sin(angle);
      pot = -G*1.0*InvDistance;  /* Direct term from star */
      pot -= IndirectTerm.x*x + IndirectTerm.y*y; /* Indirect term from star */
      Pot[l] += pot;	
    }
  }
} 

void AdvanceSystemFromDisk (force, Rho, Energy, sys, dt)
     Force *force;
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy;
     real dt;		       
{
  int NbPlanets, k;
  Pair gamma;
  real x, y, r, m, smoothing;
  NbPlanets = sys->nb;
  for (k = 0; k < NbPlanets; k++) {
    if (sys->FeelDisk[k] == YES) {
      m=sys->mass[k]*MassTaper + sys->accreted_mass[k];
      x = sys->x[k];
      y = sys->y[k];
      r = sqrt(x*x + y*y);
      if (RocheSmoothing)
	smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
      else
	smoothing = compute_smoothing (r);
      gamma = ComputeAccel (force, Rho, x, y, smoothing, m);
      sys->vx[k] += dt * gamma.x;
      sys->vy[k] += dt * gamma.y;
      sys->vx[k] += dt * IndirectTerm.x;
      sys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

void AdvanceSystemRK5 (sys, dt)
     PlanetarySystem *sys;
     real dt;
{
  extern boolean ForcedCircular;
  int i, n;
  boolean *feelothers;
  real dtheta, omega, rdot, x, y, r, new_r, vx, vy, theta, denom;
  n = sys->nb;
  if ( !ForcedCircular ) {
    for (i = 0; i < n; i++) {
      q0[i] = sys->x[i];
      q0[i+n] = sys->y[i];
      q0[i+2*n] = sys->vx[i];
      q0[i+3*n] = sys->vy[i];
      PlanetMasses[i] = sys->mass[i]*MassTaper + sys->accreted_mass[i]; /*** ##### MASS TAPER EDIT HERE  ##### ***/
    }
    feelothers = sys->FeelOthers;
    RungeKunta (q0, dt, PlanetMasses, q1, n, feelothers);
  }
  for (i = 1-(PhysicalTime >= RELEASEDATE); i < sys->nb; i++) {
    if ( !ForcedCircular ) {
      sys->x[i] = q1[i];
      sys->y[i] = q1[i+n];
      sys->vx[i] = q1[i+2*n];
      sys->vy[i] = q1[i+3*n];
    }
    else {
      x = sys->x[i];
      y = sys->y[i];
      theta = atan2(y,x);
      vx = sys->vx[i];
      vy = sys->vy[i];
      r = sqrt(x*x + y*y);
      omega = (-y*vx + x*vy)/r/r;
      dtheta = omega*dt;
      sys->x[i] = r*cos(theta+dtheta);
      sys->y[i] = r*sin(theta+dtheta);
      sys->vx[i]= vx*cos(dtheta+theta) - vy*sin(dtheta+theta);
      sys->vy[i]= vx*sin(dtheta+theta) + vy*cos(dtheta+theta);
    }
  }
  if (PhysicalTime < RELEASEDATE) {
    x = sys->x[0];
    y = sys->y[0];
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    rdot = (RELEASERADIUS-r)/(RELEASEDATE-PhysicalTime);
    omega = sqrt((1.+sys->mass[0])/r/r/r);
    new_r = r + rdot*dt;
    denom = r-new_r;
    if (denom != 0.0) {
      dtheta = 2.*dt*r*omega/denom*(sqrt(r/new_r)-1.);
    } else {
      dtheta = omega*dt;
    }
    vx = rdot;
    vy = new_r*sqrt((1.+sys->mass[0])/new_r/new_r/new_r);
    sys->x[0] = new_r*cos(dtheta+theta);
    sys->y[0] = new_r*sin(dtheta+theta);
    sys->vx[0]= vx*cos(dtheta+theta) - vy*sin(dtheta+theta); 
    sys->vy[0]= vx*sin(dtheta+theta) + vy*cos(dtheta+theta); 
  }
}

void SolveOrbits (sys)
     PlanetarySystem *sys;
{
  int i, n;
  real x, y, vx, vy;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    x = sys->x[i];
    y = sys->y[i];
    vx = sys->vx[i];
    vy = sys->vy[i];
    FindOrbitalElements (x, y, vx, vy, 1.0+sys->mass[i], i);
  }
} 

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }
  return lapl;
}

void InitGasDensity (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  real *dens;
  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  FillSigma ();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      dens[l] = SigmaMed[i];
    }
  }
}

void InitGasEnergy (Energy)
     PolarGrid *Energy;
{
  int i, j, l, nr, ns;
  real *energy;
  energy = Energy->Field;
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  FillEnergy ();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energy[l] = EnergyMed[i];
    }
  }
}

void InitGasVelocities (Vr, Vt)
     PolarGrid *Vr, *Vt;
{
  extern boolean SGZeroMode;
  extern boolean SelfGravity;
  int i, j, l, nr, ns;
  real *vr, *vt, *pres, *cs;
  real  r, omega, ri;
  real viscosity, t1, t2, r1, r2;
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Vt->Nrad;
  ns  = Vt->Nsec;
  cs = SoundSpeed->Field;
  pres = Pressure->Field;  /* Pressure is already initialized: cf initeuler in SourceEuler.c ... */
  /* --------- */
  // Initialization of azimutal velocity with exact centrifugal balance
  /* --------- */
  if ( CentrifugalBalance ) {
    /* vt_int \equiv rOmega² = grad(P)/sigma +  \partial(phi)/\partial(r)  -  acc_sg_radial */
    mpi_make1Dprofile (pres, GLOBAL_bufarray);
    /* global axisymmetric pressure field, known by all cpus*/
    for (i = 1; i < GLOBALNRAD; i++) {
      vt_int[i] = ( GLOBAL_bufarray[i] - GLOBAL_bufarray[i-1] ) /	\
	(.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1]) + \
	G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    }
    /* Case of a disk with self-gravity */
    if ( SelfGravity ) { // Better test with CL rigid!
      if ( !SGZeroMode )
	  mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      	else
	  GLOBAL_AxiSGAccr = SG_Accr;
      for (i = 1; i < GLOBALNRAD; i++)
	vt_int[i] -= ( (Radii[i] - GlobalRmed[i-1])*GLOBAL_AxiSGAccr[i] + \
		       (GlobalRmed[i] - Radii[i])*GLOBAL_AxiSGAccr[i-1] ) / (GlobalRmed[i]-GlobalRmed[i-1]);
    }
    for (i = 1; i < GLOBALNRAD; i++)
      vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
    
    t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
    r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
    t2 = vt_cent[0];
    r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    t1 = t1-r1/(r2-r1)*(t2-t1);
    vt_cent[0] = t1;
    ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[GLOBALNRAD] = vt_cent[GLOBALNRAD-1];
  }
  /* --------- */
  // Initialization with self-gravity, without exact centrifugal balance
  if (SelfGravity && !CentrifugalBalance)
    init_azimutalvelocity_withSG (Vt);
  /* --------- */
  if (ViscosityAlpha)
    mpi_make1Dprofile (cs, GLOBAL_bufarray);
  /* We calculate here the cooling time radial profile (Theo.c) */
  if (Cooling) {
    FillCoolingTime();
    /* To fill qplus, one requires to calculate viscosity, hence cs if
       one uses alpha-viscosity */
    FillQplus();
  }
  
  for (i = 0; i <= nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    }
    else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    viscosity = FViscosity (r);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* --------- */
      if (!SelfGravity) {
	omega = sqrt(G*1.0/r/r/r);
	vt[l] = omega*r*sqrt(1.0-pow(ASPECTRATIO,2.0)*			\
			     pow(r,2.0*FLARINGINDEX)*			\
			     (1.+SIGMASLOPE-2.0*FLARINGINDEX) );
      }
      /* --------- */
      vt[l] -= OmegaFrame*r;
      if (CentrifugalBalance)
	vt[l] = vt_cent[i+IMIN];
      if (i == nr)
	vr[l] = 0.0;
      else {
	vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
	if (ViscosityAlpha) {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
	} else {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
	}
      }
    }
  }
  for (j = 0; j < ns; j++)
    vr[j] = vr[j+ns*nr] = 0.0;
}
