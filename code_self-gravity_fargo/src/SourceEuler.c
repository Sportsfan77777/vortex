#include "mp.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

static PolarGrid *TemperInt;
static PolarGrid *VradNew, *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid *EnergyNew, *EnergyInt;
static real timeCRASH;  
extern boolean Corotating;
extern boolean Adiabatic;
extern boolean SelfGravity, SGZeroMode;
extern boolean ZMPlus;
real PhysicalTime=0.0, OmegaFrame, PhysicalTimeInitial;
int FirstGasStepFLAG=1;
static int AlreadyCrashed = 0, GasTimeStepsCFL;
extern boolean FastTransport, IsDisk;
Pair DiskOnPrimaryAcceleration;

boolean DetectCrash (array)
     PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0) 
	bool = YES;
    }
  }
  return bool;
}
 
void FillPolar1DArrays ()
{
  FILE *input, *output;
  int i,ii;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = PI*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
}

void InitEuler (Vr, Vt, Rho, Energy)
     PolarGrid *Vr, *Vt, *Rho, *Energy;
{
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
  Qplus        = CreatePolarGrid(NRAD, NSEC, "Qplus");
  InitComputeAccel ();
  /* Rho and Energy are already initialized: cf main.c */
  ComputeSoundSpeed (Rho, Energy);
  ComputePressureField (Rho, Energy);
  ComputeTemperatureField (Rho, Energy);
  InitGasVelocities (Vr, Vt);
}

real min2 (a,b)
     real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
     real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
     PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}


void AlgoGas (force, Rho, Vrad, Vtheta, Energy, Label, sys)
     Force *force;
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
     PlanetarySystem *sys;
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  int gastimestepcfl;
  boolean Crashed=NO;
  FirstGasStepFLAG=1;
  gastimestepcfl = 1;
  if (Adiabatic) {
    ComputeSoundSpeed (Rho, Energy);
    /* it is necesary to update computation of soundspeed if one uses
       alphaviscosity in FViscosity. It is not necesary in locally
       isothermal runs since cs is constant.  It is computed here for
       the needs of ConditionCFL. */
  }
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
    if (SloppyCFL == YES)
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;
  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    dtemp += dt;
    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      /* Indirect term star's potential computed here */
      DiskOnPrimaryAcceleration   = ComputeAccel (force, Rho, 0.0, 0.0, 0.0, 0.0);
      /* Gravitational potential from star and planet(s) is computed
	 and stored here*/
      FillForcesArrays (sys, Rho, Energy);
      /* Planets' velocities are updated here from gravitationnal
	 interaction with disk */
      AdvanceSystemFromDisk (force, Rho, Energy, sys, dt);
    }
    /* Planets' positions and velocities are updated from
       gravitational interaction with star and other planets */
    AdvanceSystemRK5 (sys, dt);
    /* Below we correct vtheta, planet's position and velocities if we
       work in a frame non-centered on the star */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES)
	CorrectVtheta (Vtheta, domega);
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    /* Now we update gas */
    if (IsDisk == YES) {
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);
      Crashed = DetectCrash (Rho);    /* test for negative density values */
      Crashed = DetectCrash (Energy);  /* test for negative energy values */
      if (Crashed == YES) {
	if (AlreadyCrashed == 0) {
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
	  WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
	  WriteDiskPolar (Vtheta, 999); /* of what happened */
	  WriteDiskPolar (Temperature, 999);
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (".");
      }
      fflush (stdout);
      if (ZMPlus)
	compute_anisotropic_pressurecoeff (sys);
      ComputePressureField (Rho, Energy);
      SubStep1 (Vrad, Vtheta, Rho, dt);
      SubStep2 (Rho, Energy, dt);
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);
      if (Adiabatic) {
	ComputeViscousTerms (Vrad, Vtheta, Rho);
	SubStep3 (Rho, dt);
	ActualiseGas (Energy, EnergyNew);
      }
      Transport (Rho, Vrad, Vtheta, Energy, Label, dt);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);
      ComputeTemperatureField (Rho, Energy);
      mdcp = CircumPlanetaryMass (Rho, sys);
      exces_mdcp = mdcp - mdcp0;
    }
    PhysicalTime += dt;
  }
  masterprint ("\n");
}

void SubStep1 (Vrad, Vtheta, Rho, dt)
     PolarGrid *Vrad, *Vtheta, *Rho;
     real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  boolean selfgravityupdate;
  extern boolean Evanescent;
  real *vrad, *vtheta, *rho;
  real *Pot, *Press;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
  Press = Pressure->Field;
  /* In this substep we take into account  */
  /* the source part of Euler equations  */
  /* We evolve velocities with pressure gradients */
  /* gravitational forces and curvature terms */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	gradp = (Press[l]-Press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
	vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
	vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
	vt2 = vt2*vt2;
	vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	gradp = (Press[l]-Press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	if ( ZMPlus )
	  gradp *= SG_aniso_coeff;
	gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
	vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
	vthetaint[l] += dt*supp_torque;
      }
    }
  }
  if ( SelfGravity ) {
    selfgravityupdate = YES;
    compute_selfgravity(Rho, VradInt, VthetaInt, dt, selfgravityupdate);
  }
  ComputeViscousTerms (VradInt, VthetaInt, Rho);
  UpdateVelocitiesWithViscosity (VradInt, VthetaInt, Rho, dt);
  if ( !Evanescent )
    ApplySubKeplerianBoundary (VthetaInt);
}

void SubStep2 (Rho, Energy, dt)
     PolarGrid *Rho, *Energy;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy;
  real *vradnew, *vthetanew, *qt, *qr, *energyint;
  real dxtheta, invdxtheta;
  real dv;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	qt[l] = 0.0;
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
    /* If gas disk is adiabatic, we add artificial viscosity as a source */
    /* term for advection of thermal energy polargrid */
    if (Adiabatic) {
#pragma omp for nowait
      for (i = 0; i < nr; i++) {
	dxtheta = 2.0*PI/(real)ns*Rmed[i];
	invdxtheta = 1.0/dxtheta;
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  lip = l+ns;
	  ljp = l+1;
	  if (j == ns-1) ljp = i*ns;
	  energyint[l] = energy[l] -				\
	    dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -	\
	    dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
	}
      }
    }
  }
}
	       
void SubStep3 (Rho, dt)
     PolarGrid *Rho;
     real dt;
{
  extern boolean Cooling;
  int i, j, l, nr, ns;
  int lip, li2p;
  real r, rip, ri2p, qpip, qpi2p;
  real *dens, *pres, *energy, *energynew;
  real *divergence, *Trr, *Trp, *Tpp, *qplus;
  real viscosity, energypred, num, den;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energy = EnergyInt->Field;
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  /* In this substep we take into account  */
  /* the source part of energy equation  */
  /* We evolve internal energy with  */
  /* compression/dilatation and heating terms */
#pragma omp parallel private(j,l,viscosity,energypred)
  {
#pragma omp for
    /* We calculate the heating source term Qplus from i=1 */
    for (i = 1; i < nr; i++) { /* Trp defined from i=1 */
      viscosity = FViscosity (Rmed[i]); /* Global_Bufarray holds cs */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	if (viscosity != 0.0) {
	  qplus[l] = 0.5/viscosity/dens[l]*(Trr[l]*Trr[l] +	\
					    Trp[l]*Trp[l] +	\
					    Tpp[l]*Tpp[l] );
	  qplus[l] += (2.0/9.0)*viscosity*dens[l]*divergence[l]*divergence[l];
	}
	else
	  qplus[l] = 0.0;
      }
    }
    /* We calculate the heating source term Qplus for i=0 */
    i = 0;
    r    = Rmed[i];
    rip  = Rmed[i+1];
    ri2p = Rmed[i+2];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      li2p = lip+ns;
      qpip = qplus[lip];   /* qplus(i=1,j) */
      qpi2p = qplus[li2p]; /* qplus(i=2,j) */
      if (viscosity != 0.0) {
	/* power-law extrapolation */
	qplus[l] = qpip*exp( log(qpip/qpi2p) * log(r/rip) / log(rip/ri2p) );
      }
      else
	qplus[l] = 0.0;
    }
    /* Now we can update energy with source terms from i=0 */
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	if (!Cooling) {
	  num = dt*qplus[l] + energy[l];
	  den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l];
	  energynew[l] = num/den;
	}
	if (Cooling) {
	  num = EnergyMed[i]*dt*dens[l]/SigmaMed[i] + CoolingTimeMed[i]*energy[l] + dt*CoolingTimeMed[i]*(qplus[l]-QplusMed[i]*dens[l]/SigmaMed[i]);
	  den = dt + CoolingTimeMed[i] + (ADIABATICINDEX-1.0)*dt*CoolingTimeMed[i]*divergence[l];
	  energynew[l] = num/den;
	}
      }
    }
  }
}
   		   
int ConditionCFL (Vrad, Vtheta, deltaT)
     PolarGrid *Vrad, *Vtheta;
     real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, cs, newdt, dt;
  int ideb, jdeb;
  real itdbg1, itdbg2, itdbg3, itdbg4, mdtdbg; /* debugging variables */
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr, visct;
  real *soundspeed;
  soundspeed = SoundSpeed->Field;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
	Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
	Vresidual[j] = vt[l];	       /* Standard algorithm */
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = soundspeed[l];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4);
      if (dt < newdt) {
	newdt = dt;
	if (debug == YES) {
	  ideb = i;
	  jdeb = j;
	  itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4;
	  mdtdbg = newdt;
	  viscr = dxrad/dvr/4.0/CVNR/CVNR;     
	  visct = dxtheta/dvt/4.0/CVNR/CVNR;
	}
      }  
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Viscosity limit                : %g\n", itdbg4);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (int)(ceil(deltaT/newdt));
}

void ComputeSoundSpeed (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns, ii;
  real *dens, *energ, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic)
	cs[l] = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      else
	cs[l] = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]/dens[l] );
    }
  }
}

void ComputePressureField (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns, ii;
  real *dens, *pres, *energ, *cs;
  real peq, dp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic) {
	pres[l] = dens[l]*cs[l]*cs[l]; /* since SoundSpeed is not updated */
                                       /* from initialization, cs remains */ 
                                       /* axisymmetric */
      }
      else
	pres[l] = (ADIABATICINDEX-1.0)*energ[l];
    }
  }
}

void ComputeTemperatureField (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *temp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  temp = Temperature->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic)
	temp[l] = MU/R* pres[l]/dens[l];
      else
	temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
    }
  }
}

real CircumPlanetaryMass (Rho, sys)
     PolarGrid *Rho;
     PlanetarySystem *sys;
{
  int i, j, l, ns;
  real xpl, ypl;
  real dist, mdcplocal, mdcptotal;
  real *dens, *abs, *ord;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  xpl = sys->x[0];
  ypl = sys->y[0];
  mdcplocal = 0.0;
  mdcptotal = 0.0;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  for ( i = Zero_or_active; i < Max_or_active; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      dist = sqrt ( (abs[l]-xpl)*(abs[l]-xpl) +		\
		    (ord[l]-ypl)*(ord[l]-ypl) );
      if ( dist < HillRadius ) {
	mdcplocal += Surf[i] * dens[l];
      }
    }
  }
  if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&mdcplocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    mdcptotal = mdcplocal;
  }
  return mdcptotal;
}
