/** \file SourceEuler.c 

Contains routines used by the hydrodynamical loop. More specifically,
it contains the main loop itself and all the source term substeps
(with the exception of the evaluation of the viscous force). The
transport substep is treated elsewhere. */

#include "fargo.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

static PolarGrid *TemperInt, *DTemperInt;
static PolarGrid *VradNew,   *VradInt, *DVradNew,   *DVradInt;
static PolarGrid *VthetaNew, *VthetaInt, *DVthetaNew, *DVthetaInt;
static real timeCRASH;  
extern boolean Corotating;

static int AlreadyCrashed = 0;
static long GasTimeStepsCFL;

extern int TimeStep;
extern boolean FastTransport, IsDisk, VortexDiffusion;
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
      if (ptr[l]<DENSITYFLOOR)ptr[l]=DENSITYFLOOR;
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
      fprintf (output, "%.18g %.18g\n", Radii[i], GlobalRmed[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
}

void InitEuler (Rho, Vr, Vt, DRho, DVr, DVt)
PolarGrid *Rho, *Vr, *Vt, *DRho, *DVr, *DVt;
{
  FillPolar1DArrays ();
  FillSigma ();
  InitTransport ();
  InitViscosity ();

  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");

  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");

  DRhoStar      = CreatePolarGrid(NRAD, NSEC, "DRhoStar"); /*initialize dust density star/int and vradint vradnew and vthetanew vthetaint */
  DRhoInt       = CreatePolarGrid(NRAD, NSEC, "DRhoInt");
  DVradNew      = CreatePolarGrid(NRAD, NSEC, "DVradNew");
  DVradInt      = CreatePolarGrid(NRAD, NSEC, "DVradInt");
  DVthetaNew    = CreatePolarGrid(NRAD, NSEC, "DVthetaNew");
  DVthetaInt    = CreatePolarGrid(NRAD, NSEC, "DVthetaInt");
  DTemperInt    = CreatePolarGrid(NRAD, NSEC, "DTemperInt");

  TStop		= CreatePolarGrid(NRAD, NSEC, "TStop"); /* dust stoping time */
  Diag1         = CreatePolarGrid(NRAD, NSEC, "Diag1"); /* three diagnostic */
  Diag2         = CreatePolarGrid(NRAD, NSEC, "Diag2");
  Diag3         = CreatePolarGrid(NRAD, NSEC, "Diag3");

  if (DustDiff == YES) {
    Fdiffrp         = CreatePolarGrid(NRAD, NSEC, "Fdiffrp");
    Fdifftp         = CreatePolarGrid(NRAD, NSEC, "Fdifftp");
  }
  InitGas (Rho, Vr, Vt, DRho, DVr, DVt);
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

void AlgoGas (Rho,Vrad,Vtheta,DRho, DVrad, DVtheta, Label,sys)
PolarGrid *Rho, *Vrad, *Vtheta, *DRho, *DVrad, *DVtheta, *Label;
PlanetarySystem *sys; 
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  long gastimestepcfl;
  boolean Crashed=NO;
  gastimestepcfl = 1;
  if (IsDisk == YES) { 
    CommunicateBoundaries (Rho,Vrad,Vtheta,DRho,DVrad,DVtheta,Label);
    if (SloppyCFL == YES)
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, DVrad, DVtheta, DT-dtemp); /*dust CFL condition included */
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;
  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho,Vrad,Vtheta,DRho,DVrad,DVtheta,Label);
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DVrad, DVtheta, DT-dtemp); /* dust CFL condition included*/
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta,DRho, DVrad, DVtheta, dt, sys); /* Dust can accrete but its mass not added to the planet and won't affect planet momentum */
    }
    dtemp += dt;
    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {	/*Dust's effect onto the planet movement is ignored */
      DiskOnPrimaryAcceleration   = ComputeAccel (Rho, 0.0, 0.0, 0.0, 0.0);
      FillForcesArrays (sys);
      AdvanceSystemFromDisk (Rho, sys, dt);  /* the planets move due to the force from the disk */
    }
    AdvanceSystemRK5 (sys, dt); /* planets move due to their interaction */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES){
	CorrectVtheta (Vtheta, domega);
	CorrectVtheta (DVtheta, domega);  /*In Corotating frame, both gas and dust vtheta needs to be changed*/ 
      }
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    if (IsDisk == YES) {
      ApplyBoundaryCondition(Vrad, Vtheta, Rho, DVrad, DVtheta, DRho,dt); /*dust boundary condition is included*/
      Crashed = DetectCrash (Rho);  /* test for negative density values */
      Crashed = DetectCrash (DRho); 
      if (Crashed == YES) {
	if (AlreadyCrashed == 0) {
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
	  WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
	  WriteDiskPolar (Vtheta, 999); /* of what happened */
	  WriteDiskPolar (DRho, 999);    /* We write the HD arrays */
          WriteDiskPolar (DVrad, 999);   /* in order to keep a track */
          WriteDiskPolar (DVtheta, 999); /* of what happened */
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (".");
      }
      fflush (stdout);
      SubStep1 (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, dt);
      SubStep2 (Rho, DRho, dt);
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      if(!SFTA){
        ActualiseGas (DVrad, DVradNew); /*actualize dust vrad and vtheta */
        ActualiseGas (DVtheta, DVthetaNew);
      }
      if(ADDM) AddMass (Rho, Vrad, Vtheta, dt);
      CommunicateBoundaries (Rho,Vrad,Vtheta,DRho,DVrad,DVtheta,Label);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, dt);
      Transport (Rho, Vrad, Vtheta, Label, dt);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, dt);
      if(SFTA) SFTAvelocity(Rho,Vrad,Vtheta, DRho, DVrad, DVtheta);
      Transportd (DRho, Rho, DVrad, DVtheta, Label, dt);
      if(DustDiff == YES) Diffd(DRho, Rho, DVrad, DVtheta,dt);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, dt);
      SFTAtest(Rho,Vrad,Vtheta, DRho, DVrad, DVtheta);
    }
    PhysicalTime += dt;
  }
  masterprint ("\n");
}

void SFTAvelocity (Rho, Vrad, Vtheta, DRho, DVrad, DVtheta)
PolarGrid *Rho, *DRho, *Vrad, *Vtheta,*DVrad,*DVtheta;
{
  int i, j, l, nr, ns, lim, ljm, lip, ljp;
  real *rho, *drho, *vrad, *vtheta, *dvrad, *dvtheta, *ts, *diag1, *diag2, *diag3;
  real cs2,cs2m,gradp, dxtheta, invdxtheta;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dvrad = DVrad->Field;
  dvtheta = DVtheta->Field;
  ts = TStop->Field;
  diag1 = Diag1->Field;
  diag2 = Diag2->Field;
  diag3 = Diag3->Field;
  for (i = 1; i < nr; i++){
      for (j = 0; j < ns; j++){
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        cs2 = SOUNDSPEED[i]*SOUNDSPEED[i];
        cs2m= SOUNDSPEED[i-1]*SOUNDSPEED[i-1];
        gradp = (cs2*rho[l]-cs2m*rho[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	dvrad[l]=vrad[l]+1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]*gradp;
//	diag1[l]=rho[l];
//	diag2[l]=(dvrad[l]-vrad[l])/(1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]*(gradp+1.e-20));
//   	diag2[l]=dvrad[l];
//	diag3[l]=vrad[l];
      }
  }

  for (i = 0; i < nr; i++) { 
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
	cs2   = SOUNDSPEED[i]*SOUNDSPEED[i];
        gradp = cs2*(rho[l]-rho[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	dvtheta[l]=vtheta[l]+1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]*gradp;
//        diag3[l]=dvtheta[l];
      }
  }
}

void SFTAtest (Rho, Vrad, Vtheta, DRho, DVrad, DVtheta)
PolarGrid *Rho, *DRho, *Vrad, *Vtheta,*DVrad,*DVtheta;
{
  int i, j, l, nr, ns, lim, ljm, lip, ljp;
  real *rho, *drho, *vrad, *vtheta, *dvrad, *dvtheta, *ts, *diag1, *diag2, *diag3;
  real cs2,cs2m,gradp, dxtheta, invdxtheta;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dvrad = DVrad->Field;
  dvtheta = DVtheta->Field;
  ts = TStop->Field;
  diag1 = Diag1->Field;
  diag2 = Diag2->Field;
  diag3 = Diag3->Field;
  for (i = 1; i < nr; i++){
      for (j = 0; j < ns; j++){
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        cs2 = SOUNDSPEED[i]*SOUNDSPEED[i];
        cs2m= SOUNDSPEED[i-1]*SOUNDSPEED[i-1];
        gradp = (cs2*rho[l]-cs2m*rho[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
//        diag2[l]=(vrad[l]-dvrad[l])/(1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]*(gradp+1.e-20));
//	diag2[l]=(diag2[l]-dvrad[l])/(dvrad[l]+1.e-20);
//	diag3[l]=(diag3[l]-dvtheta[l])/(dvtheta[l]+1.e-20);
//	diag3[l]=dvrad[l]-vrad[l];
//	diag2[l]=rho[l];
      }
  }

  for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        cs2   = SOUNDSPEED[i]*SOUNDSPEED[i];
        gradp = cs2*(rho[l]-rho[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
//        diag3[l]=(vtheta[l]-dvtheta[l])/(1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]*(gradp+1.e-20));
      }
  }
}

void Diffd (DRho, Rho, DVrad, DVtheta,dt)
PolarGrid *DRho, *Rho, *DVrad, *DVtheta;
real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *drho, *rho, *dvrad, *dvtheta, *fdiffrp, *fdifftp, *diag3;
  real dust_overdensity;
  real dtheta, invdtheta;
  real viscosityp, viscosity;
  real frp, frm, ftp, ftm;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  rho = Rho->Field;
  drho = DRho->Field;
  dvrad = DVrad->Field;
  dvtheta = DVtheta->Field;
  fdiffrp = Fdiffrp->Field;
  fdifftp = Fdifftp->Field;
  diag3 = Diag3->Field;

  for (i = 0; i < nr-1; i++){
    dtheta = 2.0*PI/(real)ns;
    invdtheta = 1.0/dtheta;
    //viscosityp = FViscosity (Rsup[i]);
    //viscosity = FViscosity (Rmed[i]);

    for (j = 0; j < ns; j++){
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      /* Adding Separate Vortex Diffusion */
      if (VortexDiffusion == YES) {
        if (Rmed[i] > VORTEXDIFFIN && Rmed[i] < VORTEXDIFFOUT) {
          dust_overdensity = drho[l] / DSigmaMed[i] 
          // check if 'current divided by initial' overdensity exceeds threshold
          if (dust_overdensity > VORTEXDIFFUSIONTHRESHOLD) {
            viscosityp = VORTEXDIFFUSIONCOEFFICIENT;
            viscosity = VORTEXDIFFUSIONCOEFFICIENT;  
          }
          else if (dust_overdensity > VORTEXDIFFUSIONLOWERTHRESHOLD) {
            // Ramp: e^(log(10^4) (sin((pi / 2)(x-1.5) / 0.5))^2) 
            viscosityp = FViscosity (Rsup[i]) * exp(log(VORTEXDIFFUSIONCOEFFICIENT / FViscosity(Rsup[i])) * pow(sin((PI / 2) * (dust_overdensity - VORTEXDIFFUSIONLOWERTHRESHOLD) / (VORTEXDIFFUSIONTHRESHOLD - VORTEXDIFFUSIONLOWERTHRESHOLD)), 2))
            viscosity = FViscosity (Rmed[i]) * exp(log(VORTEXDIFFUSIONCOEFFICIENT / FViscosity(Rmed[i])) * pow(sin((PI / 2) * (dust_overdensity - VORTEXDIFFUSIONLOWERTHRESHOLD) / (VORTEXDIFFUSIONTHRESHOLD - VORTEXDIFFUSIONLOWERTHRESHOLD)), 2))
          }
          else {
            viscosityp = FViscosity (Rsup[i]);
            viscosity = FViscosity (Rmed[i]);
          }
        }
      }
      /* End of Separate Vortex Diffusion */  

      fdiffrp[l]=-viscosityp*(rho[lip]+rho[l])/2.*(drho[lip]/rho[lip]-drho[l]/rho[l])*InvDiffRmed[i+1];
      fdifftp[l]=-viscosity*(rho[ljp]+rho[l])/2.*(drho[ljp]/rho[ljp]-drho[l]/rho[l])*invdtheta/Rmed[i];
      diag3[l]=-viscosityp*(log(drho[lip]/rho[lip])-log(drho[l]/rho[l]))/(log(Rmed[i+1])-log(Rmed[i]))/Rsup[i];
    }
  } 
  for (i = 1; i < nr-1; i++){
    if(Rmed[i]>DUSTDIFFINN){
      for (j = 0; j < ns; j++){
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        frp=fdiffrp[l];
        frm=fdiffrp[lim];
        ftp=fdifftp[l];
        ftm=fdifftp[ljm];
        drho[l]=drho[l]-(Rsup[i]*dtheta*frp-Rinf[i]*dtheta*frm+(Rsup[i]-Rinf[i])*(ftp-ftm))*dt*InvSurf[i];
      }
    }
  }
}

void SubStep1 (Vrad, Vtheta, Rho, DVrad, DVtheta, DRho, dt)
PolarGrid *Vrad, *Vtheta, *Rho, *DVrad, *DVtheta, *DRho;
real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *dvrad, *dvtheta, *drho, *ts, *diag1, *diag2;
  real *vradint, *vthetaint, *dvradint, *dvthetaint;
  real gradp, gradphi, vt2, dvt2, dxtheta, cs2, cs2m, dcs2, dcs2m, dgradp, dforcer, dforcet;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  real *Pot;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;

  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;

  drho = DRho->Field;
  dvrad = DVrad->Field; 
  dvtheta = DVtheta->Field;
  dvradint = DVradInt->Field;
  dvthetaint = DVthetaInt->Field;

  Pot = Potential->Field;
  ts = TStop->Field;
  diag1 = Diag1->Field;
  diag2 = Diag2->Field;
				/* In this substep we take into account     */
				/* the source part of Euler equations       */
				/* (i.e. the R.H.S. in classical notation). */
#pragma omp parallel private(j,l,lim,ljm,ljp,cs2,cs2m,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
    if(GasDcouple == YES){
      for (i = 0; i < nr; i++){
        for (j = 0; j < ns; j++){
	  l = j+i*ns;
	  ts[l]=2.813e-6*PSize/rho[l]*PI/2.;
        }
      }
    } 
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	cs2 = SOUNDSPEED[i]*SOUNDSPEED[i];
	cs2m= SOUNDSPEED[i-1]*SOUNDSPEED[i-1];
	gradp = (cs2*rho[l]-cs2m*rho[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
        vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
        vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
        vt2 = vt2*vt2;
	vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
	if(!SFTA){			/* dust */
	  dcs2 = DSOUNDSPEED[i]*DSOUNDSPEED[i];
          dcs2m= DSOUNDSPEED[i-1]*DSOUNDSPEED[i-1];
	  dgradp = (dcs2*drho[l]-dcs2m*drho[lim])*2.0/(drho[l]+drho[lim])*InvDiffRmed[i];
	  dvt2 = dvtheta[l]+dvtheta[ljp]+dvtheta[lim]+dvtheta[ljp-ns];
          dvt2 = dvt2/4.0+Rinf[i]*OmegaFrame;
          dvt2 = dvt2*dvt2;
	  dvradint[l] =dvrad[l]+dt*(-gradphi-dgradp+dvt2*InvRinf[i]);
	  if(GasDcouple == YES){
	    dforcer = -1./(1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]+dt/2.)*(dvrad[l]-vrad[l]);
	    dvradint[l] += dt*dforcer;
	    diag1[l] = dforcer;
	  }
	}
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	cs2m = cs2 = SOUNDSPEED[i]*SOUNDSPEED[i];
	gradp = cs2*(rho[l]-rho[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
	vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
        vthetaint[l] += dt*supp_torque; /* imposeddiskdrift cannot be used for dust */
	if(!SFTA){		/* dust */
          dcs2m = dcs2 = DSOUNDSPEED[i]*DSOUNDSPEED[i];
          dgradp=dcs2*(drho[l]-drho[ljm])*2.0/(drho[l]+drho[ljm])*invdxtheta;
	  dvthetaint[l] = dvtheta[l]-dt*(gradphi+dgradp);
	  if(GasDcouple == YES){
	    dforcet = -1./(1./sqrt(G*1.0/Rmed[i])*Rmed[i]*ts[l]+dt/2.)*(dvtheta[l]-vtheta[l]);
	    dvthetaint[l] += dt*dforcet; 
	    diag2[l] = dforcet;
	  }
	}
      }
    }
  }
  ViscousTerms (VradInt, VthetaInt, Rho, DVradInt, DVthetaInt, DRho, dt);
}

void SubStep2 (Rho, DRho, dt)
PolarGrid *Rho, *DRho;
real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *dvrad, *dvtheta, *drho;
  real *vradnew, *vthetanew, *dvradnew, *dvthetanew, *qt, *qr, *dqt, *dqr;
  real dxtheta, invdxtheta;
  real dv, ddv;
  nr = VradInt->Nrad;
  ns = VradInt->Nsec;

  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;

  drho = DRho->Field;
  dvrad = DVradInt->Field;
  dvtheta = DVthetaInt->Field;
  dqr = DRhoInt->Field;
  dvradnew = DVradNew->Field;
  dvthetanew = DVthetaNew->Field;
  dqt = DTemperInt->Field;
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lim = l-ns;
      lip = l+ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      ddv = dvrad[lip]-dvrad[l];  /*dust*/
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 

      if (ddv < 0.0)
        dqr[l] = CVNR*CVNR*drho[l]*ddv*ddv;
      else
        dqr[l] = 0.0;

      dv = vtheta[ljp]-vtheta[l];
      ddv = dvtheta[ljp]-vtheta[l]; /*dust */
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	qt[l] = 0.0;

      if (ddv < 0.0) 
        dqt[l] = CVNR*CVNR*drho[l]*ddv*ddv;
      else
        dqt[l] = 0.0;
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	lip = l+ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
	dvradnew[l] = dvrad[l]-dt*2.0/(drho[l]+drho[lim])*(dqr[l]-dqr[lim])*InvDiffRmed[i]; /* dust */
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	lip = l+ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
	dvthetanew[l] = dvtheta[l]-dt*2.0/(drho[l]+drho[ljm])*(dqt[l]-dqt[ljm])*invdxtheta; /* dust */
      }
    }
  }
}
		   		   
long ConditionCFL (Vrad, Vtheta, DVrad, DVtheta, deltaT)
PolarGrid *Vrad, *Vtheta, *DVrad, *DVtheta;
real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D], DVresidual[MAX1D], DVmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, invdt5, invdt6, invdt7, invdt8, invdt9, cs, csd, newdt, dt, dtg, dtd;
  int ideb=0, jdeb=0;
  real itdbg=0., itdbg1=0., itdbg2=0., itdbg3=0., itdbg4=0., itdbg5=0., itdbg6=0., itdbg7=0., itdbg8=0., itdbg9=0., mdtdbg=0.; /* debugging variables */
  real *vt, *vr, *dvt, *dvr, *ts, dxrad, dxtheta, devr, devt, dedvr, dedvt,viscr=0., visct=0.;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;

  vt = Vtheta->Field;
  vr = Vrad->Field;

  dvt = DVtheta->Field;
  dvr = DVrad->Field;

  ts = TStop->Field;

  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    DVmoy[i]= 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
      DVmoy[i]+= dvt[l];
    }
    Vmoy[i] /= (real)ns;
    DVmoy[i]/= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES){
	Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
	DVresidual[j] = dvt[l]-DVmoy[i];
      }else{
	Vresidual[j] = vt[l];	       /* Standard algorithm */
	DVresidual[j]= dvt[l];
      }
    }
    Vresidual[ns]=Vresidual[0];
    DVresidual[ns]=DVresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = SOUNDSPEED[i];
      csd = DSOUNDSPEED[i];
      invdt1 = cs/(min2(dxrad,dxtheta));

      if(!SFTA)invdt7 = csd/(min2(dxrad,dxtheta));
      else invdt7=1.e-10;

      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;

      if(!SFTA)invdt5 = fabs(dvr[l])/dxrad;
      else invdt5=1.e-10;
      if(!SFTA)invdt6 = fabs(DVresidual[j])/dxtheta;
      else invdt6=1.e-10;

      if(Rmed[i]>COUPLEINN){
	invdt8 = 1./(ts[l]+1.e-6)*sqrt(G*1.0/Rmed[i])/Rmed[i];
	if(DustImp==YES||GasDcouple==NO)invdt8 = 1.e-10;
 	if(SFTA)invdt8=1.e-10;
      }else{
        invdt8=1e-10;
      }

      devr = vr[lip]-vr[l];
      devt = vt[ljp]-vt[l];
      if (devr >= 0.0) devr = 1e-10;
      else devr = -devr;
      if (devt >= 0.0) devt = 1e-10;
      else devt = -devt;
      invdt4 = max2(devr/dxrad,devt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;

      dedvr = dvr[lip]-dvr[l];
      dedvt = dvt[ljp]-dvt[l];
      if (dedvr >= 0.0) dedvr = 1e-10;
      else dedvr = -dedvr;
      if (dedvt >= 0.0) dedvt = 1e-10;
      else dedvt = -dedvt;
      invdt9 = max2(dedvr/dxrad,dedvt/dxtheta);
      invdt9*= 4.0*CVNR*CVNR;


      dtg = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4);
      dtd = CFLSECURITY/sqrt(invdt7*invdt7+invdt5*invdt5+invdt6*invdt6+invdt8*invdt8+invdt9*invdt9);
      dt=min2(dtg,dtd);
      dt=max2(dt,1.e-10);
      if (dt < newdt) {
	newdt = dt;
	ideb = i;
        jdeb = j;
	itdbg8=1.0/invdt8;
	itdbg=dtg;
	if (debug == YES) {
	  ideb = i;
	  jdeb = j;
	  itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4; itdbg5=1.0/invdt5; itdbg6=1.0/invdt6 ; itdbg7=1.0/invdt7 ; itdbg8=1.0/invdt8; itdbg9=1.0/invdt9;
	  mdtdbg = newdt;
	  viscr = dxrad/devr/4.0/CVNR/CVNR;     
	  visct = dxtheta/devt/4.0/CVNR/CVNR;
	}
      }  
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }
  if(!SFTA){
    for (i = Zero_or_active; i < MaxMO_or_active; i++) {
      dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(DVmoy[i]*InvRmed[i]-DVmoy[i+1]*InvRmed[i+1]);
      if (dt < newdt) newdt = dt;
    }
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU SFTA %d %d: \n", CPU_Rank,SFTA);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Viscosity limit                : %g\n", itdbg4);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Dust Viscosity limit                : %g\n", itdbg9);
    printf ("Dust Sound speed limit	    : %g\n", itdbg7);
    printf ("Dust Radial motion limit            : %g\n", itdbg5);
    printf ("Dust Residual circular motion limit : %g\n", itdbg6);
    printf ("Dust Radial aerodrag limit     : dvr %g vr %g ts %g %g\n", dvr[jdeb+ideb*ns], vr[jdeb+ideb*ns], ts[jdeb+ideb*ns], itdbg8);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (long)(ceil(deltaT/newdt));
}

void AddMass(Rho,Vrad,Vtheta,dt)
     PolarGrid *Rho, *Vrad, *Vtheta;
     real dt;
{
  int i,j,l,nr,ns;
  real *rho, *vrad, *vtheta, addma, dvth;
  real tenvcgs, TUNIT, rgascgs, ricgs, omecgs, xgcontcgs, xmccgs, cscgs, xdmcgs, tcgs, r0cgs, rccgs, xdm;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  TUNIT=1./sqrt(6.67e-8/LUNIT/LUNIT/LUNIT*MUNIT);
#pragma
   for (i = 0; i < nr; i++){
    for (j = 0; j < ns; j++){
      l = j+i*ns;
        if(Rmed[i]>RIN&&Rmed[i]<ROUT){
          // Adds Dust!
          addma = 2.e-6*6.3e25/MUNIT*TUNIT*0.25/(pow(ROUT,0.25)-pow(RIN,0.25))/2./PI*pow(Rmed[i],-1.75);
          rho[l] += addma*dt;
        }
    }
   }
}

