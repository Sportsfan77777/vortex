/* 12/1/2006 * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This file is done for Barycentric computation.              *
 * It uses Gframeforce.c and Pframeforce1D.C .                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include "mp.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

static PolarGrid *TemperInt;
static PolarGrid *VradNew,   *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid1D *VradNew1D,   *VradInt1D;
static PolarGrid1D *VthetaNew1D, *VthetaInt1D;
static real timeCRASH;  
extern boolean Corotating;
real PhysicalTime=0.0, OmegaFrame;
int FirstGasStepFLAG=1;

static int AlreadyCrashed = 0, GasTimeStepsCFL;
extern boolean FastTransport, IsDisk;

extern real LostMass2D;
Pair Delta_H_Flux;



boolean DetectCrash2D1D (array,array1D)
PolarGrid *array;
PolarGrid1D *array1D;
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
  if (IAmTheFirst || IAmTheLast) {
    nr = array1D->Nrad;
    ptr= array1D->Field;
#pragma omp parallel for private() shared(bool)
    for (i = 0; i < nr; i++) {
      if (ptr[i] < 0.0)
	bool = YES;
    }
  }
  return bool;
}

 
void FillPolar1DArrays2D1D ()
{
  FILE *output1D, *output, *input1D, *input;
  int i,ii,j,quid,i_inner;
  real drrsep, test;
  float temporary;
  char InputName1D[256], InputName[256];
  char OutputName1D[256], OutputName[256];

  ReadPrevDim();
  ReadPrevDim1D();

  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (InputName1D, "%s%s", OUTPUTDIR, "radii1D.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  sprintf (OutputName1D, "%s%s", OUTPUTDIR, "used_rad1D.dat");
  input = fopen (InputName, "r");
  input1D = fopen (InputName1D, "r");
  test = 0.;

  if ( (input == NULL) || (input1D == NULL) ) {

    mastererr ("Warning : no `radii.dat' or 'radii1D.dat' file found. Using default.\n");
    if ( (RMIN1D >= RMIN) && (RMAX1D <= RMAX) )
      erreur("\nERROR: 1D-Grid is included into 2D-Grid. Absurd.\n");
    if (LogGrid == YES) {
      if (RMIN1D<RMIN) {
	i_inner = - log(RMIN1D/RMIN)/log(RMAX/RMIN)*GLOBALNRAD;
	if (RMAX1D>RMAX)
	  NRAD1D = log(RMAX1D/RMIN1D)/log(RMAX/RMIN)*GLOBALNRAD;
	else
	  NRAD1D = i_inner+2*CPUOVERLAP;
      } else {
	i_inner = -(NRAD-2*CPUOVERLAP);
	NRAD1D = log(RMAX1D/RMIN)/log(RMAX/RMIN)*GLOBALNRAD+i_inner;
      }
      if(NRAD1D>MAX1D) erreur("1D-grid too large. Change its size or MAX1D. Exit.\n");
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
      for (i = 0; i <= NRAD1D; i++) {
	Radii1D[i] = RMIN*exp((real)(i-i_inner)/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      if (RMIN1D<RMIN) {
	i_inner = (RMIN-RMIN1D)/drrsep;
	if (RMAX1D>RMAX) NRAD1D = (RMAX1D-RMIN1D)/drrsep;
	else NRAD1D = i_inner + 2*CPUOVERLAP;
      } else {
	i_inner = -(NRAD-2*CPUOVERLAP);
	NRAD1D = (RMAX1D-RMAX)/drrsep + 2*CPUOVERLAP;
      }
      if(NRAD1D>MAX1D) erreur("1D-grid too large. Change its size or MAX1D. Exit.\n");
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
      for (i = 0; i <= NRAD1D; i++) {
	Radii1D[i] = RMIN+drrsep*(real)(i-i_inner);
      }
    }

  } else {

    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
    mastererr ("Reading 'radii1D.dat' file.\n");
    i = 0;
    quid = 1;
    fscanf (input1D, "%f", &temporary);
    while (quid > 0 && i <= MAX1D) {
      Radii1D[i] = (real)temporary;
      quid = fscanf (input1D, "%f", &temporary);
      i++;
    }
    NRAD1D = i-1;
    i = 0;
    if (Radii1D[0] < Radii[0]) {
      while (Radii1D[i] < 0.9999999*Radii[0]) {
	i++;
      }
      i_inner=i;
      for (i = 0; i < min2(((real)GLOBALNRAD),(real)(NRAD1D-i_inner)); i++) {
	test = max2(test,fabs(Radii[i]-Radii1D[i_inner+i])/Radii[i]);
      }
    } else {
      while (Radii[i]<0.9999999*Radii1D[0]) {
	i++;
      }
      i_inner = -i;
      for (i = -i_inner; i<GLOBALNRAD; i++)
	test = max2(test,fabs(Radii[i]-Radii1D[i_inner+i])/Radii[i]);
    }
    if (fabs(test) > 1.e-8)
      erreur("The two files are not consistent. Exit.\n");
  }
  RMIN1D=Radii1D[0];
  RMAX1D=Radii1D[NRAD1D];
  IINNER=i_inner;
  masterprint("Warning : Corrected values of 1D-Grid size :\n");
  masterprint("          RMIN1D = %f, RMAX1D = %f ; IINNER = %d, NRAD1D = %d\n",RMIN1D,Radii1D[NRAD1D],IINNER,NRAD1D);
  if(IINNER>=0 && IINNER<3)
    masterprint("Warning : 2D-Grid overlaps 1D-Grid 3 first rings.\n          Risk of instabilities.\n");
  i=NRAD1D-IINNER-GLOBALNRAD;
  if(i>=0 && i<3)
    masterprint("Warning : 2D-Grid overlaps 1D-Grid 3 last rings.\n          Risk of instabilities.\n");

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
  for (i = 0; i < NRAD1D; i++) {
    Rinf1D[i] = Radii1D[i];
    Rsup1D[i] = Radii1D[i+1];
    Rmed1D[i] = 2.0/3.0*(Rsup1D[i]*Rsup1D[i]*Rsup1D[i]-Rinf1D[i]*Rinf1D[i]*Rinf1D[i]);
    Rmed1D[i] = Rmed1D[i] / (Rsup1D[i]*Rsup1D[i]-Rinf1D[i]*Rinf1D[i]);
    Surf1D[i] = PI*(Rsup1D[i]*Rsup1D[i]-Rinf1D[i]*Rinf1D[i]);
    InvRmed1D[i] = 1.0/Rmed1D[i];
    InvSurf1D[i] = 1.0/Surf1D[i];
    InvDiffRsup1D[i] = 1.0/(Rsup1D[i]-Rinf1D[i]);
    InvRinf1D[i] = 1.0/Rinf1D[i];
  }
  Rinf1D[NRAD1D]=Radii1D[NRAD1D];

  for (j = 0; j < NSEC; j++){
    Cosinus[j]=cos(2.0*PI*(real)j/(real)NSEC);
    Sinus[j] = sin(2.0*PI*(real)j/(real)NSEC);
  }

  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  for (i = 1; i < NRAD1D; i++) {
    InvDiffRmed1D[i] = 1.0/(Rmed1D[i]-Rmed1D[i-1]);
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
    output1D = fopen (OutputName1D, "w");
    if (output1D == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName1D);
      prs_exit (1);
    }
    for (i = 0; i <= NRAD1D; i++) {
      fprintf (output1D, "%.18g\n", Radii1D[i]);
    }
    fclose (output1D);
  }
  if (input != NULL) fclose (input);
  if (input1D != NULL) fclose (input1D);
}


void InitEuler2D1D (Rho, Vr, Vt, Rho1D, Vr1D, Vt1D)
PolarGrid *Rho, *Vr, *Vt;
PolarGrid1D *Rho1D, *Vr1D, *Vt1D;
{
  masterprint("InitEuler2D1D... ");
  fflush(stdout);
  FillSigma ();
  FillSigma1D ();
  masterprint("(FillSigma done. ");
  fflush(stdout);
  InitTransport ();
  InitTransport1D ();
  InitViscosity ();
  InitViscosity1D ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  Potplanet    = CreatePolarGrid(NRAD, NSEC, "Potplanet");
  RhoStar1D     = CreatePolarGrid1D(NRAD1D, "RhoStar1D");
  RhoInt1D      = CreatePolarGrid1D(NRAD1D, "RhoInt1D");
  VradNew1D     = CreatePolarGrid1D(NRAD1D, "VradNew1D");
  VradInt1D     = CreatePolarGrid1D(NRAD1D, "VradInt1D");
  VthetaNew1D   = CreatePolarGrid1D(NRAD1D, "VthetaNew1D");
  VthetaInt1D   = CreatePolarGrid1D(NRAD1D, "VthetaInt1D");
  Potential1D   = CreatePolarGrid1D(NRAD1D, "Potential1D");
  masterprint("PolarGrids creation done.) ");
  fflush(stdout);
  InitGas (Rho, Vr, Vt);
  InitGas1D (Rho1D, Vr1D, Vt1D);
  masterprint("done.\n");
  fflush(stdout);
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

void ActualiseGas1D (array, newarray)
PolarGrid1D *array, *newarray;
{
  int i,nr;
  real *old, *new;
  nr = array->Nrad;
  old= array->Field;
  new= newarray->Field;
  for (i = 0; i <= nr; i++) {
    old[i] = new[i];
  }
}


void AlgoGas2D1D (Rho,Vrad,Vtheta,Label,Rho1D,Vrad1D,Vtheta1D,Label1D,sys)
PolarGrid *Rho, *Vrad, *Vtheta, *Label;
PolarGrid1D *Rho1D, *Vrad1D, *Vtheta1D, *Label1D;
PlanetarySystem *sys;
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  int gastimestepcfl, gastimestepcfl1D;
  boolean Crashed=NO;
  FirstGasStepFLAG=1;
  gastimestepcfl = gastimestepcfl1D = 1;
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho,Vrad,Vtheta,Label);
    if (SloppyCFL == YES) {
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
      if (IAmTheFirst || IAmTheLast)
	gastimestepcfl = (int)(max2((real) gastimestepcfl, (real)ConditionCFL1D (Vrad1D, Vtheta1D, DT-dtemp) ));
      }
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;

  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho,Vrad,Vtheta,Label);
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
      if (IAmTheFirst || IAmTheLast) {
	  gastimestepcfl1D = ConditionCFL1D (Vrad1D, Vtheta1D, DT-dtemp);
	  gastimestepcfl = (int)(max2( (real)gastimestepcfl, (real)gastimestepcfl1D));
	}
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }

      /* Order of the procedures changed for the radial velocities
         calculated in the ghosts to correspond to the real mass flux.  */
      Delta_H_Flux.x = Delta_H_Flux.y = 0.;
      // TO TEST THE 2D GRID ALONE, COMMENT FOLLOWING LINE
      Communicate1D2D (Rho1D,Vrad1D,Vtheta1D,Label1D,Rho,Vrad,Vtheta,Label,dt);
      Transport (Rho, Vrad, Vtheta, Label, dt);
      ApplyBoundaryCondition (Vrad, Rho);
      if (IAmTheFirst || IAmTheLast) {
	Transport1D (Rho1D, Vrad1D, Vtheta1D, Label1D, dt);
	ApplyBoundaryCondition1D (Vrad1D, Rho1D);
      }
      BarycenterConservation(sys,Rho,Vrad,Vtheta);
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    dtemp += dt;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      FillForcesArrays (sys, Rho);
      if (IAmTheFirst || IAmTheLast)
	FillForcesArrays1D (sys);
      AdvanceSystemFromDisk (Rho, sys, dt);
    }
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES)
	CorrectVtheta (Vtheta, domega);
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    if (IsDisk == YES) {
      ApplyBoundaryCondition(Vrad, Rho);
      if (IAmTheFirst || IAmTheLast)
	ApplyBoundaryCondition1D(Vrad1D, Rho1D);
      Crashed = DetectCrash2D1D (Rho,Rho1D);  /* test for negative density values */
      if (Crashed == YES) {
	if (AlreadyCrashed == 0) {
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
	  WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
	  WriteDiskPolar (Vtheta, 999); /* of what happened */
	  WriteDiskPolar1D (Rho1D, 999);    /* We write the HD arrays */
	  WriteDiskPolar1D (Vrad1D, 999);   /* in order to keep a track */
	  WriteDiskPolar1D (Vtheta1D, 999); /* of what happened */
	  Write1DSummary(Rho,Rho1D,999);
	  Write1DWeightedSummary(Vrad,Rho,Vrad1D,999);
	  Write1DWeightedSummary(Vtheta,Rho,Vtheta1D,999);
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (".");
      }
      fflush (stdout);

      SubStep1 (Vrad, Vtheta, Rho, dt);
      SubStep2 (Rho, dt);
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      ApplyBoundaryCondition (Vrad, Rho);
      fflush (stdout);

      if (IAmTheFirst || IAmTheLast) {
	SubStep1_1D (Vrad1D, Vtheta1D, Rho1D, dt, sys);
      }
      AdvanceSystemFromDisk1D (Rho1D, sys, dt);
      if (IAmTheFirst || IAmTheLast) {
	SubStep2_1D (Rho1D, dt);
	ActualiseGas1D (Vrad1D, VradNew1D);
	ActualiseGas1D (Vtheta1D, VthetaNew1D);
	ApplyBoundaryCondition1D (Vrad1D, Rho1D);
     }

    }
    AdvanceSystemRK5 (sys, dt);
    PhysicalTime += dt;
  }
  masterprint ("\n");
}


void SubStep1 (Vrad, Vtheta, Rho, dt)
PolarGrid *Vrad, *Vtheta, *Rho;
real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta, cs2, cs2m;
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
  Pot = Potential->Field;
				/* In this substep we take into account     */
				/* the source part of Euler equations       */
				/* (i.e. the R.H.S. in classical notation). */
#pragma omp parallel private(j,l,lim,ljm,ljp,cs2,cs2m,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
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
	vthetaint[l] += dt*supp_torque;

      }
    }
  }
  ViscousTerms (VradInt, VthetaInt, Rho, dt);
}


void SubStep1_1D (Vrad, Vtheta, Rho, dt, sys)
PolarGrid1D *Vrad, *Vtheta, *Rho;
real dt;
PlanetarySystem *sys;
{
  int i, nr, k;
  real *vrad, *vtheta, *rho;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, cs2, cs2m;
  real supp_torque=0.0;		/* for imposed disk drift */
  real *Pot;
  real xplanet,yplanet,mplanet,rplanet,xstar,ystar,RRoche,tg;
  real PlanetDistance,Delta,Delta4,e;
  nr = Vrad->Nrad;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt1D->Field;
  vthetaint = VthetaInt1D->Field;
  Pot = Potential1D->Field;
				/* In this substep we take into account     */
				/* the source part of Euler equations       */
				/* (i.e. the R.H.S. in classical notation). */
#pragma omp parallel private(cs2,cs2m,vt2,gradp,gradphi,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      cs2 = SOUNDSPEED1D[i]*SOUNDSPEED1D[i];
      cs2m= SOUNDSPEED1D[i-1]*SOUNDSPEED1D[i-1];
      gradp = (cs2*rho[i]-cs2m*rho[i-1])*2.0/(rho[i]+rho[i-1])*InvDiffRmed1D[i];
      gradphi = (Pot[i]-Pot[i-1])*InvDiffRmed1D[i];
      vt2 = vtheta[i]+vtheta[i-1];
      vt2 = vt2/2.0+Rinf1D[i]*OmegaFrame;
      vt2 = vt2*vt2;
      vradint[i] = vrad[i]+dt*(-gradp-gradphi+vt2*InvRinf1D[i]);
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*.5*pow(Rmed1D[i],-2.5+SIGMASLOPE);
      vthetaint[i] = vtheta[i];
      vthetaint[i] += dt*supp_torque;
      /* Planetary torque tg added on 1D rings : */
      for (k = 1; k < sys->nb; k++) {
	if (sys->FeelDisk[k] == YES) {
	  xplanet = sys->x[k];
	  yplanet = sys->y[k];
	  mplanet = sys->mass[k];
	  rplanet = sqrt(xplanet*xplanet+yplanet*yplanet);
	  xstar=sys->x[0];
	  ystar=sys->y[0];
	  PlanetDistance = sqrt( (xplanet-xstar)*(xplanet-xstar)+\
				 (yplanet-ystar)*(yplanet-ystar) );
	  RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
	  Delta = Rmed1D[i] - rplanet;
	  if ( fabs(Delta) > 2.*RRoche ) {
	    e = Delta / fabs(Delta);
	    Delta4 = Delta*Delta*Delta*Delta;
	    tg = 0.35*mplanet*mplanet*rplanet*rplanet*Rmed1D[i]/Delta4*e;
	    vthetaint[i] += dt*tg/Rmed1D[i];
	  }
	}
      }
    }
  }
  ViscousTerms1D (VradInt1D, VthetaInt1D, Rho, dt);
}


void SubStep2 (Rho, dt)
PolarGrid *Rho;
real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho;
  real *vradnew, *vthetanew, *qt, *qr;
  real dxtheta, invdxtheta;
  real dv;
  nr = VradInt->Nrad;
  ns = VradInt->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
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
	lip = l+ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
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
      }
    }
  }
}

void SubStep2_1D (Rho, dt)
PolarGrid1D *Rho;
real dt;
{
  int i, nr;
  real *vrad, *vtheta, *rho;
  real *vradnew, *vthetanew, *qr;
  real dv;
  nr = VradInt1D->Nrad;
  rho = Rho->Field;
  vrad = VradInt1D->Field;
  qr = RhoInt1D->Field;
  vradnew = VradNew1D->Field;
  vtheta = VthetaInt1D->Field;
  vthetanew = VthetaNew1D->Field;
#pragma omp parallel for private(dv)
  for (i = 0; i < nr; i++) {
    dv = vrad[i+1]-vrad[i];
    if (dv < 0.0)
      qr[i] = CVNR*CVNR*rho[i]*dv*dv;
    else 
      qr[i] = 0.0; 
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      vradnew[i] = vrad[i]-dt*2.0/(rho[i]+rho[i-1])*(qr[i]-qr[i-1])*InvDiffRmed1D[i];
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      vthetanew[i] = vtheta[i];
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
      cs = SOUNDSPEED[i];
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
	   		   
int ConditionCFL1D (Vrad, Vtheta, deltaT)
PolarGrid1D *Vrad, *Vtheta;
real deltaT;
{
  static real Vresidual, Vmoy[MAX1D];
  int i, nr, MinLoop, MaxLoop;
  real invdt1, invdt2, invdt4, cs, newdt, dt;
  int ideb;
  real itdbg1, itdbg2, itdbg4, mdtdbg; /* debugging variables */
  real *vt, *vr, dxrad, dvr, viscr;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    Vmoy[i] = vt[i];
  }
  MinLoop = ( IAmTheFirst ? 0 : IINNER-CPUOVERLAP+GLOBALNRAD);
  MaxLoop = ( IAmTheLast ? NRAD1D : IINNER+CPUOVERLAP);	/* Correct even in sequential version */
  for (i = MinLoop; i < MaxLoop; i++) {
    dxrad = Rsup1D[i]-Rinf1D[i];
    cs = SOUNDSPEED1D[i];
    invdt1 = cs/dxrad;
    invdt2 = fabs(vr[i])/dxrad;
    dvr = vr[i+1]-vr[i];
    if (dvr >= 0.0) dvr = 1e-10;
    else dvr = -dvr;
    invdt4 = max2(dvr/dxrad,0.);
    invdt4*= 4.0*CVNR*CVNR;
    dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt4*invdt4);
    if (dt < newdt) {
      newdt = dt;
      if (debug == YES) {
	ideb = i;
	itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg4=1.0/invdt4;
	mdtdbg = newdt;
	viscr = dxrad/dvr/4.0/CVNR/CVNR;     
      }
    }  
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("  ConditionCFL1D :\n");
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d in 1D-grid\n", ideb);
    printf ("located at radius Rmed         : %g\n", Rmed1D[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : none\n");
    printf ("Viscosity limit                : %g\n", itdbg4);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg)
      printf ("Discrepancy arise from the imposed DT interval.\n");
  }
  return (int)(ceil(deltaT/newdt));
}
