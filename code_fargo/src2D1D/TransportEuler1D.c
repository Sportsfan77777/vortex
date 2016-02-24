#include "mp.h"

extern real OmegaFrame;

static real *dq1D;
static PolarGrid1D *RadMomP1D, *RadMomM1D, *ThetaMom1D, *VthetaRes1D;
static PolarGrid1D *Work1D, *QRStar1D, *ExtLabel1D, *Elongations1D;

extern int TimeStep;
extern boolean OpenInner, FastTransport;
extern real LostMass, LostLabel;
extern Pair Delta_H_Flux;
static boolean H = NO;

void Transport1D (Rho, Vrad, Vtheta, Label, dt)
PolarGrid1D *Rho, *Vrad, *Vtheta, *Label;
real dt;
{
  ComputeLRMomenta1D (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES)
     ComputeExtQty1D (Rho, Label, ExtLabel1D);
				/* No-Alternate Directionnal Splitting */
  OneWindRad1D (Rho, Vrad, dt);
  ComputeVelocities1D (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES) 
     ComputeSpeQty1D (Rho, Label, ExtLabel1D);
}


void OneWindRad1D (Rho, Vrad, dt)
PolarGrid1D *Rho, *Vrad;
real dt;
{
  real *vr;
  vr = Vrad->Field;
  H = NO;
  ComputeStarRad1D (Rho, Vrad, RhoStar1D, dt);
  ActualiseGas1D (RhoInt1D, Rho);
  VanLeerRadial1D (Vrad, RadMomP1D, dt); 
  VanLeerRadial1D (Vrad, RadMomM1D, dt);  
  H = YES;
  VanLeerRadial1D (Vrad, ThetaMom1D, dt);   
  Tartinage (ThetaMom1D,dt,Delta_H_Flux);
  H = NO;
  if (AdvecteLabel == YES)
    LostLabel += VanLeerRadial1D (Vrad, ExtLabel1D, dt);
  LostMass += VanLeerRadial1D (Vrad, Rho, dt); /* MUST be the last line */
}


   /* Hereafter are the new specific procedures to the fast algorithm */

void ComputeExtQty1D (Rho, Label, ExtLabel1D)
PolarGrid1D *Rho, *Label, *ExtLabel1D;
{
  int i,nr;
  real *rho, *lab, *extlab;
  nr = Rho->Nrad;
  rho= Rho->Field;
  lab= Label->Field;
  extlab= ExtLabel1D->Field;
  for (i = 0; i < nr; i++) {
    extlab[i] = rho[i]*lab[i]; /* compressive flow  if line commentarized
				    extlab[l] = lab[l]; */
  }
}

void ComputeSpeQty1D (Rho, Label, ExtLabel1D)
PolarGrid1D *Rho, *Label, *ExtLabel1D;
{
  int i,nr;
  real *rho, *lab, *extlab;
  nr = Rho->Nrad;
  rho= Rho->Field;
  lab= Label->Field;
  extlab= ExtLabel1D->Field;
  for (i = 0; i < nr; i++) {
    lab[i] = extlab[i]/rho[i]; /* compressive flow if line commentarized
				    lab[l] = extlab[l]; */
  }
}

void InitTransport1D () 
{
  RadMomP1D    = CreatePolarGrid1D(NRAD1D, "RadMomP1D");
  RadMomM1D    = CreatePolarGrid1D(NRAD1D, "RadMomM1D");
  ThetaMom1D   = CreatePolarGrid1D(NRAD1D, "ThetaMom1D");
  Work1D       = CreatePolarGrid1D(NRAD1D, "Work1DGrid");
  QRStar1D     = CreatePolarGrid1D(NRAD1D, "QRStar1D");
  ExtLabel1D   = CreatePolarGrid1D(NRAD1D, "ExtLabel1D");
  VthetaRes1D  = CreatePolarGrid1D(NRAD1D, "VthetaRes1D");
  Elongations1D= CreatePolarGrid1D(NRAD1D, "Elongations1D");
  dq1D         = (real *)malloc(NRAD1D*sizeof(real));
}

void ComputeStarRad1D (Qbase, Vrad, QStar, dt)
PolarGrid1D *Qbase, *Vrad, *QStar;
real dt;
{
  int i,nr;
  real *qb, *qs, *vr;
  real dqp, dqm;
  nr = Qbase->Nrad;
  qb = Qbase->Field;
  qs = QStar->Field;
  vr = Vrad->Field;
#pragma omp parallel for private (i,dqm,dqp)   

  for (i = 0; i < nr; i++) {
    if ((i == 0) || (i == nr-1)) dq1D[i] = 0.0;
    else {
      dqm = (qb[i]-qb[i-1])*InvDiffRmed1D[i];
      dqp = (qb[i+1]-qb[i])*InvDiffRmed1D[i+1];
      if (dqp * dqm > 0.0)
	dq1D[i] = 2.0*dqp*dqm/(dqp+dqm);
      else
	dq1D[i] = 0.0;
    }
  }
  for (i = 0; i < nr; i++) {
    if (vr[i] > 0.0)
      qs[i] = qb[i-1]+(Rmed1D[i]-Rmed1D[i-1]-vr[i]*dt)*0.5*dq1D[i-1];
    else
      qs[i] = qb[i]-(Rmed1D[i+1]-Rmed1D[i]+vr[i]*dt)*0.5*dq1D[i];
  }
  qs[0] = qs[nr] = 0.0;
}


void ComputeLRMomenta1D (Rho, Vrad, Vtheta)
PolarGrid1D *Rho, *Vrad, *Vtheta;
{
  int i,nr;
  real *vr,*vt,*rho;
  real *rp, *rm, *t;
  nr = Vrad->Nrad;
  rho= Rho->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  rp = RadMomP1D->Field;
  rm = RadMomM1D->Field;
  t  = ThetaMom1D->Field;
#pragma omp parallel for private()
  for (i = 0; i < nr; i++) {
    rp[i] = rho[i]*vr[i+1];
    rm[i] = rho[i]*vr[i];
    t[i]  = (rho[i])*(vt[i]+Rmed1D[i]*OmegaFrame)*Rmed1D[i]; /* it is the angular momentum */
  }
}


void ComputeVelocities1D (Rho, Vrad, Vtheta)
PolarGrid1D *Rho, *Vrad, *Vtheta;
{
  int i,nr;
  real *vr,*vt,*rho;
  real *rp, *rm, *t;
  nr = Vrad->Nrad;
  rho= Rho->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  rp = RadMomP1D->Field;
  rm = RadMomM1D->Field;
  t  = ThetaMom1D->Field;
#pragma omp parallel for private()
  for (i = 0; i < nr; i++) {
    if (i == 0)
      vr[i] = 0.0;
    else
      vr[i] = (rp[i-1]*(Rinf1D[i]-Rmed1D[i-1])+rm[i]*(Rmed1D[i]-Rinf1D[i]))/(rho[i]*(Rmed1D[i]-Rinf1D[i])+rho[i-1]*(Rinf1D[i]-Rmed1D[i-1]));
    vt[i] = (t[i])/(rho[i])/Rmed1D[i]-Rmed1D[i]*OmegaFrame;
    /* It was the angular momentum */
  }
}


real VanLeerRadial1D (Vrad, Qbase, dt)
PolarGrid1D *Vrad, *Qbase;
real dt;
{
  int i,nr;
  MPI_Request req;
  real *qrs, *rhos, *vr, *qb;
  real dtheta, varq;
  real lostByDisk=0.0, LostByDisk;
  DivisePolarGrid1D (Qbase, RhoInt1D, Work1D);
  ComputeStarRad1D (Work1D, Vrad, QRStar1D, dt);
  nr = Qbase->Nrad;
  qrs  = QRStar1D->Field;
  rhos = RhoStar1D->Field;
  vr = Vrad->Field;
  qb = Qbase->Field;
  dtheta = 2.0*PI;
  if(vr[nr]!=0.) fprintf(stderr,"Error : vr[NRAD1D]=%f",vr[nr]);
#pragma omp parallel for private(varq)
  for (i = 0; i < nr; i++) {
    varq =dt*dtheta*Rinf1D[i]*qrs[i]*rhos[i]*vr[i];
    varq-=dt*dtheta*Rsup1D[i]*qrs[i+1]*rhos[i+1]*vr[i+1];
    qb[i] += varq*InvSurf1D[i];
    if ( ((i==0)&&(RMIN1D<=RMIN)&&(IAmTheFirst)) ||\
	 ((i==nr-1)&&(RMAX1D>=RMAX)&&(IAmTheLast)) )
#pragma omp atomic
      lostByDisk += varq;
  }
  LostByDisk = lostByDisk;
  if (CPU_Number > 1) {
    if (IAmTheLast) {
      MPI_Isend (&lostByDisk,1,MPI_DOUBLE,0,34,MPI_COMM_WORLD,&req);
      MPI_Wait (&req, &stat);
    }
    if (IAmTheFirst) {
      MPI_Irecv (&lostByDisk,1,MPI_DOUBLE,CPU_Number-1,34,MPI_COMM_WORLD,&req);
      MPI_Wait (&req, &stat);
    }
    LostByDisk += lostByDisk;
  }
  if (H == YES) {
    if (IAmTheFirst) {
      i=IINNER+CPUOVERLAP-1;
      Delta_H_Flux.y += dt*dtheta*Rsup1D[i]*qrs[i+1]*rhos[i+1]*vr[i+1];
    }
    if (IAmTheLast) {
      i=IINNER+GLOBALNRAD-CPUOVERLAP;
      Delta_H_Flux.x -= dt*dtheta*Rinf1D[i]*qrs[i]*rhos[i]*vr[i];
    }
  }
  return LostByDisk;
}

/* The following procedure spreads into the 1D grid the 
 * angular momentum flux through the boundaries of the 2D grid.
 */
void Tartinage (Qbase, dt, Flux)
PolarGrid1D *Qbase;
real dt;
Pair Flux;
{
  int i,nr;
  real *qb;
  real R, l, varq;
  nr = Qbase->Nrad;
  qb = Qbase->Field;
  l = 0.5; /* Damping length of the wave. */

#pragma omp parallel for private(varq)
  for (i = 0; i < nr; i++) {
    varq = 0.;
    if ( (i<IINNER+CPUOVERLAP) && (Flux.y!=0.) && IAmTheFirst ) {
      R = Rinf1D[IINNER+CPUOVERLAP];
      varq += Flux.y*exp(-2.*(R-Rsup1D[i])/l);
      varq -= Flux.y*exp(-2.*(R-Rinf1D[i])/l);
    }
    if ((i>=IINNER+GLOBALNRAD-CPUOVERLAP) && (Flux.x!=0.) && IAmTheLast ) {
      R = Rinf1D[IINNER+GLOBALNRAD-CPUOVERLAP];
      varq += Flux.x*exp(2.*(R-Rinf1D[i])/l);
      varq -= Flux.x*exp(2.*(R-Rsup1D[i])/l);
    }
    qb[i] += varq*InvSurf1D[i];
  }
}
