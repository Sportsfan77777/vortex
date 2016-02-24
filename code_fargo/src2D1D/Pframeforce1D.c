#include "mp.h"

extern boolean AllowAccretion, Corotating;
static real vt_int_1D[MAX1D], vt_cent1D[MAX1D];


				/* Below : work in non-rotating frame */
				/* centered on the primary */
void FillForcesArrays1D (sys)
PlanetarySystem *sys;
{
  int i,nr,k,NbPlanets;
  real *Pot;
  real InvDistance, Xp, Yp, mplanet, PlanetDistance;
  real massPlanets, mstar;
  Pot= Potential1D->Field;
  nr = Potential1D->Nrad;
  NbPlanets = sys->nb;
  mstar = sys->mass[0];
#pragma omp parallel for private(InvDistance)
  for (i = 0; i < nr; i++) {
    massPlanets=0.0000000000;
#pragma omp parallel for
    for (k = 1; k < NbPlanets; k++) {
      Xp = sys->x[k];
      Yp = sys->y[k];
      mplanet = sys->mass[k];
      if (k>0) mplanet *= MassTaper;
      PlanetDistance = sqrt(Xp*Xp+Yp*Yp);
      if(Rinf1D[i]>PlanetDistance)
	massPlanets+= mplanet;
      else if (Rsup1D[i]>PlanetDistance)
	massPlanets+= 0.5*mplanet;
    }
    InvDistance = 1.0/Rmed1D[i];
    Pot[i] = -G*(mstar+massPlanets)*InvDistance;
  }
} 

void AdvanceSystemFromDisk1D (Rho, sys, dt)
PolarGrid1D *Rho;
PlanetarySystem *sys;
real dt;		       
{
  int NbPlanets, k, i, nr;
  real *rho;
  real x,y,rp,m,dvtp,thetap;
  real xstar,ystar,PlanetDistance,RRoche;
  real tg,dhp,dHp,Delta,Delta4,e;
  NbPlanets = sys->nb;
  rho = Rho->Field;
  nr = Rho->Nrad;
  for (k = 1; k < NbPlanets; k++) {
    if (sys->FeelDisk[k] == YES) {
      m=sys->mass[k];
      x=sys->x[k];
      y=sys->y[k];
      rp=sqrt(x*x+y*y);
      xstar=sys->x[0];
      ystar=sys->y[0];
      PlanetDistance = sqrt( (x-xstar)*(x-xstar)+\
			     (y-ystar)*(y-ystar) );
      RRoche = PlanetDistance*pow((1.0/3.0*m),1.0/3.0);
      dhp = 0.;

      if (IAmTheFirst && (RMIN1D<RMIN)) {
	for (i = 1; i < IINNER+CPUOVERLAP; i++) {
	  Delta = Rmed1D[i] - rp;
	  if ( fabs(Delta) > 2.*RRoche ) {
	    e = Delta / fabs(Delta);
	    Delta4 = Delta*Delta*Delta*Delta;
	    tg = 0.35*m*m*rp*rp*Rmed1D[i]/Delta4*e;
	    dhp -= tg*rho[i]*Surf1D[i];
	  }
	}
      }
      if (IAmTheLast && (RMAX1D>RMAX)) {
	for (i = IINNER+GLOBALNRAD-CPUOVERLAP; i < NRAD1D-1; i++) {
	  Delta = Rmed1D[i] - rp;
	  if ( fabs(Delta) > 2.*RRoche ) {
	    e = Delta / fabs(Delta);
	    Delta4 = Delta*Delta*Delta*Delta;
	    tg = 0.35*m*m*rp*rp*Rmed1D[i]/Delta4*e;
	    dhp -= tg*rho[i]*Surf1D[i];
	  }
	}
      }
      MPI_Allreduce (&dhp,&dHp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      thetap = atan2(y,x);
      dvtp = dHp/m/rp*dt;
      sys->vx[k] += - dvtp*sin(thetap);
      sys->vy[k] +=   dvtp*cos(thetap);
    }
  }
}


void InitGas1D (Rho, Vr, Vt)
PolarGrid1D *Rho, *Vr, *Vt;
{
  int i, nr;
  real *dens, *vr, *vt;
  float temporary;
  FILE *CS;
  char csfile[512];
  real  r, rg, omega, ri, viscosity, t1, t2, r1, r2;
  dens= Rho->Field;
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Rho->Nrad;
  sprintf (csfile, "%s%s", OUTPUTDIR, "soundspeed.dat");
  CS = fopen (csfile, "r");
  if (CS == NULL) {
    for (i = 0; i < nr; i++) {
      SOUNDSPEED1D[i] = ASPECTRATIO * sqrt(G*1.0/Rmed1D[i]) * pow(Rmed1D[i], FLARINGINDEX);
    }
  } else {
    masterprint ("Reading soundspeed.dat file\n");
    for (i = 0; i < NRAD1D; i++) {
      fscanf (CS, "%f", &temporary);
      SOUNDSPEED1D[i] = (real)temporary;
    }
  }
  for (i = 1; i < NRAD1D; i++) {
    vt_int_1D[i]=(SOUNDSPEED1D[i]*SOUNDSPEED1D[i]*Sigma(Rmed1D[i])-\
	       SOUNDSPEED1D[i-1]*SOUNDSPEED1D[i-1]*Sigma(Rmed1D[i-1]))/\
      (.5*(Sigma(Rmed1D[i])+Sigma(Rmed1D[i-1])))/(Rmed1D[i]-Rmed1D[i-1])+\
      G*(1.0/Rmed1D[i-1]-1.0/Rmed1D[i])/(Rmed1D[i]-Rmed1D[i-1]);
    vt_int_1D[i] = sqrt(vt_int_1D[i]*Radii1D[i])-Radii1D[i]*OmegaFrame;
  }
  t1 = vt_cent1D[0] = vt_int_1D[1]+.75*(vt_int_1D[1]-vt_int_1D[2]);
  r1 = ConstructSequence (vt_cent1D, vt_int_1D, NRAD1D);
  vt_cent1D[0] += .25*(vt_int_1D[1]-vt_int_1D[2]);
  t2 = vt_cent1D[0];
  r2 = ConstructSequence (vt_cent1D, vt_int_1D, NRAD1D);
  t1 = t1-r1/(r2-r1)*(t2-t1);
  vt_cent1D[0] = t1;
  ConstructSequence (vt_cent1D, vt_int_1D, NRAD1D);
  vt_cent1D[NRAD1D]=vt_cent1D[NRAD1D-1];

  for (i = 0; i <= nr; i++) {
    if (i == nr) {
      r = Rmed1D[nr-1];
      ri= Rinf1D[nr-1];
    }
    else {
      r = Rmed1D[i];
      ri= Rinf1D[i];
    }

    viscosity = VISCOSITY;
    if (ViscosityAlpha)
      viscosity = ALPHAVISCOSITY*ASPECTRATIO*ASPECTRATIO*pow(r,.5+2.0*FLARINGINDEX);

    rg = r;
    omega = sqrt(G*1.0/rg/rg/rg);
    vt[i] = omega*r*\
      sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	   pow(r,2.0*FLARINGINDEX)*\
	   (1.+SIGMASLOPE-2.0*FLARINGINDEX));
    vt[i] -= OmegaFrame*r;
    if (CentrifugalBalance)
      vt[i] = vt_cent1D[i];
    if (i == nr) 
      vr[i] = 0.0;
    else {
      vr[i] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf1D[i]/ri;
      if (ViscosityAlpha) {
	vr[i] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
      } else {
	vr[i] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
      }
    }
    dens[i] = SigmaMed1D[i];
  }
  vr[0] = vr[nr] = 0.0;
}
