#include "mp.h"

static PolarGrid1D *DivergenceVelocity1D;
static PolarGrid1D *DRR1D, *DRP1D, *DPP1D, *TAURR1D, *TAUPP1D, *TAURP1D;
extern real Trp_in,Trp_out;


void InitViscosity1D ()
{
  DivergenceVelocity1D = CreatePolarGrid1D(NRAD1D, "DivV1D");
  DRR1D                = CreatePolarGrid1D(NRAD1D, "Drr1D");
  DRP1D                = CreatePolarGrid1D(NRAD1D, "Drp1D");
  DPP1D                = CreatePolarGrid1D(NRAD1D, "Dpp1D");
  TAURR1D              = CreatePolarGrid1D(NRAD1D, "TAUrr1D");
  TAURP1D              = CreatePolarGrid1D(NRAD1D, "TAUrp1D");
  TAUPP1D              = CreatePolarGrid1D(NRAD1D, "TAUpp1D");
}

void ViscousTerms1D (RadialVelocity, AzimuthalVelocity, Rho, DeltaT)
PolarGrid1D *RadialVelocity, *AzimuthalVelocity, *Rho;
real DeltaT;
{
  int i,nr;
  real *rho, *vr, *vt;
  real *Drr, *Drp, *Dpp, *divergence;
  real *Trr, *Trp, *Tpp;
  real VKepIn, VKepOut, onethird;
  real viscosity;
  nr  = Rho->Nrad;
  rho = Rho->Field;
  vr  = RadialVelocity->Field;
  vt  = AzimuthalVelocity->Field;
  divergence = DivergenceVelocity1D->Field;
  Drr = DRR1D->Field;
  Drp = DRP1D->Field;
  Dpp = DPP1D->Field;
  Trr = TAURR1D->Field;
  Trp = TAURP1D->Field;
  Tpp = TAUPP1D->Field;
  onethird = 1.0/3.0;
#pragma omp parallel private()
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* Drr, Dpp and divV computation */
      Drr[i] = (vr[i+1]-vr[i])*InvDiffRsup1D[i];
      Dpp[i] = 0.5*(vr[i+1]+vr[i])*InvRmed1D[i];
      divergence[i]  = (vr[i+1]*Rsup1D[i]-vr[i]*Rinf1D[i])*InvDiffRsup1D[i]*InvRmed1D[i];
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* Drp computation */
      Drp[i] = 0.5*Rinf1D[i]*(vt[i]*InvRmed1D[i]-vt[i-1]*InvRmed1D[i-1])*InvDiffRmed1D[i];
    }
  }
  
#pragma omp parallel private(viscosity)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* TAUrr and TAUpp computation */
      viscosity = VISCOSITY;
      if (ViscosityAlpha)
        viscosity = ALPHAVISCOSITY*SOUNDSPEED1D[i]*SOUNDSPEED1D[i]*pow(Rmed1D[i], 1.5);
      Trr[i] = 2.0*rho[i]*viscosity*(Drr[i]-onethird*divergence[i]);
      Tpp[i] = 2.0*rho[i]*viscosity*(Dpp[i]-onethird*divergence[i]);
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* TAUrp computation */
      viscosity = VISCOSITY;
      if (ViscosityAlpha)
	viscosity = ALPHAVISCOSITY*SOUNDSPEED1D[i]*SOUNDSPEED1D[i]*pow(Rmed1D[i], 1.5);
      Trp[i] = 2.0*0.5*(rho[i]+rho[i-1])*viscosity*Drp[i];
    }
  }
  Trp[1]=Trp[nr-1]=0.; /* For angular momentum conservation. */
  if (IAmTheFirst)
    Trp[IINNER+CPUOVERLAP] = Trp_in;
  if (IAmTheLast)
    Trp[IINNER+GLOBALNRAD-CPUOVERLAP] = Trp_out;

				/* Now we can update velocities */
				/* with the viscous source term */
				/* of Navier-Stokes equations.  */
#pragma omp parallel private()
  {
#pragma omp for nowait
    for (i = 1; i < nr-1; i++) {/* vtheta first */
      vt[i] += DeltaT*InvRmed1D[i]*( (Rsup1D[i]*Trp[i+1]-Rinf1D[i]*Trp[i])*InvDiffRsup1D[i]+\
				     0.5*(Trp[i]+Trp[i+1]) )/rho[i];
    }
#pragma omp for nowait
    for (i = 1; i < nr; i++) {	/* and now vrad */
      vr[i] += DeltaT*InvRinf1D[i]*( (Rmed1D[i]*Trr[i]-Rmed1D[i-1]*Trr[i-1])*InvDiffRmed1D[i]-\
				   0.5*(Tpp[i]+Tpp[i-1]) )/(0.5*(rho[i]+rho[i-1]));
    }
				/* Now we impose the velocity at the boundaries  */
#pragma omp single
    {
      VKepIn  = sqrt(G*1.0/Rmed1D[0])*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(Rmed1D[0],2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      VKepOut = sqrt(G*1.0/Rmed1D[nr-1])*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(Rmed1D[nr-1],2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      i = 0;
      vt[i] = VKepIn-Rmed1D[0]*OmegaFrame;
      // if (CPU_Rank == 0) vt[i] = VKepIn-Rmed1D[0]*OmegaFrame;
      i = nr-1;
      vt[i] = VKepOut-Rmed1D[nr-1]*OmegaFrame;
      // if (CPU_Rank == CPU_Number-1) vt[i] = VKepOut-Rmed1D[nr-1]*OmegaFrame;
    }
  }
}

