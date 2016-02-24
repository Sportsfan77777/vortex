#include "mp.h"

static PolarGrid *DivergenceVelocity, *DRR, *DRP, *DPP, *TAURR, *TAUPP, *TAURP;
real Trp_in, Trp_out;

void InitViscosity ()
{
  DivergenceVelocity = CreatePolarGrid(NRAD, NSEC, "DivV");
  DRR                = CreatePolarGrid(NRAD, NSEC, "Drr");
  DRP                = CreatePolarGrid(NRAD, NSEC, "Drp");
  DPP                = CreatePolarGrid(NRAD, NSEC, "Dpp");
  TAURR              = CreatePolarGrid(NRAD, NSEC, "TAUrr");
  TAURP              = CreatePolarGrid(NRAD, NSEC, "TAUrp");
  TAUPP              = CreatePolarGrid(NRAD, NSEC, "TAUpp");
}

void ViscousTerms (RadialVelocity, AzimuthalVelocity, Rho, DeltaT)
PolarGrid *RadialVelocity, *AzimuthalVelocity, *Rho;
real DeltaT;
{
  int i,j,l,nr,ns,lip,ljp,ljm,lim,ljmim;
  real *rho, *vr, *vt;
  real *Drr, *Drp, *Dpp, *divergence;
  real *Trr, *Trp, *Tpp;
  real dphi, VKepIn, VKepOut, onethird, invdphi;
  real viscosity;
  nr  = Rho->Nrad;
  ns  = Rho->Nsec;
  rho = Rho->Field;
  vr  = RadialVelocity->Field;
  vt  = AzimuthalVelocity->Field;
  divergence = DivergenceVelocity->Field;
  Drr = DRR->Field;
  Drp = DRP->Field;
  Dpp = DPP->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  dphi = 2.0*M_PI/(real)ns;
  invdphi = 1.0/dphi;
  onethird = 1.0/3.0;
#pragma omp parallel private(l,lip,ljp,j,ljm,lim)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* Drr, Dpp and divV computation */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lip = l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	Drr[l] = (vr[lip]-vr[l])*InvDiffRsup[i];
	Dpp[l] = (vt[ljp]-vt[l])*invdphi*InvRmed[i]+0.5*(vr[lip]+vr[l])*InvRmed[i];
	divergence[l]  = (vr[lip]*Rsup[i]-vr[l]*Rinf[i])*InvDiffRsup[i]*InvRmed[i];
	divergence[l] += (vt[ljp]-vt[l])*invdphi*InvRmed[i];
      }
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* Drp computation */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	lip = l+ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	lim = l-ns;
	Drp[l] = 0.5*(Rinf[i]*(vt[l]*InvRmed[i]-vt[lim]*InvRmed[i-1])*InvDiffRmed[i]+\
		      (vr[l]-vr[ljm])*invdphi*InvRinf[i]);
      }
    }
  }
#pragma omp parallel private(l,ljmim,j,ljm,lim,viscosity)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* TAUrr and TAUpp computation */
      viscosity = VISCOSITY;
      if (ViscosityAlpha)
	viscosity = ALPHAVISCOSITY*SOUNDSPEED[i]*SOUNDSPEED[i]*pow(Rmed[i], 1.5);
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	Trr[l] = 2.0*rho[l]*viscosity*(Drr[l]-onethird*divergence[l]);
	Tpp[l] = 2.0*rho[l]*viscosity*(Dpp[l]-onethird*divergence[l]);
      }
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* TAUrp computation */
      viscosity = VISCOSITY;
      if (ViscosityAlpha)
	viscosity = ALPHAVISCOSITY*SOUNDSPEED[i]*SOUNDSPEED[i]*pow(Rmed[i], 1.5);
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljmim=ljm-ns;
	Trp[l] = 2.0*0.25*(rho[l]+rho[lim]+rho[ljm]+rho[ljmim])*viscosity*Drp[l];
      }
    }
    Trp_in=0.;
    Trp_out=0.;
    i=CPUOVERLAP;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      Trp_in += Trp[l];
    }
    i=NRAD-CPUOVERLAP;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      Trp_out += Trp[l];
    }
    Trp_in /= (real)ns;
    Trp_out /= (real)ns;
  }

				/* Now we can update velocities */
				/* with the viscous source term */
				/* of Navier-Stokes equations   */
#pragma omp parallel private(l,j,lip,ljp,ljm,lim)
  {
#pragma omp for nowait
    for (i = 1; i < nr-1; i++) {	/* vtheta first */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lip = l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vt[l] += DeltaT*InvRmed[i]*( (Rsup[i]*Trp[lip]-Rinf[i]*Trp[l])*InvDiffRsup[i]+\
				     (Tpp[l]-Tpp[ljm])*invdphi+\
                                     0.5*(Trp[l]+Trp[lip]) )\
                                  /(0.5*(rho[l]+rho[ljm]));
      }
    }
#pragma omp for nowait
    for (i = 1; i < nr; i++) {	/* and now vrad */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lip = l+ns;
	lim = l-ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vr[l] += DeltaT*InvRinf[i]*((Rmed[i]*Trr[l]-Rmed[i-1]*Trr[lim])*InvDiffRmed[i]+\
				    (Trp[ljp]-Trp[l])*invdphi-\
				    0.5*(Tpp[l]+Tpp[lim]))/(0.5*(rho[l]+rho[lim]));
      }
    }
				/* Now we impose the velocity at the boundaries */
#pragma omp single
    {
      VKepIn  = sqrt(G*1.0/Rmed[0])*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(Rmed[0],2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      VKepOut = sqrt(G*1.0/Rmed[nr-1])*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(Rmed[nr-1],2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      i = 0;
      if (CPU_Rank == 0) {
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  vt[l] = VKepIn-Rmed[0]*OmegaFrame;
	}
      }
      i = nr-1;
      if (CPU_Rank == CPU_Number-1) {
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  vt[l] = VKepOut-Rmed[nr-1]*OmegaFrame;
	} 
      }
    }
  }
}

