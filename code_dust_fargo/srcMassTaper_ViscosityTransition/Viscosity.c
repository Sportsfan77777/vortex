/** \file Viscosity.c

Calculation of the viscous %force.  The function FViscosity() returns
the (kinematic) viscosity as a function of the radius (it handles all
case: alpha or uniform viscosity, and inner cavity with a different
viscosity).  The update of the velocity is done in ViscousTerm(),
which properly evaluate the stress tensor in 2D cylindrical
coordinates. This file also contains the function AspectRatio(), which
gives the aspect ratio as a function of the radius, in the case of a
temperature jump in the disk (much in the manner as cavities arising
from a viscosity jump are handled, hence the location of this
function). Note that AspectRatio() does not feature the FLARINGINDEX,
which is taken into account by the calling function.

*/

#include "fargo.h"

static PolarGrid *DivergenceVelocity, *DRR, *DRP, *DPP, *TAURR, *TAUPP, *TAURP, *DDivergenceVelocity, *DDRR, *DDRP, *DDPP, *DTAURR, *DTAUPP, *DTAURP;

real FViscosity (rad)
     real rad;
{
  real viscosity, rmin, rmax, scale;
  int i=0;
  viscosity = VISCOSITY;
  if (ViscosityAlpha) {
    while (GlobalRmed[i] < rad) i++;
    viscosity = ALPHAVISCOSITY*GLOBAL_SOUNDSPEED[i]*\
      GLOBAL_SOUNDSPEED[i]*pow(rad, 1.5);
  }
  rmin = CAVITYRADIUS-CAVITYWIDTH*ASPECTRATIO;
  rmax = CAVITYRADIUS+CAVITYWIDTH*ASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) viscosity *= CAVITYRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    viscosity *= exp((rmax-rad)/(rmax-rmin)*log(CAVITYRATIO));
  }
  return viscosity;
}
 
real DFViscosity (rad)
     real rad;
{
  real dviscosity, rmin, rmax, scale;
  int i=0;
  dviscosity = DVISCOSITY;
  if (DViscosityAlpha) {
    while (GlobalRmed[i] < rad) i++;
    dviscosity = DALPHAVISCOSITY*DGLOBAL_SOUNDSPEED[i]*\
      DGLOBAL_SOUNDSPEED[i]*pow(rad, 1.5);
  }
  rmin = CAVITYRADIUS-CAVITYWIDTH*DASPECTRATIO;
  rmax = CAVITYRADIUS+CAVITYWIDTH*DASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) dviscosity *= CAVITYRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    dviscosity *= exp((rmax-rad)/(rmax-rmin)*log(CAVITYRATIO));
  }
  return dviscosity;
}

 
real AspectRatio (rad)
     real rad;
{
  real aspectratio, rmin, rmax, scale;
  aspectratio = ASPECTRATIO;
  rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*ASPECTRATIO;
  rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*ASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) aspectratio *= TRANSITIONRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    aspectratio *= exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
  }
  return aspectratio;
}

real DAspectRatio (rad) /*dust aspect ratio all other parameters are the same as gas only DASPECRATIO is different*/
     real rad;
{
  real aspectratio, rmin, rmax, scale;
  aspectratio = DASPECTRATIO;
  rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*DASPECTRATIO;
  rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*DASPECTRATIO;
  scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
  rmin *= scale;
  rmax *= scale;
  if (rad < rmin) aspectratio *= TRANSITIONRATIO;
  if ((rad >= rmin) && (rad <= rmax)) {
    aspectratio *= exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
  }
  return aspectratio;
}

  
void InitViscosity ()
{
  DivergenceVelocity = CreatePolarGrid(NRAD, NSEC, "DivV");
  DRR                = CreatePolarGrid(NRAD, NSEC, "Drr");
  DRP                = CreatePolarGrid(NRAD, NSEC, "Drp");
  DPP                = CreatePolarGrid(NRAD, NSEC, "Dpp");
  TAURR              = CreatePolarGrid(NRAD, NSEC, "TAUrr");
  TAURP              = CreatePolarGrid(NRAD, NSEC, "TAUrp");
  TAUPP              = CreatePolarGrid(NRAD, NSEC, "TAUpp");
  DDivergenceVelocity = CreatePolarGrid(NRAD, NSEC, "DDivV");
  DDRR                = CreatePolarGrid(NRAD, NSEC, "DDrr");
  DDRP                = CreatePolarGrid(NRAD, NSEC, "DDrp");
  DDPP                = CreatePolarGrid(NRAD, NSEC, "DDpp");
  DTAURR              = CreatePolarGrid(NRAD, NSEC, "DTAUrr");
  DTAURP              = CreatePolarGrid(NRAD, NSEC, "DTAUrp");
  DTAUPP              = CreatePolarGrid(NRAD, NSEC, "DTAUpp");

}

void ViscousTerms (RadialVelocity, AzimuthalVelocity, Rho, DRadialVelocity, DAzimuthalVelocity, DRho, DeltaT)
PolarGrid *RadialVelocity, *AzimuthalVelocity, *Rho, *DRadialVelocity, *DAzimuthalVelocity, *DRho;
real DeltaT;
{
  int i,j,l,nr,ns,lip,ljp,ljm,lim,ljmim;
  real *rho, *vr, *vt, *drho, *dvr, *dvt;
  real *Drr, *Drp, *Dpp, *divergence, *dDrr, *dDrp, *dDpp, *ddivergence;
  real *Trr, *Trp, *Tpp, *dTrr, *dTrp, *dTpp;
  real dphi, VKepIn, VKepOut, DVKepIn, DVKepOut, onethird, invdphi;
  real viscosity,dviscosity;
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

/*dust initialize*/
  drho = DRho->Field;
  dvr  = DRadialVelocity->Field;
  dvt = DAzimuthalVelocity->Field;
  ddivergence = DDivergenceVelocity->Field;
  dDrr = DDRR->Field;
  dDrp = DDRP->Field;
  dDpp = DDPP->Field;
  dTrr = DTAURR->Field;
  dTrp = DTAURP->Field;
  dTpp = DTAUPP->Field;

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
	dDrr[l] = (dvr[lip]-dvr[l])*InvDiffRsup[i];
        dDpp[l] = (dvt[ljp]-dvt[l])*invdphi*InvRmed[i]+0.5*(dvr[lip]+dvr[l])*InvRmed[i];
        ddivergence[l]  = (dvr[lip]*Rsup[i]-dvr[l]*Rinf[i])*InvDiffRsup[i]*InvRmed[i];
        ddivergence[l] += (dvt[ljp]-dvt[l])*invdphi*InvRmed[i];
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
	dDrp[l] = 0.5*(Rinf[i]*(dvt[l]*InvRmed[i]-dvt[lim]*InvRmed[i-1])*InvDiffRmed[i]+\
                      (dvr[l]-dvr[ljm])*invdphi*InvRinf[i]);
      }
    }
  }
#pragma omp parallel private(l,ljmim,j,ljm,lim,viscosity)
  {
#pragma omp for nowait
    for (i = 0; i < nr; i++) {	/* TAUrr and TAUpp computation */
      viscosity = FViscosity (Rmed[i]);
      dviscosity = DFViscosity (Rmed[i]);
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	Trr[l] = 2.0*rho[l]*viscosity*(Drr[l]-onethird*divergence[l]);
	Tpp[l] = 2.0*rho[l]*viscosity*(Dpp[l]-onethird*divergence[l]);
	dTrr[l] = 2.0*drho[l]*dviscosity*(dDrr[l]-onethird*ddivergence[l]);
        dTpp[l] = 2.0*drho[l]*dviscosity*(dDpp[l]-onethird*ddivergence[l]);
      }
    }
#pragma omp for
    for (i = 1; i < nr; i++) {	/* TAUrp computation */
      viscosity = FViscosity (Rmed[i]);
      dviscosity = DFViscosity (Rmed[i]);
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljmim=ljm-ns;
	Trp[l] = 2.0*0.25*(rho[l]+rho[lim]+rho[ljm]+rho[ljmim])*viscosity*Drp[l];
	dTrp[l] = 2.0*0.25*(drho[l]+drho[lim]+drho[ljm]+drho[ljmim])*dviscosity*dDrp[l];
      }
    }
  }
				/* Now we can update velocities */
				/* with the viscous source term */
				/* of Navier-Stokes equations */
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
	vt[l] += DeltaT*InvRmed[i]*((Rsup[i]*Trp[lip]-Rinf[i]*Trp[l])*InvDiffRsup[i]+\
				    (Tpp[l]-Tpp[ljm])*invdphi+\
				    0.5*(Trp[l]+Trp[lip]))/(0.5*(rho[l]+rho[ljm]));
	dvt[l] += DeltaT*InvRmed[i]*((Rsup[i]*dTrp[lip]-Rinf[i]*dTrp[l])*InvDiffRsup[i]+\
                                    (dTpp[l]-dTpp[ljm])*invdphi+\
                                    0.5*(dTrp[l]+dTrp[lip]))/(0.5*(drho[l]+drho[ljm]));
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
	dvr[l] += DeltaT*InvRinf[i]*((Rmed[i]*dTrr[l]-Rmed[i-1]*dTrr[lim])*InvDiffRmed[i]+\
                                    (dTrp[ljp]-dTrp[l])*invdphi-\
                                    0.5*(dTpp[l]+dTpp[lim]))/(0.5*(drho[l]+drho[lim]));
      }
    }
				/* Now we impose the velocity at the boundaries */
#pragma omp single
    {
      VKepIn  = sqrt(G*1.0/Rmed[0])*\
	sqrt(1.0-pow(AspectRatio(Rmed[0]),2.0)*	\
	     pow(Rmed[0],2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      VKepOut = sqrt(G*1.0/Rmed[nr-1])*\
	sqrt(1.0-pow(AspectRatio(Rmed[nr-1]),2.0)*	\
	     pow(Rmed[nr-1],2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      DVKepIn  = sqrt(G*1.0/Rmed[0])*\
        sqrt(1.0-pow(DAspectRatio(Rmed[0]),2.0)* \
             pow(Rmed[0],2.0*DFLARINGINDEX)*\
             (1.+SIGMASLOPE-2.0*DFLARINGINDEX));
      DVKepOut = sqrt(G*1.0/Rmed[nr-1])*\
        sqrt(1.0-pow(DAspectRatio(Rmed[nr-1]),2.0)*      \
             pow(Rmed[nr-1],2.0*DFLARINGINDEX)*\
             (1.+SIGMASLOPE-2.0*DFLARINGINDEX));

      i = 0;
      if (CPU_Rank == 0) {
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  vt[l] = VKepIn-Rmed[0]*OmegaFrame;
//	  vt[l] = vt[l+ns]+(Rmed[1]-Rmed[0])*OmegaFrame;
	  dvt[l] = DVKepIn-Rmed[0]*OmegaFrame;
	}
      }
      i = nr-1;
      if (CPU_Rank == CPU_Number-1) {
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  vt[l] = VKepOut-Rmed[nr-1]*OmegaFrame;
//	  vt[l]=vt[l-ns]+(Rmed[nr-2]-Rmed[nr-1])*OmegaFrame;
	  dvt[l] = DVKepOut-Rmed[nr-1]*OmegaFrame;
	} 
      }
    }
  }
}

