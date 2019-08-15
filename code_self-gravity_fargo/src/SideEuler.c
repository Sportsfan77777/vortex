#include "mp.h"

extern boolean OpenInner, NonReflecting, OuterSourceMass, Evanescent;
extern boolean SelfGravity, SGZeroMode, Adiabatic;
extern Pair DiskOnPrimaryAcceleration;
extern int dimfxy;
real Hp0, Hg0, Ht0;

real GasTotalMass (array)
     PolarGrid *array;
{
  int i, j, ns;
  real *density, total = 0.0, fulltotal=0.0;
  ns = array->Nsec;
  density = array->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      total += Surf[i]*density[j+i*ns];
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}

real GasMomentum (Density, Vtheta)
     PolarGrid *Density, *Vtheta;
{
  int i,j,ns;
  real *density, *vtheta, total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vtheta = Vtheta->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &stat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 1; j < ns; j++) {
      total += Surf[i]*(density[j+i*ns]+density[j-1+i*ns])*Rmed[i]*(vtheta[j+i*ns]+OmegaFrame*Rmed[i]);
    }
    total += Surf[i]*(density[i*ns]+density[i*ns+ns-1])*Rmed[i]*(vtheta[i*ns]+OmegaFrame*Rmed[i]);
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return 0.5*fulltotal;
}

real GasTotalEnergy (Density, Vrad, Vtheta, Energy)
     PolarGrid *Density, *Vrad, *Vtheta, *Energy;
{
  int i, j, l, ns;
  real *density, *vrad, *vtheta, *energy, *pot;
  real vr_cent, vt_cent;
  real total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  energy = Energy->Field;
  pot = Potential->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &stat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* centered-in-cell radial velocity */
      vr_cent = (Rmed[i]-Rinf[i])*vrad[l+ns] + (Rsup[i]-Rmed[i])*vrad[l];
      vr_cent /= (Rsup[i]-Rinf[i]);
      /* centered-in-cell azimuthal velocity */
      if (j < ns-1)
	vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
	vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      total += 0.5*Surf[i]*density[l]*(vr_cent*vr_cent + vt_cent*vt_cent) + \
	Surf[i]*energy[l] -						\
	Surf[i]*density[l]*pot[l];
      /* Gas total energy is the sum of its kinematic energy, internal energy */
      /* and gravitational potential energy, including self-gravity */
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}


void CheckMomentumConservation (Density, Vtheta, sys)
     PolarGrid *Density, *Vtheta;
     PlanetarySystem *sys;
{
  FILE *fichmom;
  char name[80];
  int k;
  real totalmomentum, plmom;
  real xplanet, yplanet, vxplanet, vyplanet;
  real rpl, thetapl, vazimpl, masspl;
  real gasmom, planetsmom;
  gasmom = GasMomentum (Density, Vtheta);
  planetsmom = 0.;
  
  for ( k = 0; k < sys->nb; k++ ) {
    xplanet     = sys->x[k];
    yplanet     = sys->y[k];
    rpl         = sqrt( xplanet*xplanet + yplanet*yplanet );
    thetapl     = atan2 (yplanet, xplanet);
    vxplanet    = sys->vx[k];
    vyplanet    = sys->vy[k];
    vazimpl     = -vxplanet*sin(thetapl) + vyplanet*cos(thetapl);
    masspl      = sys->mass[k];
    plmom       = masspl*rpl*vazimpl;
    planetsmom += plmom;
  }
  totalmomentum = gasmom + planetsmom;
  if ( PhysicalTime < 1e-10 ) {
    Hp0 = plmom;
    Hg0 = gasmom;
    Ht0 = totalmomentum;
    printf("time = %lg, Hp0 = %lg, Hg0 = %lg et Ht0 = %lg\n", PhysicalTime, Hp0, Hg0, Ht0);
  }
  if (!CPU_Master) return;
  sprintf (name, "%s%s.dat", OUTPUTDIR, "Momentum");
  fichmom = fopen(name, "a");
  if (fichmom == NULL) {
    fprintf (stderr, "Can't write 'Momentum.dat' file. Aborting.\n");
    prs_exit (1);
  }
  plmom = fabs (plmom - Hp0);
  gasmom = fabs (gasmom - Hg0);
  totalmomentum = fabs (totalmomentum - Ht0);
  fprintf (fichmom, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", PhysicalTime, plmom, gasmom, totalmomentum, totalmomentum / Ht0);
  fclose (fichmom);
}


void DivisePolarGrid (Num, Denom, Res)
     PolarGrid *Num, *Denom, *Res;
{
  int i,j,l,nr,ns;
  real *num, *denom, *res;
  num = Num->Field;
  denom=Denom->Field;
  res = Res->Field;
  ns = Res->Nrad;
  nr = Res->Nsec;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      res[l] = num[l]/(denom[l]+1e-20);
    }
  }
}

void InitComputeAccel ()
{
  int i, j, l, nr, ns;
  real *abs, *ord;
  CellAbscissa = CreatePolarGrid (NRAD,NSEC,"abscissa");
  CellOrdinate = CreatePolarGrid (NRAD,NSEC,"ordinate");
  nr = CellAbscissa->Nrad;
  ns = CellAbscissa->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      abs[l] = Rmed[i] * cos(2.0*PI*(real)j/(real)ns);
      ord[l] = Rmed[i] * sin(2.0*PI*(real)j/(real)ns);
    }
  }
}
  
Pair ComputeAccel (force, Rho, x, y, rsmoothing, mass)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
{
  Pair acceleration;
  ComputeForce (force, Rho, x, y, rsmoothing, mass, dimfxy);
  if (ExcludeHill) {
    acceleration.x = force->fx_ex_inner+force->fx_ex_outer;
    acceleration.y = force->fy_ex_inner+force->fy_ex_outer;
  } else {
    acceleration.x = force->fx_inner+force->fx_outer;
    acceleration.y = force->fy_inner+force->fy_outer;
  }
  return acceleration;
}

void OpenBoundary (Vrad, Rho, Energy)
     PolarGrid *Vrad, *Rho, *Energy;
{
  int i,j,l,ns,nr;
  real *rho, *vr, *energy;
  if (CPU_Rank != 0) return;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  energy = Energy->Field;
  i = 1;
#pragma omp parallel for private(l)
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l-ns] = rho[l] ;		// copy first ring into ghost ring
    energy[l-ns] = energy[l];
    if ((vr[l+ns] > 0.0) || (rho[l] < SigmaMed[0]))
      vr[l] = 0.0; // we just allow outflow [inwards]
    else
      vr[l] = vr[l+ns];
  }
  /*
  i = nr-1;
#pragma omp parallel for private(l)
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l] = rho[l-ns];		// copy first ring into ghost ring
    energy[l] = energy[l-ns];
    if ((vr[l-ns] < 0.0) || (rho[l] < SigmaMed[nr-2]))
      vr[l] = 0.0; // we just allow outflow [outwards]
    else
      vr[l] = vr[l-ns];
  }
  */
}

void NonReflectingBoundary (Vrad, Rho, Energy)
     PolarGrid *Vrad, *Rho, *Energy;
{
  int i, j, l, ns, nr, jp, lp, i_angle;
  real *rho, *vr, *cs, *energy;
  real dangle, mean, vr_med;
  real cs0, cs1, csnrm1, csnrm2;
  cs = SoundSpeed->Field;
  energy = Energy->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  if (CPU_Rank == 0) {
    cs0 = 0.0;
    cs1 = 0.0;
    for (j=0; j<ns; j++) {
      cs0 += cs[j];
      cs1 += cs[ns+j];
    }
    cs0 /= (real)ns;
    cs1 /= (real)ns;
    i = 1;		 /* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[1],-1.5)-1.0)/(.5*(cs0+cs1));
    dangle *= (Rmed[1]-Rmed[0]);
    i_angle = (int)(dangle/2.0/PI*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j+i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp;
      rho[lp] = rho[l];		/* copy first ring into ghost ring */
      energy[lp] = energy[l];   /* copy first ring into ghost ring */
      vr_med = -cs[l]*(rho[l]-SigmaMed[1])/SigmaMed[1];
      vr[l] = 2.*vr_med-vr[l+ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j] += SigmaMed[0]-mean;
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += energy[j];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      energy[j] += EnergyMed[0]-mean;
    }
  }
  if (CPU_Rank == CPU_Highest) {
    csnrm2 = 0.0;
    csnrm1 = 0.0;
    for (j=0; j<ns; j++) {
      csnrm2 += cs[(nr-2)*ns+j];
      csnrm1 += cs[(nr-1)*ns+j];
    }
    csnrm1 /= (real)ns;
    csnrm2 /= (real)ns;
    i = nr-1;		 /* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[nr-2],-1.5)-1.0)/(.5*(csnrm1+csnrm2));
    dangle *= (Rmed[nr-1]-Rmed[nr-2]);
    i_angle = (int)(dangle/2.0/PI*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j-i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp+(i-1)*ns;
      rho[l] = rho[lp];		/* copy first ring into ghost ring */
      energy[l] = energy[lp];	/* copy first ring into ghost ring */
      vr_med = cs[l]*(rho[l-ns]-SigmaMed[nr-2])/SigmaMed[nr-2];
      vr[l] = 2.*vr_med-vr[l-ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j+(nr-1)*ns] += SigmaMed[nr-1]-mean;
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += energy[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      energy[j+(nr-1)*ns] += EnergyMed[nr-1]-mean;
    }
  }
}

void EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real vrad0, vtheta0, viscosity, dens0, energ0;
  real DRMIN, DRMAX, damping, Tin, Tout, lambda;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* Orbital period at inner and outer boundary */
  Tin = 2.0*PI*pow(GlobalRmed[0],3./2);
  Tout = 2.0*PI*pow(GlobalRmed[GLOBALNRAD-1],3./2);
  /* DRMIN AND DRMAX are global Radii boundaries of killing wave zones */
  DRMIN = GlobalRmed[0]*1.25;
  DRMAX = GlobalRmed[GLOBALNRAD-1]*0.84;

  lambda = 0.0;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    if ( (Rmed[i] < DRMIN) || (Rmed[i] > DRMAX) ) {
      /* Damping operates only inside the wave killing zones */
      if (Rmed[i] < DRMIN) {
	damping = (Rmed[i]-DRMIN)/(GlobalRmed[0]-DRMIN);
	lambda = damping*damping*10.0*step/Tin;
      }
      if (Rmed[i] > DRMAX) {
	damping = (Rmed[i]-DRMAX)/(GlobalRmed[GLOBALNRAD-1]-DRMAX);
	lambda = damping*damping*10.0*step/Tout;
      }
      viscosity = FViscosity (Rmed[i]);
      if (!SelfGravity) {
	vtheta0 = sqrt ( G*1.0/Rmed[i] *				\
			 ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			   pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) );
      }
      if (SelfGravity) {
	vtheta0 = sqrt (  G*1.0/Rmed[i] *				\
			  ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			    pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) - \
			  Rmed[i]*GLOBAL_AxiSGAccr[i+IMIN] );
      }
      /* this could be refined if CentrifugalBalance is used... */
      vtheta0 -= Rmed[i]*OmegaFrame;
      vrad0 = -3.0*viscosity/Rmed[i]*(-SIGMASLOPE+.5);
      dens0 = SigmaMed[i];
      energ0 = EnergyMed[i];
      
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vrad[l]   = (vrad[l]+lambda*vrad0)/(1.0+lambda);
	vtheta[l] = (vtheta[l]+lambda*vtheta0)/(1.0+lambda);
	dens[l]   = (dens[l]+lambda*dens0)/(1.0+lambda);
	if (Adiabatic)
	  energ[l]  = (energ[l]+lambda*energ0)/(1.0+lambda);
      }
    }
  }
}

void ApplyOuterSourceMass (Rho, Vrad)
     PolarGrid *Rho, *Vrad;
{
  int i, j, l, nr, ns;
  real *rho, average_rho = 0.0, *vr, penul_vr;
  if (CPU_Rank != CPU_Highest) return;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    average_rho += rho[l];
  }
  average_rho /= (real)ns;
  average_rho = SigmaMed[nr-1]-average_rho;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l] += average_rho;
  }
  i = nr-1;
  penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[nr-1]/1.0),-SIGMASLOPE);
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    vr[l] = penul_vr;
  }
}

void ApplySubKeplerianBoundary (Vtheta)
     PolarGrid *Vtheta;
{
  int i, j, l, nr, ns;
  real VKepIn, VKepOut;
  real *vt;
  vt = Vtheta->Field;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
#pragma omp single
  {
    if ( !SelfGravity ) {
      VKepIn = sqrt (  G*1.0/Rmed[0] *					\
		       ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			 pow(AspectRatio(Rmed[0]),2.0)*pow(Rmed[0],2.0*FLARINGINDEX) ) );
      VKepOut = sqrt (  G*1.0/Rmed[nr-1] *				\
			( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			  pow(AspectRatio(Rmed[nr-1]),2.0)*pow(Rmed[nr-1],2.0*FLARINGINDEX) ) );
    }
    else {
      if ( !SGZeroMode )
	mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	GLOBAL_AxiSGAccr = SG_Accr;
      VKepIn = sqrt (  G*1.0/Rmed[0] *					\
		       ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			 pow(AspectRatio(Rmed[0]),2.0)*pow(Rmed[0],2.0*FLARINGINDEX) ) - \
		       Rmed[0]*GLOBAL_AxiSGAccr[0] );
      VKepOut = sqrt (  G*1.0/Rmed[nr-1] *				\
			( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			  pow(AspectRatio(Rmed[nr-1]),2.0)*pow(Rmed[nr-1],2.0*FLARINGINDEX) ) - \
			Rmed[nr-1]*GLOBAL_AxiSGAccr[nr-1+IMIN] );
    }
    /* ----- */
    /* i = 0 */
    /* ----- */
    if ( CPU_Rank == 0 ) {
      i = 0;
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vt[l] = VKepIn-Rmed[i]*OmegaFrame;
      }
    }
    /* ---------- */
    /* i = nr - 1 */
    /* ---------- */
    if ( CPU_Rank == CPU_Highest ) {
      i = nr - 1;
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vt[l] = VKepOut-Rmed[i]*OmegaFrame;
      }
    }
  }
}

void ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, step)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
{
  if (OpenInner == YES) OpenBoundary (Vrad, Rho, Energy);
  if (NonReflecting == YES) {
    if (Adiabatic)
      ComputeSoundSpeed (Rho, Energy);
    NonReflectingBoundary (Vrad, Rho, Energy);
  }
  if (Evanescent == YES) EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step);
  if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
}

void CorrectVtheta (vtheta, domega)
     PolarGrid *vtheta;
     real domega;
{
  int i, j, l, nr, ns;
  real *vt;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  vt = vtheta->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vt[l] -= domega*Rmed[i];
    }
  }
}
