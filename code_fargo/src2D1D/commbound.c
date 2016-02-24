#include "mp.h"

static real *SendInnerBoundary;
static real *SendOuterBoundary;
static real *RecvInnerBoundary;
static real *RecvOuterBoundary;

static int allocated_com = 0;
static int size_com;

extern real OmegaFrame;

static real *dq1D_;
real rhos_test=0.;
real F_test=0.;


void InitCompute1DGhostsVr () 
{
  dq1D_         = (real *)malloc(NRAD1D*sizeof(real));
}


void Compute1DGhostsVr (Density1D,Vrad1D,Density2D,Vrad2D,dt)
     PolarGrid1D *Density1D, *Vrad1D;
     PolarGrid *Density2D, *Vrad2D;
     real dt;
{
  int nr,ns,i,i1D,i2D,j,l,nr1D;
  real* dens1D, *dens2D, *vr1D, *vr2D, *dens2DStar;
  real F,a,b,c,Delta,u,dqm,dqp,dq1D_m;
  ComputeStarRad (Density2D, Vrad2D, RhoStar, dt);
  dens1D  = Density1D->Field;
  vr1D    = Vrad1D->Field;
  dens2D  = Density2D->Field;
  dens2DStar = RhoStar->Field;
  vr2D    = Vrad2D->Field;
  nr = Density2D->Nrad;
  ns = Density2D->Nsec;
  nr1D = Density1D->Nrad;
  if (IAmTheFirst || IAmTheLast) {
    for (i1D = 0; i1D < nr1D; i1D++) {    /* to be improved. */
      if ((i1D == 0) || (i1D == nr-1)) dq1D_[i1D] = 0.0;
      else {
	dqm = (dens1D[i1D]-dens1D[i1D-1])*InvDiffRmed1D[i1D];
	dqp = (dens1D[i1D+1]-dens1D[i1D])*InvDiffRmed1D[i1D+1];
	if (dqp * dqm > 0.0)
	  dq1D_[i1D] = 2.0*dqp*dqm/(dqp+dqm);
	else
	  dq1D_[i1D] = 0.0;
      }
    }
  }
  if ((RMAX1D>=RMAX) && IAmTheLast) {
    for(i = 0; i <= CPUOVERLAP; i++) {
      i2D=i+nr-2*CPUOVERLAP;
      i1D=i2D+IINNER+IMIN;
      if(Rmed[i2D]!=Rmed1D[i1D]) { /* little test */
	fprintf(stderr,"Error in Compute1DGhostsVr:\n");
	fprintf(stderr,"i2D=%d, i1D=%d, Rmed[i2D]=%.6f, Rmed1D[i1D]=%.6f\n",i2D,i1D,Rmed[i2D],Rmed1D[i1D]);
	exit (1);
      }
      F=0.;
      for(j = 0; j < ns; j++) {
	l=j+i2D*ns;
	F+=dens2DStar[l]*vr2D[l];
      }
      F/=ns;
      if(F>0.) {
	dq1D_m = (i1D>0) ? dq1D_[i1D-1] : 0. ;
	a = -dt*0.5*dq1D_m;
	b = dens1D[i1D-1]+(Rmed1D[i1D]-Rmed1D[i1D-1])*0.5*dq1D_m;
	c = -F;
      } else {
	a = -dt*0.5*dq1D_[i1D];
	b = dens1D[i1D]-(Rmed1D[i1D+1]-Rmed1D[i1D])*0.5*dq1D_[i1D];
	c = -F;
      }
      Delta = b*b-4.*a*c;
      u = -4.*a*c/b/b;
      if (Delta<0.) {
	message (" Delta<0 in Compute1DGhostsVr. Exited.\n");
	prs_exit (1);
      }
      if ( (fabs(u)<4.e-16) || ( (fabs(a/b)<1.e-14) && (fabs(u)<1.e-14) ) ) {
	vr1D[i1D] = -c/b;                   /* first order equation */
      }
      else if ( fabs(u) < 1.e-5 ) {
	vr1D[i1D] = b/2./a*( 1./2.*u - 1./8.*u*u + 1./16.*u*u*u);
      }
      else if ( fabs(u) < 1.2e-3 ) {
	vr1D[i1D] = b/2./a*( 1./2.*u - 1./8.*u*u + 1./16.*u*u*u - 5./128.*u*u*u*u + 7./256.*u*u*u*u*u );
      }
      else {
	vr1D[i1D] = (-b+sqrt(Delta))/2./a;  /* the only physical solution. */
      }
      /* The above physical solution is expanded if u<<1
	 to preserve numerical precision of order 10^16 */
    }
  }
  if ((RMIN1D<=RMIN) && IAmTheFirst) {
    for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
      i2D=i;
      i1D=i2D+IINNER+IMIN;
      F=0.;
      for(j = 0; j < ns; j++) {
	l=j+i2D*ns;
	F+=dens2DStar[l]*vr2D[l];
      }
      F/=ns;
      if(F>0.) {
	a = -dt*0.5*dq1D_[i1D-1];
	b = dens1D[i1D-1]+(Rmed1D[i1D]-Rmed1D[i1D-1])*0.5*dq1D_[i1D-1];
	c = -F;
      } else {
	a = -dt*0.5*dq1D_[i1D];
	b = dens1D[i1D]-(Rmed1D[i1D+1]-Rmed1D[i1D])*0.5*dq1D_[i1D];
	c = -F;
      }
      Delta = b*b-4*a*c;
      u = -4.*a*c/b/b;
      if (Delta<0.) {
	message ("Delta<0 in Compute1DGhostsVr. Exited.\n");
	prs_exit (1);
      }
      if ( (fabs(u)<4.e-16) || ( (fabs(a/b)<1.e-14) && (fabs(u)<1.e-14) ) )
	vr1D[i1D] = -c/b;                   /* first order equation */
      else if ( fabs(u) < 1.e-5 )
	vr1D[i1D] = b/2./a*( 1./2.*u - 1./8.*u*u + 1./16.*u*u*u);
      else if ( fabs(u) < 1.2e-3 )
	vr1D[i1D] = b/2./a*( 1./2.*u - 1./8.*u*u + 1./16.*u*u*u - 5./128.*u*u*u*u + 7./256.*u*u*u*u*u );
      else
	vr1D[i1D] = (-b+sqrt(Delta))/2./a;  /* the only physical solution. */
      /* The above physical solution is expanded if u<<1
	 to preserve numerical precision of order 10^16 */
    }
  }
}


void Communicate1D2D (Density1D,Vrad1D,Vtheta1D,Label1D,Density2D,Vrad2D,Vtheta2D,Label2D,dt)
     PolarGrid1D *Density1D, *Vrad1D, *Vtheta1D, *Label1D;
     PolarGrid *Density2D, *Vrad2D, *Vtheta2D, *Label2D;
     real dt;
{
  //  int ii;
  int i,i1D,i2D,j,l,ljm,ljp,lim,lip,nr,ns;
  real *dens1D, *dens2D, *vr1D, *vr2D, *vt1D, *vt2D, *lab1D, *lab2D;
  real Dens2D_Moy[2*CPUOVERLAP], Labdens2D_Moy[2*CPUOVERLAP]; 
  real ThetaMom2D_Moy[2*CPUOVERLAP], RadMomP2D_Moy[2*CPUOVERLAP], RadMomM2D_Moy[2*CPUOVERLAP];
  real Dens_true[2*CPUOVERLAP], Labdens_true[2*CPUOVERLAP];
  real ThetaMom_true[2*CPUOVERLAP], RadMomP_true[2*CPUOVERLAP], RadMomM_true[2*CPUOVERLAP];
  real dDens;
  nr=Density2D->Nrad;
  ns=Density2D->Nsec;
  dens1D = Density1D->Field;
  vr1D   = Vrad1D->Field;
  vt1D   = Vtheta1D->Field;
  lab1D  = Label1D->Field;
  dens2D = Density2D->Field;
  vr2D   = Vrad2D->Field;
  vt2D   = Vtheta2D->Field;
  lab2D  = Label2D->Field;

  if ((RMAX1D >= RMAX) && IAmTheLast) {  /* OUTER BOUNDARY */
    /* Calculation of the azimuthally-averaged quantities in the 2D-Grid
     and calculation of the physically pertinent quantities. */
    for(i = 0; i < 2*CPUOVERLAP; i++) {
      Dens2D_Moy[i]=ThetaMom2D_Moy[i]=RadMomM2D_Moy[i]=RadMomP2D_Moy[i]=Labdens2D_Moy[i]=0.;
      i2D=i+nr-2*CPUOVERLAP;
      i1D=i2D+IINNER+IMIN;
      for(j = 0; j < ns; j++) {
	l=j+i2D*ns;
	lip=l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i2D*ns;
	Dens2D_Moy[i] += dens2D[l];
	ThetaMom2D_Moy[i] += dens2D[l]*((vt2D[l]+vt2D[ljp])/2.+Rmed[i2D]*OmegaFrame)*Rmed[i2D];
	//ThetaMom2D_Moy[i] += (vt2D[ljp]+Rmed[i2D]*OmegaFrame)*(dens2D[l]+dens2D[ljp])*0.5*Rmed[i2D];
	RadMomM2D_Moy[i] += dens2D[l]*vr2D[l];
	RadMomP2D_Moy[i] += dens2D[l]*vr2D[lip];
	if(AdvecteLabel==YES) Labdens2D_Moy[i] += dens2D[l]*lab2D[l];
      }
      Dens2D_Moy[i] /= ns;
      ThetaMom2D_Moy[i] /= ns;
      RadMomM2D_Moy[i] /= ns;
      RadMomP2D_Moy[i] /= ns;
      if(AdvecteLabel==YES) Labdens2D_Moy[i] /= ns;
      if(i < CPUOVERLAP) {
	Dens_true[i] = Dens2D_Moy[i];
	ThetaMom_true[i] = ThetaMom2D_Moy[i];
	RadMomM_true[i] = RadMomM2D_Moy[i];
	RadMomP_true[i] = RadMomP2D_Moy[i];
	if(AdvecteLabel==YES) Labdens_true[i] = Labdens2D_Moy[i];
      } else {
	Dens_true[i] = dens1D[i1D];
	ThetaMom_true[i] = dens1D[i1D]*(vt1D[i1D]+Rmed1D[i1D]*OmegaFrame)*Rmed1D[i1D];
	RadMomM_true[i] = dens1D[i1D]*vr1D[i1D];
	RadMomP_true[i] = dens1D[i1D]*vr1D[i1D+1];
	if(AdvecteLabel==YES) Labdens_true[i] = dens1D[i1D]*lab1D[i1D];
      }
    }
    /* 1D-Grid ghosts : */
    for(i = 0; i < CPUOVERLAP; i++) {
      i2D=i+nr-2*CPUOVERLAP;
      i1D=i2D+IINNER+IMIN;
      if(Rmed[i2D]!=Rmed1D[i1D]) { /* little test */
	fprintf(stderr,"Error in Communicate1D2D:\n");
	fprintf(stderr,"i2D=%d, i1D=%d, Rmed[i2D]=%.6f, Rmed1D[i1D]=%.6f\n",i2D,i1D,Rmed[i2D],Rmed1D[i1D]);
	exit (1);
      }
      dens1D[i1D] = Dens_true[i];
      vt1D[i1D] = ThetaMom_true[i]/dens1D[i1D]/Rmed1D[i1D]-Rmed1D[i1D]*OmegaFrame;
    }
    if(AdvecteLabel==YES) {
      for(i = 0; i < CPUOVERLAP; i++) {
	i2D=i+nr-2*CPUOVERLAP;
	i1D=i2D+IINNER+IMIN;
	lab1D[i1D] = Labdens_true[i]/dens1D[i1D];
      }
    }
    /* 2D-Grid ghost rings : */
    for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
      i2D=i+nr-2*CPUOVERLAP;
      for(j=0; j < ns; j++) {
	l=j+i2D*ns;
	lim=l-ns;
	ljm = l-1;
	if (j == 0) ljm = i2D*ns+ns-1;
	if (i>CPUOVERLAP) 
	  vr2D[l] = vr2D[l]*(dens2D[l]+dens2D[lim])/2.\
	    - (RadMomM2D_Moy[i]+RadMomP2D_Moy[i-1])/2.\
	    + (RadMomM_true[i]+RadMomP_true[i-1])/2.;
   	vt2D[l] = (vt2D[l]+Rmed[i2D]*OmegaFrame)*Rmed[i2D]*0.5*(dens2D[l]+dens2D[ljm])\
	  - ThetaMom2D_Moy[i] + ThetaMom_true[i];
      }
    }
    for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
      i2D=i+nr-2*CPUOVERLAP;
      dDens = - Dens2D_Moy[i] + Dens_true[i];
      for(j=0; j < ns; j++) {
	l=j+i2D*ns;
	if (dens2D[l]>-dDens)
	  dens2D[l] += dDens;
	else {              
	  dens2D[l] *= 0.001;    // to avoid a crash.
	}
	// Old versions :
	//dens2D[l] = dens2D[l] - Dens2D_Moy[i] + Dens_true[i];
	//dens2D[l] = dens2D[l] / Dens2D_Moy[i] * Dens_true[i];
      }
    }
    for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
      i2D=i+nr-2*CPUOVERLAP;
      for(j=0; j < ns; j++) {
	l=j+i2D*ns;
	lim=l-ns;
	ljm = l-1;
	if (j == 0) ljm = i2D*ns+ns-1;
	vt2D[l] = vt2D[l]/(0.5*(dens2D[l]+dens2D[ljm]))/Rmed[i2D]- Rmed[i2D]*OmegaFrame;
	if (i>CPUOVERLAP) 
	  vr2D[l] /= 0.5*(dens2D[l]+dens2D[lim]);
      }
    }
    if(AdvecteLabel==YES) {
      for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
	i2D=i+nr-2*CPUOVERLAP;
	for(j=0; j < ns; j++) {
	  l=j+i2D*ns;
	  lab2D[l] = lab2D[l]*dens2D[l] - Labdens2D_Moy[i] + Labdens_true[i];
	  lab2D[l] /= dens2D[l];
	}
      }
    }
    /* Boundary conditions for the 2D-Grid */
    i1D=nr-1+IINNER+IMIN;
    if(Rmed[nr-1]!=Rmed1D[i1D]) fprintf(stderr," ErroR in Communicate1D2D\n");
    SigmaMed[nr-1] = dens1D[i1D];
    SigmaMed[nr-2] = dens1D[i1D-1];
  }

  if ((RMIN1D <= RMIN ) && IAmTheFirst) { /* INNER BOUNDARY */

    /* Calculation of the azimuthally-averaged quantities in the 2D-Grid
     and calculation of the physically pertinent quantities. */
    for(i = 0; i < 2*CPUOVERLAP; i++) {
      Dens2D_Moy[i]=ThetaMom2D_Moy[i]=RadMomM2D_Moy[i]=RadMomP2D_Moy[i]=Labdens2D_Moy[i]=0.;
      i1D=i+IINNER+IMIN;
      if(Rmed[i]!=Rmed1D[i1D]) {
	fprintf(stderr,"Error in Communicate1D2D:\n");
	fprintf(stderr,"i=%d, i1D=%d, Rmed[i]=%.6f, Rmed1D[i1D]=%.6f\n",i,i1D,Rmed[i],Rmed1D[i1D]);
      }
      for(j = 0; j < ns; j++) {
	l=j+i*ns;
	lip=l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	Dens2D_Moy[i] += dens2D[l];
	ThetaMom2D_Moy[i] += dens2D[l]*((vt2D[l]+vt2D[ljp])/2.+Rmed[i]*OmegaFrame)*Rmed[i];
	RadMomM2D_Moy[i] += dens2D[l]*(Rmed[i]-Rinf[i])*vr2D[l];
	RadMomP2D_Moy[i] += dens2D[l]*(Rsup[i]-Rmed[i])*vr2D[lip];
	if(AdvecteLabel==YES) Labdens2D_Moy[i] += dens2D[l]*lab2D[l];
      }
      Dens2D_Moy[i] /= ns;
      ThetaMom2D_Moy[i] /= ns;
      RadMomM2D_Moy[i] /= ns;
      RadMomP2D_Moy[i] /= ns;
      if(AdvecteLabel==YES) Labdens2D_Moy[i] /= ns;
      if(i >= CPUOVERLAP) {
	Dens_true[i] = Dens2D_Moy[i];
	ThetaMom_true[i] = ThetaMom2D_Moy[i];
	RadMomM_true[i] = RadMomM2D_Moy[i];
	RadMomP_true[i] = RadMomP2D_Moy[i];
	if(AdvecteLabel==YES) Labdens_true[i] = Labdens2D_Moy[i];
      } else {
	Dens_true[i] = dens1D[i1D];
	ThetaMom_true[i] = dens1D[i1D]*(vt1D[i1D]+Rmed[i]*OmegaFrame)*Rmed1D[i1D];
	RadMomM_true[i] = dens1D[i1D]*(Rmed1D[i1D]-Rinf1D[i1D])*vr1D[i1D];
	RadMomP_true[i] = dens1D[i1D]*(Rsup1D[i1D]-Rmed1D[i1D])*vr1D[i1D+1];
	if(AdvecteLabel==YES) Labdens_true[i] = dens1D[i1D]*lab1D[i1D];
      }
    }

    /* 2D-Grid ghost rings : */
    for(i = 0; i < CPUOVERLAP; i++) {
      for(j=0; j < ns; j++) {
	l=j+i*ns;
	lim=l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	if (i>0)
	  vr2D[l] = vr2D[l]*(dens2D[l]*(Rmed[i]-Rinf[i])+dens2D[lim]*(Rinf[i]-Rmed[i-1]))/2.\
		    - (RadMomM2D_Moy[i]+RadMomP2D_Moy[i-1])/2.\
		    + (RadMomM_true[i]+RadMomP_true[i-1])/2.;
	vt2D[l] = (vt2D[l]+Rmed[i]*OmegaFrame)*Rmed[i]*0.5*(dens2D[l]+dens2D[ljm])\
		  - ThetaMom2D_Moy[i] + ThetaMom_true[i];
      }
    }
    for(i = 0; i < CPUOVERLAP; i++) {
      dDens = - Dens2D_Moy[i] + Dens_true[i];
      for(j=0; j < ns; j++) {
	l=j+i*ns;
	if (dens2D[l]>-dDens)
	  dens2D[l] += dDens;
	else {               
	  dens2D[l] *= 0.001;
	}
	// Old versions :
	//dens2D[l] = dens2D[l] - Dens2D_Moy[i] + Dens_true[i];
	//dens2D[l] = dens2D[l] / Dens2D_Moy[i] * Dens_true[i];
      }
    }
    for(i = 0; i < CPUOVERLAP; i++) {
      for(j=0; j < ns; j++) {
	l=j+i*ns;
	lim=l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vt2D[l] = vt2D[l]/(0.5*(dens2D[l]+dens2D[ljm]))/Rmed[i]- Rmed[i]*OmegaFrame;
	if (i>0)
	  vr2D[l] /= 0.5*(dens2D[l]*(Rmed[i]-Rinf[i])+dens2D[lim]*(Rinf[i]-Rmed[i-1]));
      }
    }
    if(AdvecteLabel==YES) {
      for(i = 0; i < CPUOVERLAP; i++) {
	for(j=0; j < ns; j++) {
	  l=j+i*ns;
	  lab2D[l] = lab2D[l]*dens2D[l] - Labdens2D_Moy[i] + Labdens_true[i];
	  lab2D[l] /= dens2D[l];
	}
      }
    }
    /* 1D-Grid ghosts : */
    for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
      i1D=i+IINNER+IMIN;
      dens1D[i1D] = Dens_true[i];
      vt1D[i1D] = ThetaMom_true[i]/dens1D[i1D]/Rmed1D[i1D]-Rmed1D[i1D]*OmegaFrame;
    }
    if(AdvecteLabel==YES) {
      for(i = CPUOVERLAP; i < 2*CPUOVERLAP; i++) {
	i1D=i+IINNER+IMIN;
	lab1D[i1D] = Labdens_true[i]/dens1D[i1D];
      }
    }

    /* Boundary conditions for the 2D-Grid */
    SigmaMed[0] = dens1D[IINNER];
    SigmaMed[1] = dens1D[IINNER+1];
  }
  Compute1DGhostsVr (Density1D,Vrad1D,Density2D,Vrad2D,dt);
}




void AllocateComm () {
  size_com = 3;
  if (AdvecteLabel == YES) size_com = 4;
  size_com *= NSEC * CPUOVERLAP;
  SendInnerBoundary = malloc (size_com * sizeof(real));
  SendOuterBoundary = malloc (size_com * sizeof(real));
  RecvInnerBoundary = malloc (size_com * sizeof(real));
  RecvOuterBoundary = malloc (size_com * sizeof(real));
  if ((SendInnerBoundary == NULL) ||\
      (SendOuterBoundary == NULL) ||\
      (RecvInnerBoundary == NULL) ||\
      (RecvOuterBoundary == NULL)) {
    fprintf (stderr, "CPU %d had not enough memory to allocate communicators.\n", CPU_Rank);
    prs_exit(0);
  }
  allocated_com = 1;
}

void CommunicateBoundaries (Density, Vrad, Vtheta, Label)
     PolarGrid *Density, *Vrad, *Vtheta, *Label;
{
  MPI_Request req1, req2, req3, req4;
  int l, prev, next, oo, o, nr;
  if (!allocated_com) AllocateComm ();
  prev = CPU_Rank-1;
  next = CPU_Rank+1;
  l = CPUOVERLAP*NSEC;
  nr = Density->Nrad;
  oo = (nr-CPUOVERLAP)*NSEC;
  o = (nr-2*CPUOVERLAP)*NSEC;
  memcpy (SendInnerBoundary, Density->Field+l, l*sizeof(real));
  memcpy (SendInnerBoundary+l, Vrad->Field+l, l*sizeof(real));
  memcpy (SendInnerBoundary+2*l, Vtheta->Field+l, l*sizeof(real));
  memcpy (SendOuterBoundary, Density->Field+o, l*sizeof(real));
  memcpy (SendOuterBoundary+l, Vrad->Field+o, l*sizeof(real));
  memcpy (SendOuterBoundary+2*l, Vtheta->Field+o, l*sizeof(real));
  if (AdvecteLabel == YES) {
    memcpy (SendInnerBoundary+3*l, Label->Field+l, l*sizeof(real));
    memcpy (SendOuterBoundary+3*l, Label->Field+o, l*sizeof(real));
  }
  if (CPU_Rank % 2 == 0) {
    if (CPU_Rank > 0) {
      MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req2);
    }
    if (CPU_Rank < CPU_Number-1) {
      MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req3);
      MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req4);
    }
  } else {
    if (CPU_Rank < CPU_Number-1) {
      MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req3);
      MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req4);
    }
    if (CPU_Rank > 0) {
      MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req2);
    }
  }
  if (CPU_Rank > 0) {
    MPI_Wait (&req1, &stat);
    MPI_Wait (&req2, &stat);
    memcpy (Density->Field, RecvInnerBoundary, l*sizeof(real));
    memcpy (Vrad->Field, RecvInnerBoundary+l, l*sizeof(real));
    memcpy (Vtheta->Field, RecvInnerBoundary+2*l, l*sizeof(real));
    if (AdvecteLabel == YES) {
      memcpy (Label->Field, RecvInnerBoundary+3*l, l*sizeof(real));
    }
  }
  if (CPU_Rank < CPU_Number-1) {
    MPI_Wait (&req3, &stat);
    MPI_Wait (&req4, &stat);
    memcpy (Density->Field+oo, RecvOuterBoundary, l*sizeof(real));
    memcpy (Vrad->Field+oo, RecvOuterBoundary+l, l*sizeof(real));
    memcpy (Vtheta->Field+oo, RecvOuterBoundary+2*l, l*sizeof(real));
    if (AdvecteLabel == YES) {
      memcpy (Label->Field+oo, RecvOuterBoundary+3*l, l*sizeof(real));
    }
  }
}
