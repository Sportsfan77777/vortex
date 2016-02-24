/* This file replaces the Pframeforce.c of the standard version. *
 * The forces are computed in the frame of the barycenter.       */


#include "mp.h"

extern boolean AllowAccretion, Corotating;
extern real pfermi;
static real q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static real vt_int[MAX1D], vt_cent[MAX1D];

/* The following procedure ensures that the centre of mass
 * of the system (star + planets + gas) remains at r=0, v=0.
 */
void BarycenterConservation(sys,gasdens,gasvrad,gasvtheta)
PlanetarySystem *sys;
PolarGrid *gasdens, *gasvrad, *gasvtheta;
{
  real *dens, *vrad, *vtheta, *abs, *ord;
  real m_disk,x_disk,y_disk,vx_disk,vy_disk;
  real loc_m_disk,loc_x_disk,loc_y_disk,loc_vx_disk,loc_vy_disk;
  real m_sys,x_sys,y_sys,vx_sys,vy_sys;
  real dx, dy, dvx, dvy;
  real cellmass,mp,vx,vy;
  int i,j,k,l,ljp,lip,ns,n;
  ns = gasdens->Nsec;
  n = sys->nb;
  dens = gasdens->Field;
  vrad = gasvrad->Field;
  vtheta = gasvtheta->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  loc_m_disk = loc_x_disk = loc_y_disk = loc_vx_disk = loc_vy_disk = 0.;
  m_sys = x_sys = y_sys = vx_sys = vy_sys = 0.;
  if (n > 1) {
    for (i = CPUOVERLAP; i < NRAD-CPUOVERLAP; i++) { /* Valid in MPI Version as well... */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lip = l+ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	cellmass = dens[l]*Surf[i];
	loc_m_disk += cellmass;
	loc_x_disk += abs[l]*cellmass;
	loc_y_disk += ord[l]*cellmass;
	vx = (vrad[l]+vrad[lip])/2.*Cosinus[j];
	vx -= (vtheta[l]+vtheta[ljp])/2.*Sinus[j];
	vy = (vrad[l]+vrad[lip])/2.*Sinus[j];
	vy += (vtheta[l]+vtheta[ljp])/2.*Cosinus[j];
	loc_vx_disk += vx*cellmass;
	loc_vy_disk += vy*cellmass;
      }
    }
    MPI_Allreduce (&loc_m_disk,&m_disk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (&loc_x_disk,&x_disk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (&loc_y_disk,&y_disk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (&loc_vx_disk,&vx_disk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce (&loc_vy_disk,&vy_disk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (k = 0; k < n; k++) {
      mp = sys->mass[k];
      m_sys += mp;
      x_sys += sys->x[k]*mp;
      y_sys += sys->y[k]*mp;
      vx_sys += sys->vx[k]*mp;
      vy_sys += sys->vy[k]*mp;
    }
    dx = x_disk+x_sys;
    dy = y_disk+y_sys;
    dvx = vx_disk+vx_sys;
    dvy = vy_disk+vy_sys;
    for (k = 0; k < n; k++) {
      sys->x[k] -= dx/m_sys;
      sys->y[k] -= dy/m_sys;
      sys->vx[k] -= dvx/m_sys;
      sys->vy[k] -= dvy/m_sys;
    }
  }
  if (n == 1)
    sys->x[0] = sys->y[0] = sys->vx[0] = sys->vy[0] = 0.;
}


void FillForcesArrays (sys, gasdens)
PlanetarySystem *sys;
PolarGrid *gasdens;
{
  int i,ii,j,l,nr,ns,k,NbPlanets;
  real x, y, distance, distancesmooth, iplanet;
  real xplanet, yplanet, rplanet, RRoche,smooth, mplanet, frac;
  real PlanetDistance, *Pot, pot, smoothing, cs;
  real *abs, *ord;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  Pot = Potential->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  NbPlanets = sys->nb;
#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) Pot[i] = 0.0;
  for (k = 0; k < NbPlanets; k++) {
    xplanet = sys->x[k];
    yplanet = sys->y[k];
    mplanet = sys->mass[k];
    if (k>0) mplanet *= MassTaper;
    if (RocheSmoothing) {
      PlanetDistance = sqrt((xplanet-sys->x[0])*(xplanet-sys->x[0])+(yplanet-sys->y[0])*(yplanet-sys->y[0]));
      RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
      smoothing = RRoche*ROCHESMOOTHING;
    } else {
      if (k==0) rplanet=0.;
      else      rplanet = sqrt(xplanet*xplanet+yplanet*yplanet);
      iplanet = GetGlobalIFrac (rplanet);
      frac = iplanet-floor(iplanet);
      ii = (int)iplanet;
      cs = GLOBAL_SOUNDSPEED[ii]*(1.0-frac)+\
	GLOBAL_SOUNDSPEED[ii+1]*frac;
      smoothing = cs * rplanet * sqrt(rplanet) * THICKNESSSMOOTHING;
    }
    smooth = smoothing*smoothing;
#pragma omp parallel for private(j,l,x,y,distance,distancesmooth,pot)
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	x = abs[l];
	y = ord[l];
	distance = (x-xplanet)*(x-xplanet)+(y-yplanet)*(y-yplanet);
	distancesmooth = sqrt(distance+smooth);
	pot = -G*mplanet/distancesmooth;
	Pot[l] += pot;
      }
    }
  }
}


void FillPotplanet (sys, k, gasdens)
PlanetarySystem *sys;
int k;
PolarGrid *gasdens;
{
  int i,ii,j,l,nr,ns;
  real x, y, distance, distancesmooth, iplanet;
  real xplanet, yplanet, RRoche,smooth, mplanet, frac;
  real PlanetDistance, rplanet, *Potp, pot, smoothing, cs;
  real *abs, *ord;
  real xstar, ystar;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  Potp= Potplanet->Field;
  nr = Potplanet->Nrad;
  ns = Potplanet->Nsec;
  xplanet = sys->x[k];
  yplanet = sys->y[k];
  mplanet = sys->mass[k];
  xstar=sys->x[0];
  ystar=sys->y[0];
  if (RocheSmoothing) {
    PlanetDistance = sqrt( (xplanet-xstar)*(xplanet-xstar)+\
			   (yplanet-ystar)*(yplanet-ystar) );
    RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
    smoothing = RRoche*ROCHESMOOTHING;
  } else {
    if (k==0) rplanet=0.;
    else      rplanet = sqrt(xplanet*xplanet+yplanet*yplanet);
    iplanet = GetGlobalIFrac (rplanet);
    frac = iplanet-floor(iplanet);
    ii = (int)iplanet;
    cs = GLOBAL_SOUNDSPEED[ii]*(1.0-frac)+\
         GLOBAL_SOUNDSPEED[ii+1]*frac;
    smoothing = cs * rplanet * sqrt(rplanet) * THICKNESSSMOOTHING;
  }
  smooth = smoothing*smoothing;

#pragma omp parallel for private(j,l,x,y,distance,distancesmooth,pot)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      x = abs[l];
      y = ord[l];
      distance = (x-xplanet)*(x-xplanet)+(y-yplanet)*(y-yplanet);
      distancesmooth = sqrt(distance+smooth);
      pot = -G*mplanet/distancesmooth;
      Potp[l] = pot;
    }
  }
}


Pair FeltByDisk (Rho, xp, yp, dt, Mp)
PolarGrid *Rho;
real xp, yp, dt, Mp;		       
{
  int i,j,l,ljm,ljp,lim,lip,nr,ns;
  real *dens, *Potp, *abs, *ord;
  real thetap,costhetap,sinthetap,cosu,sinu;
  real dP_rp,dH,dP_tp,dvt,dvr,dvrp,dvtp,cellmass;
  real dxtheta,invdxtheta;
  real xcell,ycell,dist2,RRoche,PlanetDistance,excl;
  Pair dHdPrp;
  Potp = Potplanet->Field;
  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  thetap = atan2(yp,xp);
  costhetap = cos(thetap);
  sinthetap = sin(thetap);
  dH = dP_rp = dP_tp = 0.;
  if (ExcludeHill) {
    abs = CellAbscissa->Field;
    ord = CellOrdinate->Field;
    PlanetDistance = sqrt(xp*xp+yp*yp);    // Could be improved in taking a_p.
    RRoche = PlanetDistance*pow((1.0/3.0*Mp),1.0/3.0);
  }

  //  for (i = 1; i < NRAD-CPUOVERLAP; i++) { // OK even in parallel version
  for (i = CPUOVERLAP; i < NRAD-CPUOVERLAP; i++) { // OK even in parallel version
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    invdxtheta = 1.0/dxtheta;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lim = l-ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      lip = l+ns;
      lim = l-ns;
      cellmass = dens[l]*Surf[i];
      if (ExcludeHill) {
	/* In the standard version, a Gaussian function is used.
	 * We use a fermi function, of parameter pfermi,
         * set in Force.c. */
	xcell = abs[l];
	ycell = ord[l];
	dist2 = (xp-xcell)*(xp-xcell)+(yp-ycell)*(yp-ycell);
	excl  = 1/(exp(-(sqrt(dist2)/RRoche-pfermi)/pfermi*10)+1);
	cellmass *= excl;
      }
      cosu = Cosinus[j]*costhetap + Sinus[j]*sinthetap;
      sinu = Sinus[j]*costhetap - Cosinus[j]*sinthetap;
      dvt = -(Potp[ljp]-Potp[ljm])*invdxtheta/2.*dt;
      dvr = -0.5*( (Potp[l]-Potp[lim])*InvDiffRmed[i]
                 + (Potp[lip]-Potp[l])*InvDiffRmed[i+1] )*dt;
      dvrp = cosu*dvr - sinu*dvt;
      dvtp = sinu*dvr + cosu*dvt;
      dH += Rmed[i]*dvt*cellmass;
      dP_rp += cellmass*dvrp;
      dP_tp += cellmass*dvtp;
    }
  }
  dHdPrp.x = dH;
  dHdPrp.y = dP_rp;
  return dHdPrp;
}


void AdvanceSystemFromDisk (Rho, sys, dt)
PolarGrid *Rho;
PlanetarySystem *sys;
real dt;		       
{
  int NbPlanets, k;
  Pair dHdP;
  Pair dhdp;
  real x,y,rp,m,dvtp,dvrp,thetap,costhetap,sinthetap;
  NbPlanets = sys->nb;
  for (k = 0; k < NbPlanets; k++) {
    if (sys->FeelDisk[k] == YES) {
      m=sys->mass[k];
      x=sys->x[k];
      y=sys->y[k];
      FillPotplanet(sys, k, Rho);
      dhdp = FeltByDisk (Rho, x, y, dt, m);
      MPI_Allreduce (&(dhdp.x),&(dHdP.x),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce (&(dhdp.y),&(dHdP.y),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      rp=sqrt(x*x+y*y);
      thetap = atan2(y,x);
      costhetap = cos(thetap);
      sinthetap = sin(thetap);
      dvrp = -dHdP.y/m;
      dvtp = -dHdP.x/m/rp;
      sys->vx[k] += dvrp*costhetap - dvtp*sinthetap;
      sys->vy[k] += dvrp*sinthetap + dvtp*costhetap;
    }
  }
}



void AdvanceSystemRK5 (sys, dt)
PlanetarySystem *sys;
real dt;
{
  int i, n;
  boolean *feelothers;
  real dtheta, omega, rdot, x, y, r, new_r, vx, vy, theta, denom;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    q0[i] = sys->x[i];
    q0[i+n] = sys->y[i];
    q0[i+2*n] = sys->vx[i];
    q0[i+3*n] = sys->vy[i];
    PlanetMasses[i] = sys->mass[i];
  }
  feelothers = sys->FeelOthers;
  RungeKunta (q0, dt, PlanetMasses, q1, n, feelothers);
  for (i = 0; i < n; i++) {
    sys->x[i] = q1[i];
    sys->y[i] = q1[i+n];
    sys->vx[i] = q1[i+2*n];
    sys->vy[i] = q1[i+3*n];
    if ( (i==0) && (PhysicalTime < RELEASEDATE)) i++;
  }
  /* The RELEASE stuff concerns planet number 1 (#0 is the star) */
  if (PhysicalTime < RELEASEDATE) {
    x = sys->x[1];
    y = sys->y[1];
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    rdot = (RELEASERADIUS-r)/(RELEASEDATE-PhysicalTime);
    omega = sqrt((1.+sys->mass[1])/r/r/r);
    new_r = r + rdot*dt;
    denom = r-new_r;
    if (denom != 0.0) {
      dtheta = 2.*dt*r*omega/denom*(sqrt(r/new_r)-1.);
    } else {
      dtheta = omega*dt;
    }
    vx = rdot;
    vy = new_r*sqrt((sys->mass[0]+sys->mass[1])/new_r/new_r/new_r);
    sys->x[1] = new_r*cos(dtheta+theta);
    sys->y[1] = new_r*sin(dtheta+theta);
    sys->vx[1]= vx*cos(dtheta+theta)-vy*sin(dtheta+theta); 
    sys->vy[1]= vx*sin(dtheta+theta)+vy*cos(dtheta+theta); 
  }
}


void SolveOrbits (sys)
PlanetarySystem *sys;
{
  int i, n;
  real x,y,vx,vy;
  real xs,ys,vxs,vys,m;
  n = sys->nb;
  xs = sys->x[0];
  ys = sys->y[0];
  vxs = sys->vx[0];
  vys = sys->vy[0];
  for (i = 1; i < n; i++) {
    x = sys->x[i]-xs;        // use heliocentric coordinates.
    y = sys->y[i]-ys;        // use heliocentric coordinates.
    vx = sys->vx[i]-vxs;     // use heliocentric coordinates.
    vy = sys->vy[i]-vys;     // use heliocentric coordinates.
    m = sys->mass[0]+sys->mass[i];
    FindOrbitalElements (x,y,vx,vy,m,i);
  }
} 

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }
  return lapl;
}

void
InitGas (Rho, Vr, Vt)
PolarGrid *Rho, *Vr, *Vt;
{
  int i, j, l, nr, ns;
  real *dens, *vr, *vt;
  float temporary;
  FILE *CS;
  char csfile[512];
  real  r, rg, omega, ri, viscosity, t1, t2, r1, r2;
  dens= Rho->Field;
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Rho->Nrad;
  ns  = Rho->Nsec;
  sprintf (csfile, "%s%s", OUTPUTDIR, "soundspeed.dat");
  CS = fopen (csfile, "r");
  if (CS == NULL) {
    for (i = 0; i < nr; i++) {
      SOUNDSPEED[i] = ASPECTRATIO * sqrt(G*1.0/Rmed[i]) * pow(Rmed[i], FLARINGINDEX);
    }
    for (i = 0; i < GLOBALNRAD; i++) {
      GLOBAL_SOUNDSPEED[i] = ASPECTRATIO * sqrt(G*1.0/GlobalRmed[i]) * pow(GlobalRmed[i], FLARINGINDEX);
    }
  } else {
    masterprint ("Reading soundspeed.dat file\n");
    for (i = 0; i < GLOBALNRAD; i++) {
      fscanf (CS, "%f", &temporary);
      GLOBAL_SOUNDSPEED[i] = (real)temporary;
    }
    for (i = 0; i < nr; i++) {
      SOUNDSPEED[i] = GLOBAL_SOUNDSPEED[i+IMIN];
    }
  }
  for (i = 1; i < GLOBALNRAD; i++) {
    vt_int[i]=(GLOBAL_SOUNDSPEED[i]*GLOBAL_SOUNDSPEED[i]*Sigma(GlobalRmed[i])-\
	       GLOBAL_SOUNDSPEED[i-1]*GLOBAL_SOUNDSPEED[i-1]*Sigma(GlobalRmed[i-1]))/\
      (.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1])+\
      G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
  }
  t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
  r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
  t2 = vt_cent[0];
  r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  t1 = t1-r1/(r2-r1)*(t2-t1);
  vt_cent[0] = t1;
  ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  vt_cent[GLOBALNRAD]=vt_cent[GLOBALNRAD-1];
  for (i = 0; i <= nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    }
    else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    viscosity = VISCOSITY;
    if (ViscosityAlpha)
      viscosity = ALPHAVISCOSITY*ASPECTRATIO*ASPECTRATIO*pow(r,.5+2.0*FLARINGINDEX);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rg = r;
      omega = sqrt(G*1.0/rg/rg/rg);
      vt[l] = omega*r*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(r,2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      vt[l] -= OmegaFrame*r;
      if (CentrifugalBalance)
	vt[l] = vt_cent[i+IMIN];
      if (i == nr) 
	vr[l] = 0.0;
      else {
	vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
	if (ViscosityAlpha) {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
	} else {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
	}
      }
      dens[l] = SigmaMed[i];
    }
  }
  for (j = 0; j < ns; j++)
    vr[j] = vr[j+ns*nr] = 0.0;
}
