#include "mp.h"

extern boolean OpenInner, NonReflecting;
extern Pair DiskOnPrimaryAcceleration;

Force ComputeForceStockholm (Rho, x, y, rsmoothing, mass)
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
{
  int i, j, l, ns;
  real localforce[8]={0.,0.,0.,0.,0.,0.,0.,0.}, globalforce[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  real xc, yc, cellmass, dx, dy, distance, dist2, rh, a;
  real InvDist3, fxi, fyi, fxhi, fyhi, fxo, fyo, fxho, fyho, outside_hill, inside_hill;
  real *dens, *abs, *ord;
  Force Force;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  fxi = fyi = fxhi = fyhi = fxo = fyo = fxho = fyho = 0.0;
  a = sqrt(x*x+y*y);
  rh = pow(mass/3., 1./3.)*a+1e-15;
  if (FakeSequential && (CPU_Rank > 0)) {
    MPI_Recv (&globalforce, 8, MPI_DOUBLE, CPU_Rank-1, 27, MPI_COMM_WORLD, &stat);
    fxi = globalforce[0];
    fyi = globalforce[1];
    fxhi= globalforce[2];
    fyhi= globalforce[3];
    fxo = globalforce[4];
    fyo = globalforce[5];
    fxho= globalforce[6];
    fyho= globalforce[7];
  }
#pragma omp parallel for private(j,outside_hill,inside_hill,cellmass,l,xc,yc,dist2,distance,InvDist3,dx,dy) shared(fxi,fyi,fxhi,fyhi,fxo,fyo,fxho,fyho)
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      xc = abs[l];
      yc = ord[l];
      cellmass = Surf[i]*dens[l];
      dx = xc-x;
      dy = yc-y;
      dist2 = dx*dx+dy*dy;
      outside_hill = (dist2 >= rh*rh ? 1.0 : 0.0);
      inside_hill = (((dist2 >= 0.25*rh*rh) && (dist2 < rh*rh)) ? 1.0 : 0.0);
      dist2 += rsmoothing*rsmoothing;
      distance = sqrt(dist2);
      InvDist3 = 1.0/dist2/distance;
      if (Rmed[i] < a) {
#pragma omp atomic
	fxhi += G*cellmass*dx*InvDist3*inside_hill;
#pragma omp atomic
	fyhi += G*cellmass*dy*InvDist3*inside_hill;
#pragma omp atomic
	fxi+= G*cellmass*dx*InvDist3*outside_hill;
#pragma omp atomic
	fyi+= G*cellmass*dy*InvDist3*outside_hill;
      } else {
#pragma omp atomic
	fxho += G*cellmass*dx*InvDist3*inside_hill;
#pragma omp atomic
	fyho += G*cellmass*dy*InvDist3*inside_hill;
#pragma omp atomic
	fxo+= G*cellmass*dx*InvDist3*outside_hill;
#pragma omp atomic
	fyo+= G*cellmass*dy*InvDist3*outside_hill;
      }
    }
  }
  if (FakeSequential) {
    globalforce[0]=fxi;
    globalforce[1]=fyi;
    globalforce[2]=fxhi;
    globalforce[3]=fyhi;
    globalforce[4]=fxo;
    globalforce[5]=fyo;
    globalforce[6]=fxho;
    globalforce[7]=fyho;
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&globalforce, 8, MPI_DOUBLE, CPU_Rank+1, 27, MPI_COMM_WORLD);
  } else {
    localforce[0] = fxi;
    localforce[1] = fyi;
    localforce[2] = fxhi;
    localforce[3] = fyhi;
    localforce[4] = fxo;
    localforce[5] = fyo;
    localforce[6] = fxho;
    localforce[7] = fyho;
    MPI_Allreduce (&localforce, &globalforce, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential)
    MPI_Bcast (&globalforce, 8, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  Force.fx_inner    = globalforce[0];
  Force.fy_inner    = globalforce[1];
  Force.fx_ex_inner = globalforce[2];
  Force.fy_ex_inner = globalforce[3];
  Force.fx_outer    = globalforce[4];
  Force.fy_outer    = globalforce[5];
  Force.fx_ex_outer = globalforce[6];
  Force.fy_ex_outer = globalforce[7];
  return Force;
}

Pair MassInOut (Rho, a)
PolarGrid *Rho;
real a;
{
  int i, j, l, ns;
  real localmass[2]={0.,0.}, globalmass[2]={0.,0.};
  real massin=0.0, massout=0.0, outside, inside;
  real *dens, cellmass;
  Pair massinout;
  ns = Rho->Nsec;
  dens = Rho->Field;
  if (FakeSequential && (CPU_Rank > 0)) {
    MPI_Recv (&globalmass, 2, MPI_DOUBLE, CPU_Rank-1, 38, MPI_COMM_WORLD, &stat);
    massin  = globalmass[0];
    massout = globalmass[1];
  }
#pragma omp parallel for private(j,cellmass,inside,outside,l) shared(massin, massout)
  for (i = Zero_or_active; i < Max_or_active; i++) {
    if (Rinf[i] > a) {
      outside = 1.0;
      inside = 0.0;
    }
    if (Rsup[i] < a) {
      outside = 0.0;
      inside = 1.0;
    }
    if ((Rinf[i] < a) && (Rsup[i] > a)) {
      outside = (Rsup[i]*Rsup[i]-a*a)/(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
      inside = (-Rinf[i]*Rinf[i]+a*a)/(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    }
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      cellmass = Surf[i]*dens[l];
#pragma omp atomic
      massin += cellmass*inside;
#pragma omp atomic
      massout += cellmass*outside;
    }
  }
  if (FakeSequential) {
    globalmass[0]=massin;
    globalmass[1]=massout;
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&globalmass, 2, MPI_DOUBLE, CPU_Rank+1, 38, MPI_COMM_WORLD);
  } else {
    localmass[0] = massin;
    localmass[1] = massout;
    MPI_Allreduce (&localmass, &globalmass, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential)
    MPI_Bcast (&globalmass, 2, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  massinout.x    = globalmass[0];
  massinout.y    = globalmass[1];
  return massinout;
}

void UpdateLogStockholm (psys, Rho, Energy, outputnb, time)
     PolarGrid *Rho, *Energy;
     PlanetarySystem *psys;
     int outputnb;
     real time;
{
  int i, ii, nb;
  real x, y, r, m, vx, vy, smoothing, iplanet, cs, frac;
  Force fc;
  Pair massinout;
  FILE *out;
  char filename[MAX1D];
  nb = psys->nb;
  for (i = 0; i < nb; i++) {
    x = psys->x[i];
    y = psys->y[i];
    vx = psys->vx[i];
    vy = psys->vy[i];
    r = sqrt(x*x+y*y);
    m = psys->mass[i];
    if (RocheSmoothing)
      smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
    else
      smoothing = compute_smoothing (r);
    fc = ComputeForceStockholm (Rho, x, y, smoothing, m);
    massinout = MassInOut (Rho, r);
    if (CPU_Rank == CPU_Number-1) {
      sprintf (filename, "%storque%d.dat", OUTPUTDIR, i);
      out = fopen (filename, "a");
      if (out == NULL) {
	fprintf (stderr, "Can't open %s\n", filename);
	fprintf (stderr, "Aborted.\n");
      }
      fprintf (out, "%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", time,\
	       massinout.x, massinout.y,\
	       x*fc.fy_inner-y*fc.fx_inner,\
	       x*fc.fy_outer-y*fc.fx_outer,\
	       x*fc.fy_ex_inner-y*fc.fx_ex_inner,\
	       x*fc.fy_ex_outer-y*fc.fx_ex_outer);
      fclose (out);
    }
  }
}
