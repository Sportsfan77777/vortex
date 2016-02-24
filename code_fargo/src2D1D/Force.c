#include "mp.h"

extern boolean OpenInner, NonReflecting;
real pfermi=0.6;

Force ComputeForce (Rho, x, y, rsmoothing, mass)
PolarGrid *Rho;
real x, y, rsmoothing, mass;
{
  int i, j, l, ns;
  real localforce[8]={0.,0.,0.,0.,0.,0.,0.,0.}, globalforce[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  real xc, yc, cellmass, dx, dy, distance, dist2, rh, a;
  real InvDist3, fxi, fyi, fxhi, fyhi, fxo, fyo, fxho, fyho, hill_cut;
  real *dens, *abs, *ord;
  Force Force;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  fxi = fyi = fxhi = fyhi = fxo = fyo = fxho = fyho = 0.0;
  a = sqrt(x*x+y*y);
  rh = pow((1./3.*mass), 1./3.)*a;
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
#pragma omp parallel for private(j,hill_cut,cellmass,l,xc,yc,dist2,distance,InvDist3,dx,dy) shared(fxi,fyi,fxhi,fyhi,fxo,fyo,fxho,fyho)
  //  for (i = Zero_or_active; i < Max_or_active; i++) {
  for (i = CPUOVERLAP; i < NRAD-CPUOVERLAP; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      xc = abs[l];
      yc = ord[l];
      cellmass = Surf[i]*dens[l];
      dx = xc-x;
      dy = yc-y;
      dist2 = dx*dx+dy*dy;
      if (rh == 0.)
	hill_cut=0.;
      else
	hill_cut = 1/(exp(-(sqrt(dist2)/rh-pfermi)/pfermi*10)+1);
        /* In the standard version, a Gaussian function is used.
	 * We use a fermi function, of parameter pfermi, set above.
         * (see also Gframeforce.c) */
      dist2 += rsmoothing*rsmoothing;
      distance = sqrt(dist2);
      InvDist3 = 1.0/dist2/distance;
      if (Rmed[i] < a) {
#pragma omp atomic
	fxi += G*cellmass*dx*InvDist3;
#pragma omp atomic
	fyi += G*cellmass*dy*InvDist3;
#pragma omp atomic
	fxhi+= G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	fyhi+= G*cellmass*dy*InvDist3*hill_cut;
      } else {
#pragma omp atomic
	fxo += G*cellmass*dx*InvDist3;
#pragma omp atomic
	fyo += G*cellmass*dy*InvDist3;
#pragma omp atomic
	fxho+= G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	fyho+= G*cellmass*dy*InvDist3*hill_cut;
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

void UpdateLog (psys, Rho, outputnb, time)
     PolarGrid *Rho;
     PlanetarySystem *psys;
     int outputnb;
     real time;
{
  int i, ii, nb;
  real x, y, PlanetDistance, rplanet, m, RRoche, vx, vy, smoothing, iplanet, cs, frac;
  Force fc;
  FILE *out;
  char filename[MAX1D];
  nb = psys->nb;
  for (i = 0; i < nb; i++) {
    x = psys->x[i];
    y = psys->y[i];
    vx = psys->vx[i];
    vy = psys->vy[i];
    m = psys->mass[i];
    if (RocheSmoothing) {
      PlanetDistance = sqrt((x-psys->x[0])*(x-psys->x[0])+(y-psys->y[0])*(y-psys->y[0]));
      RRoche = PlanetDistance*pow((1.0/3.0*m),1.0/3.0);
      smoothing = RRoche*ROCHESMOOTHING;
    } else {
    if (i==0) rplanet=0.;
    else      rplanet = sqrt(x*x+y*y);
      iplanet = GetGlobalIFrac (rplanet);
      frac = iplanet-floor(iplanet);
      ii = (int)iplanet;
      cs = GLOBAL_SOUNDSPEED[ii]*(1.0-frac)+\
	GLOBAL_SOUNDSPEED[ii+1]*frac;
      smoothing = cs * rplanet * sqrt(rplanet) * THICKNESSSMOOTHING;
    }
    fc = ComputeForce (Rho, x, y, smoothing, m);
    if (CPU_Rank == CPU_Number-1) {
      sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
      out = fopen (filename, "a");
      if (out == NULL) {
	fprintf (stderr, "Can't open %s\n", filename);
	fprintf (stderr, "Aborted.\n");
      }
      fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb,\
	       x*fc.fy_inner-y*fc.fx_inner,\
	       x*fc.fy_outer-y*fc.fx_outer,\
	       x*fc.fy_ex_inner-y*fc.fx_ex_inner,\
	       x*fc.fy_ex_outer-y*fc.fx_ex_outer,\
	       vx*fc.fx_inner+vy*fc.fy_inner,\
	       vx*fc.fx_outer+vy*fc.fy_outer,\
	       vx*fc.fx_ex_inner+vy*fc.fy_ex_inner,\
	       vx*fc.fx_ex_outer+vy*fc.fy_ex_outer, time);
      fclose (out);
    }
  }
}
