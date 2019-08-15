#include "mp.h"

extern boolean OpenInner, NonReflecting;
extern Pair DiskOnPrimaryAcceleration;

Force *AllocateForce (dimfxy)
int dimfxy;
{
  int i;
  Force *force;
  real *globalforce;
  force  = (Force *) prs_malloc (sizeof(Force));
  globalforce = (real *) prs_malloc (sizeof(real) * 4 * dimfxy);
  for (i = 0; i < 4*dimfxy; i++)
    globalforce[i] = 0.;
  force->GlobalForce = globalforce;
  return force;
}

void FreeForce (force)
     Force *force;
{
  free (force->GlobalForce);
}

void ComputeForce (force, Rho, x, y, rsmoothing, mass, dimfxy)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
     int dimfxy;
{
  int i, j, l, ns, k;
  real xc, yc, cellmass, dx, dy, distance, d2, dist2, rh, a;
  real InvDist3, hill_cut, hillcutfactor;
  real *fxi, *fxo, *fyi, *fyo;
  real *localforce, *globalforce;
  real *dens, *abs, *ord;
  fxi = (real *) prs_malloc (sizeof(real) * dimfxy);
  fxo = (real *) prs_malloc (sizeof(real) * dimfxy);
  fyi = (real *) prs_malloc (sizeof(real) * dimfxy);
  fyo = (real *) prs_malloc (sizeof(real) * dimfxy);
  localforce = (real *) prs_malloc (sizeof(real) * 4 * dimfxy);
  globalforce = force->GlobalForce;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  a = sqrt(x*x+y*y);
  rh = pow(mass/3., 1./3.)*a+1e-15;
  for ( k = 0; k < dimfxy; k++ ) {
    fxi[k] = 0.;
    fxo[k] = 0.;
    fyi[k] = 0.;
    fyo[k] = 0.;
  }
  for ( k = 0; k < 4*dimfxy; k++ ) {
    localforce[k] = 0.;
    globalforce[k] = 0.;
  }
  if (FakeSequential && (CPU_Rank > 0)) {
    MPI_Recv (&globalforce[0], 4*dimfxy, MPI_DOUBLE, CPU_Rank-1, 27, MPI_COMM_WORLD, &stat);
    for ( k = 0; k < dimfxy; k++ ) {
      fxi[k] = globalforce [k];
      fxo[k] = globalforce [k + dimfxy];
      fyi[k] = globalforce [k + 2*dimfxy];
      fyo[k] = globalforce [k + 3*dimfxy];
    }
  }
#pragma omp parallel for private(j,hill_cut,cellmass,l,xc,yc,dist2,distance,InvDist3,dx,dy) shared(fxi,fyi,fxhi,fyhi,fxo,fyo,fxho,fyho)
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      xc = abs[l];
      yc = ord[l];
      cellmass = Surf[i]*dens[l];
      dx = xc-x;
      dy = yc-y;
      d2 = dx*dx+dy*dy;
      dist2 = d2 + rsmoothing*rsmoothing;
      distance = sqrt(dist2);
      InvDist3 = 1.0/dist2/distance;
      for ( k = 0; k < dimfxy; k++ ) {
	rh = pow(mass/3., 1./3.)*a+1e-15;
	hillcutfactor = (real)k / (real)(dimfxy-1);
	if ( k != 0 ) {
	  rh *= hillcutfactor;
	  hill_cut = 1.-exp(-d2/(rh*rh));
	}
	else
	  hill_cut = 1.;
	if (Rmed[i] < a) {
#pragma omp atomic
	  fxi[k] += G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	  fyi[k] += G*cellmass*dy*InvDist3*hill_cut;
	} else {
#pragma omp atomic
	  fxo[k] += G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	  fyo[k] += G*cellmass*dy*InvDist3*hill_cut;
	}
      }
    }
  }
  if (FakeSequential) {
    for ( k = 0; k < dimfxy; k++ ) {
      globalforce [k]            = fxi[k];
      globalforce [k + dimfxy]   = fxo[k];
      globalforce [k + 2*dimfxy] = fyi[k];
      globalforce [k + 3*dimfxy] = fyo[k];
    }
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&globalforce[0], 4*dimfxy, MPI_DOUBLE, CPU_Rank+1, 27, MPI_COMM_WORLD);
  } else {
    for ( k = 0; k < dimfxy; k++ ) {
      localforce [k]            = fxi[k];
      localforce [k + dimfxy]   = fxo[k];
      localforce [k + 2*dimfxy] = fyi[k];
      localforce [k + 3*dimfxy] = fyo[k];
    }
    MPI_Allreduce (localforce, globalforce, 4*dimfxy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential)
    MPI_Bcast (globalforce, 4*dimfxy, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  force->fx_inner    = globalforce[0];
  force->fx_ex_inner = globalforce[dimfxy-1];
  force->fx_outer    = globalforce[dimfxy];
  force->fx_ex_outer = globalforce[2*dimfxy-1];
  force->fy_inner    = globalforce[2*dimfxy];
  force->fy_ex_inner = globalforce[3*dimfxy-1];
  force->fy_outer    = globalforce[3*dimfxy];
  force->fy_ex_outer = globalforce[4*dimfxy-1];
  force->GlobalForce = globalforce;
  free (localforce);
  free (fxi);
  free (fxo);
  free (fyi);
  free (fyo);
}

real compute_smoothing (r)
     real r;
{
  real smooth;
  smooth = THICKNESSSMOOTHING * AspectRatio(r) * pow(r, 1.0+FLARINGINDEX);
  return smooth;
}


void UpdateLog (fc, psys, Rho, Energy, outputnb, time, dimfxy)
     Force *fc;
     PolarGrid *Rho, *Energy;
     PlanetarySystem *psys;
     int outputnb;
     int dimfxy;
     real time;
{
  extern boolean SelfGravity, SGZeroMode;
  int i, nb, k;
  real x, y, r, m, vx, vy, smoothing;
  real *globalforce;
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
      smoothing = compute_smoothing(r);
    ComputeForce (fc, Rho, x, y, smoothing, m, dimfxy);
    globalforce = fc->GlobalForce;
    if (CPU_Rank == CPU_Number-1) {
      sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
      out = fopen (filename, "a");
      if (out == NULL) {
	fprintf (stderr, "Can't open %s\n", filename);
	fprintf (stderr, "Aborted.\n");
      }
      fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb, \
	       x*fc->fy_inner-y*fc->fx_inner,				\
	       x*fc->fy_outer-y*fc->fx_outer,				\
	       x*fc->fy_ex_inner-y*fc->fx_ex_inner,			\
	       x*fc->fy_ex_outer-y*fc->fx_ex_outer,			\
	       vx*fc->fx_inner+vy*fc->fy_inner,				\
	       vx*fc->fx_outer+vy*fc->fy_outer,				\
	       vx*fc->fx_ex_inner+vy*fc->fy_ex_inner,			\
	       vx*fc->fx_ex_outer+vy*fc->fy_ex_outer, time);
      fclose (out);
      if ( SGZeroMode || !SelfGravity ) {
	for ( k = 0; k < dimfxy; k++ ) {
	  sprintf (filename, "%stqwk%d_%d.dat", OUTPUTDIR, i, k);
	  out = fopen (filename, "a");
	  if (out == NULL) {
	    fprintf (stderr, "Can't open %s\n", filename);
	    fprintf (stderr, "Aborted.\n");
	  }
	  fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb, \
		   x*globalforce[2*dimfxy+k]-y*globalforce[k],		\
		   x*globalforce[3*dimfxy+k]-y*globalforce[dimfxy+k],	\
		   vx*globalforce[k]+vy*globalforce[2*dimfxy+k],	\
		   vx*globalforce[dimfxy+k]+vy*globalforce[3*dimfxy+k], time);
	  fclose (out);
	}
      }
    }
  }
}
