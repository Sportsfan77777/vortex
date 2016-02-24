#include "mp.h"

static real k1[MAX1D], k2[MAX1D], k3[MAX1D], k4[MAX1D], k5[MAX1D], k6[MAX1D];


void DerivMotionRK5 (q_init, masses, deriv, n, dt, feelothers)
real *q_init, *deriv, *masses, dt;
boolean *feelothers;
int n;
{
  real *x,*y,*vx,*vy, dist;
  real *derivx, *derivy, *derivvx, *derivvy;
  int i, j;
  x = q_init;
  y = x+n;
  vx = y+n;
  vy = vx+n;
  derivx = deriv;
  derivy = deriv+n;
  derivvx = derivy+n;
  derivvy = derivvx+n;

  for (i = 0; i < n; i++) {
    derivx[i] = vx[i];
    derivy[i] = vy[i];
    derivvx[i]= 0.;
    derivvy[i]= 0.;
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if ((j != i) && ( (feelothers[i] == YES) || (j == 0) ) ) {
	dist = (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]);
	dist = sqrt(dist);
	derivvx[i] += G*masses[j]/dist/dist/dist*(x[j]-x[i]);
	derivvy[i] += G*masses[j]/dist/dist/dist*(y[j]-y[i]);
      }
    }
  }
  for (i = 0; i < 4*n; i++)
    deriv[i] *= dt;
}


void TranslatePlanetRK5 (qold, c1, c2, c3, c4, c5, qnew, n)
real *qold, *qnew;
real c1, c2, c3, c4, c5;
int n;
{
  int i;
  for (i = 0; i < 4*n; i++)
    qnew[i] = qold[i]+c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i];
}

void RungeKunta (q0, dt, masses, q1, n, feelothers)
real *q0, *q1;
real dt, *masses;
boolean *feelothers;
int n;
{
  int i;
  real timestep;
  timestep = dt;
  for (i = 0; i < n*4; i++) {
    k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = k6[i];
  }
  DerivMotionRK5 (q0, masses, k1, n, timestep, feelothers);
  TranslatePlanetRK5 (q0, 0.2, 0.0, 0.0, 0.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k2, n, timestep, feelothers);
  TranslatePlanetRK5 (q0, 0.075, 0.225, 0.0, 0.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k3, n, timestep, feelothers);
  TranslatePlanetRK5 (q0, 0.3, -0.9, 1.2, 0.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k4, n, timestep, feelothers);
  TranslatePlanetRK5 (q0, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0, q1, n);
  DerivMotionRK5 (q1, masses, k5, n, timestep, feelothers);
  TranslatePlanetRK5 (q0, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, q1, n);
  DerivMotionRK5 (q1, masses, k6, n, timestep, feelothers);
  for (i = 0; i < 4*n; i++) {
    q1[i]=q0[i]+37.0/378.0*k1[i]+250.0/621.0*k3[i]+125.0/594.0*k4[i]+512.0/1771.0*k6[i];
  }
}
 
