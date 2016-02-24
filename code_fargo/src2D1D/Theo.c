#include "mp.h"

extern real      ScalingFactor;


real Sigma(r)
real r;
{
  real s;
  s = ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
  return s;
}

/***************************************************************
 * This version of the Sigma procedure reads the density profile
 * in the file 'sigma_init.dat' located in the output directory.
 * The file must be in ascii in two columns : r Sigma,
 * r should be increasing from one line to the next.
 * A linear fit is performed between the data points.
 * The density is then multiplied by the ScalingFactor.

real Sigma(r)
     real r;
{
  FILE * sigmafile;
  char sigmafilename[100];
  float rr,ss,rrp,ssp;
  real lambda,S;
  int i,NN;
  NN=16384;

  sprintf(sigmafilename,"%ssigma_init.dat",OUTPUTDIR);
  sigmafile=fopen(sigmafilename,"r");
  if (sigmafile == NULL) {
    mastererr ("Can't read %s file. Exit.\n",sigmafilename);
    prs_exit (1);
  }
  rr=ss=rrp=ssp=0.;
  for(i=0;(rrp<r && i<=NN);i++) {
    rr=rrp;
    ss=ssp;
    fscanf(sigmafile,"%f %f",&rrp,&ssp);
  }
  fclose(sigmafile);
  if(rr>r)
    fprintf(stderr,"Erreur sur sigmainit rr=%.3f > r=%.3f (i=%d)\n",rr,r,i);
  if(rrp<r)
    fprintf(stderr,"Erreur sur sigmainit rrp=%.3f < r=%.3f (i=%d)\n",rrp,r,i);
  lambda=(r-rr)/(rrp-rr);
  S=lambda*ssp+(1.-lambda)*ss;
  return S*SIGMA0;
}
 */


void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}
