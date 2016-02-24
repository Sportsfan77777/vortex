#include "mp.h"

static real OldRadii1D[MAX1D], OldRmed1D[MAX1D], New_r[MAX1D];
static int OldNRAD1D;

void ReadPrevDim1D () {
  FILE *DIM, *RAD;
  char name_dim[1024], name_rad[1024];
  int i, foo;
  float Foo;
  double value;
  if (!CPU_Master) return;
  OldNRAD1D = 0;
  sprintf (name_dim, "%sdims1D.dat", OUTPUTDIR);
  sprintf (name_rad, "%sused_rad1D.dat", OUTPUTDIR);
  DIM = fopen (name_dim, "r");
  if (DIM == NULL) return;
  RAD = fopen (name_rad, "r");
  if (RAD == NULL) return;
  fscanf (DIM,"%d %d %d %d %f %d %d\n",&foo,&foo,&foo,&foo,&Foo,&foo,&OldNRAD1D);
  for (i = 0; i <= OldNRAD1D; i++) {
    fscanf (RAD, "%lf", &value);
    OldRadii1D[i] = (real)value;
  }
  fclose (DIM);
  fclose (RAD);
  for (i = 0; i < OldNRAD1D; i++) {
    OldRmed1D[i] = 2.0/3.0*(OldRadii1D[i+1]*OldRadii1D[i+1]*OldRadii1D[i+1]-OldRadii1D[i]*OldRadii1D[i]*OldRadii1D[i]);
    OldRmed1D[i] = OldRmed1D[i] / (OldRadii1D[i+1]*OldRadii1D[i+1]-OldRadii1D[i]*OldRadii1D[i]);
  }
}

void CheckRebin1D (nb) 
int nb;
{
  boolean RebinNeeded = NO, found;
  char radix[1024], filename[1024];
  int i, iold, type, raderr;
  real ifrac, r, dangle;
  FILE *ARR;
  real *Old_r, *OldArray, *NewArray;
  if (!CPU_Master) return;
  if (NRAD1D != OldNRAD1D) {
    RebinNeeded = YES;
    fprintf (stderr, " NRAD1D!=OldNRAD1D : %d!=%d.\n",NRAD1D,OldNRAD1D);
  }
  raderr=0;
  for (i = 0; i <= NRAD1D; i++) {
    if (fabs((Radii1D[i]-OldRadii1D[i])/Radii1D[i]) > 1e-9) {
      RebinNeeded = YES;
      raderr=1;
    }
  }
  if (raderr==1) fprintf (stderr, " Radii error.\n");
  if (!RebinNeeded) return;
  printf ("Restart/Old 1Dmesh mismatch. Rebin needed.\n");
  OldArray = (real *)malloc(OldNRAD1D*sizeof(real));
  NewArray = (real *)malloc(NRAD1D*sizeof(real));
  if ((OldArray == NULL) || (NewArray == NULL)) {
    mastererr ("Not enough memory left for rebining.\n");
    mastererr ("Aborted.\n");
    prs_exit (1);
  }
  for (type = 0; type < 4; type++) {
    Old_r = OldRmed1D;
    memcpy (New_r, Rmed1D, (NRAD1D+1)*sizeof(double));
    dangle = 0.0;
    switch (type) {
    case 0: sprintf (radix, "dens1D"); break;
    case 1: sprintf (radix, "vrad1D"); 
      Old_r = OldRadii1D; 
      memcpy (New_r, Radii1D, (NRAD1D+1)*sizeof(double));
      break;
    case 2: sprintf (radix, "vtheta1D"); dangle = 0.5; break;
    case 3: sprintf (radix, "label1D"); break;
    }
    for (i = 0; i < NRAD1D; i++) {
      if (New_r[i] < Old_r[0]) New_r[i] = Old_r[0];
      if (New_r[i] > Old_r[OldNRAD1D-1]) New_r[i] = Old_r[OldNRAD1D-1];
    }
    sprintf (filename, "%sgas%s%d.dat", OUTPUTDIR, radix, nb);
    ARR = fopen (filename, "r");
    if (ARR != NULL) {
      fread (OldArray, sizeof(real), OldNRAD1D, ARR);
      fclose (ARR);
      for (i = 0; i < NRAD1D; i++) {
	r = New_r[i];
	iold = 0;
	found = NO;
	while ((iold < OldNRAD1D) && (!found)) {
	  if (Old_r[iold+1] <= r) iold++;
	  else found = YES;
	}
	if (r <= Old_r[0]) {
	  iold = 0;
	  ifrac = 0.0;
	} else if (r >= Old_r[OldNRAD1D-1]) {
	  iold = OldNRAD1D-2;
	  ifrac = 1.0;
	} else
	  ifrac = (r-Old_r[iold])/(Old_r[iold+1]-Old_r[iold]);

	NewArray[i] = OldArray[iold]*(1.0-ifrac)+OldArray[iold+1]*ifrac;  // 15/12
      }
      ARR = fopen (filename, "w");
      if (ARR != NULL) {
	fwrite (NewArray, sizeof(real), NRAD1D, ARR);
	fclose (ARR);
      }
    }
    else {
      mastererr("Could not rebin %s. File not found\n", filename);
    }
  }
}
