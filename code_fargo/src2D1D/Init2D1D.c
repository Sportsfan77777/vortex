#include "mp.h"

extern boolean Restart;
extern int     NbRestart;

void ReadfromFile (array, fileprefix, filenumber)
PolarGrid *array;
char *fileprefix;
int filenumber;
{
  int nr,ns,c, foo=0;
  real *field;
  char name[256];
  FILE *input;
/* Simultaneous read access to the same file have been observed to give wrong results. */
/* A sequential reading is imposed below. */
		      /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &stat);

  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    return;
  }
  else {
    fprintf (stdout, "> reading %s\n", name);
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message that it expects */
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything meanwhile */
}


void ReadfromFile1D (array, fileprefix, filenumber)
PolarGrid1D *array;
char *fileprefix;
int filenumber;
{
  int nr, foo=0;
  real *field;
  char name[256];
  FILE *input;
/* Simultaneous read access to the same file have been observed to give wrong results. */
/* A sequential reading is imposed below. */
		      /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &stat);

  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    return;
  }
  else {
    fprintf (stdout, "> reading %s\n", name);
  }
  field = array->Field;
  nr = array->Nrad;
  fread (field, sizeof(real), nr, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message that it expects */
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything meanwhile */
}



void InitLabel (array)
PolarGrid *array;
{
  int nr,ns,i,j,l;
  real *field;
  real rm,RM,Dr;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  rm=min2(RMIN,RMIN1D);
  RM=max2(RMAX,RMAX1D);
  Dr=RM-rm;
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      field[l] = (Rmed[i]-rm)/Dr;
    }
  }
}

void InitLabel1D (array)
PolarGrid1D *array;
{
  int nr,i;
  real *field;
  real rm,RM,Dr;
  field = array->Field;
  nr = array->Nrad;
  rm=min2(RMIN,RMIN1D);
  RM=max2(RMAX,RMAX1D);
  Dr=RM-rm;
  for (i = 0; i <= nr; i++) {
    field[i] = (Rmed1D[i]-rm)/Dr;
  }
}



void Initialization2D1D (gas_density, gas_v_rad, gas_v_theta, gas_label, gas_density_1D, gas_v_rad_1D, gas_v_theta_1D, gas_label_1D)
PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_label;
PolarGrid1D *gas_density_1D, *gas_v_rad_1D, *gas_v_theta_1D, *gas_label_1D;
{
  //  ReadPrevDim ();    Has to be done before to erase used_rad.dat.
  //  ReadPrevDim1D ();  Moved at the begining of FillPolar1DArrays2D1D.
  InitEuler2D1D (gas_density, gas_v_rad, gas_v_theta, gas_density_1D, gas_v_rad_1D, gas_v_theta_1D);
  InitLabel (gas_label);
  InitLabel1D (gas_label_1D);
  InitCompute1DGhostsVr ();
  if (Restart == YES) {
    CheckRebin (NbRestart);
    CheckRebin1D (NbRestart);
    MPI_Barrier (MPI_COMM_WORLD); /* Don't start reading before master has finished rebining... */
				  /* It shouldn't be a problem though since a sequential read is */
                                  /* imposed in the ReadfromFile function below */
    mastererr ("Reading restart files...");
    fflush (stderr);
    ReadfromFile (gas_density, "gasdens", NbRestart);
    ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
    ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
    ReadfromFile (gas_label, "gaslabel", NbRestart);
    ReadfromFile1D (gas_density_1D, "gasdens1D", NbRestart);
    ReadfromFile1D (gas_v_rad_1D, "gasvrad1D", NbRestart);
    ReadfromFile1D (gas_v_theta_1D, "gasvtheta1D", NbRestart);
    ReadfromFile1D (gas_label_1D, "gaslabel1D", NbRestart);
    if (StoreSigma) {
      RefillSigma (gas_density);
      RefillSigma1D (gas_density_1D);
    }
    mastererr (" done\n");
    fflush (stderr);
  }
  WriteDim (); 
  WriteDim1D (); 
}
