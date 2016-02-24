#include "mp.h"
#include <stdarg.h>


real GetGlobalIFrac (r)
     real r;
{
  int i=0;
  real ifrac;
  if (r < GlobalRmed[0]) return 0.0;
  if (r > GlobalRmed[GLOBALNRAD-1]) return (real)GLOBALNRAD-1.0;
  while (GlobalRmed[i] <= r) i++;
  ifrac = (real)i+(r-GlobalRmed[i-1])/(GlobalRmed[i]-GlobalRmed[i-1])-1.0;
  return ifrac;
}

void prs_exit (numb)
     int numb;
{
  MPI_Finalize ();
  exit (numb);
}

void masterprint (const char *template, ...)
{
  va_list ap;
  if (!CPU_Master) return;
  va_start (ap, template);
  vfprintf (stdout, template, ap);
  va_end (ap);
}

void mastererr (const char *template, ...)
{
  va_list ap;
  if (!CPU_Master) return;
  va_start (ap, template);
  vfprintf (stderr, template, ap);
  va_end (ap);
}

void
erreur(string)
char           *string;
{
  fprintf(stderr, "%s\n", string);
  prs_exit(1);
}

void message (msg) 
char *msg;
{
	fprintf (stdout, "%s", msg);
}

PolarGrid    *
CreatePolarGrid(Nr, Ns, name)
int             Nr, Ns;
char           *name;
{
  PolarGrid           *array;
  real           *field;
  char           *string;
  int             i, j, l;

  array = (PolarGrid *) malloc(sizeof(PolarGrid));
  if (array == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  field = (real *) malloc(sizeof(real) * (Nr + 3) * (Ns + 1) + 5);
  if (field == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  sprintf(string, "gas%s", name);
  array->Field = field;
  array->Name = string;
  array->Nrad = Nr;
  array->Nsec = Ns;
  for (i = 0; i <= Nr; i++) {
    for (j = 0; j < Ns; j++) {
      l = j + i*Ns;
      field[l] = 0.;
    }
  }
  return array;
}


PolarGrid1D    *
CreatePolarGrid1D(Nr, name)
int             Nr;
char           *name;
{
  PolarGrid1D    *array;
  real           *field;
  char           *string;
  int             i;

  array = (PolarGrid1D *) malloc(sizeof(PolarGrid1D));
  if (array == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  field = (real *) malloc(sizeof(real) * (Nr + 3) * 2 + 5);
  if (field == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL)
    erreur("Insufficient memory for PolarGrid creation");
  sprintf(string, "gas%s", name);
  array->Field = field;
  array->Name = string;
  array->Nrad = Nr;
   for (i = 0; i <= Nr; i++) {
    field[i] = 0.;
  }
  return array;
}

void MultiplyPolarGridbyConstant (arraysrc, constant)
PolarGrid *arraysrc;
real constant;
{
  int i, nr, ns;
  real *fieldsrc;

  nr = arraysrc->Nrad;
  ns = arraysrc->Nsec;

  fieldsrc  =  arraysrc->Field;
#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) {
    fieldsrc[i] *= constant;
  }
}

void MultiplyPolarGridbyConstant1D (arraysrc, constant)
PolarGrid1D *arraysrc;
real constant;
{
  int i, nr;
  real *fieldsrc;

  nr = arraysrc->Nrad;

  fieldsrc  =  arraysrc->Field;
#pragma omp parallel for
  for (i = 0; i < (nr+1); i++) {
    fieldsrc[i] *= constant;
  }
}


void DumpSources (argc, argv)
     int argc;
     char *argv[];
{
  char CommandLine[1024];
  char filecom[1024];
  int i;
  FILE *COM;
  if (!CPU_Master) return;
  sprintf (CommandLine, "cp .source.tar.bz2 %ssrc.tar.bz2", OUTPUTDIR);
  system (CommandLine);
  sprintf (filecom, "%srun.commandline", OUTPUTDIR);
  COM = fopen (filecom, "w");
  if (COM == NULL) {
    mastererr ("Could not open %s\nAborted.\n", filecom);
    prs_exit(1);
  }
  for (i = 0; i < argc; i++) {
    fprintf (COM, "%s ",argv[i]);
  }
  fclose (COM);
}
