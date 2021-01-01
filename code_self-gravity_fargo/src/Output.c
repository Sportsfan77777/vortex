#include "mp.h"

static real     Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, AccretedMass;
extern real     LostMass, OmegaFrame;
extern boolean  Write_Density, Write_Velocity, Write_Energy, IsDisk;
extern boolean  Write_Temperature, Write_DivV, Write_Qplus;

void EmptyPlanetSystemFile (sys)
     PlanetarySystem *sys;
{
  FILE *output;
  char name[256];
  int i, n;
  n = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < n; i++) {
    sprintf (name, "%splanet%d.dat", OUTPUTDIR, i);
    output = fopen (name, "w");
    if (output == NULL) {
      fprintf (stderr, "Can't write %s file. Aborting.\n", name);
      prs_exit (1);
    }
    fclose (output);
  }
}

void WritePlanetFile (TimeStep, n)
     int TimeStep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  printf ("Updating 'planet%d.dat'...", n);
  fflush (stdout);
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'planet%d.dat' file. Aborting.\n", n);
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", \
    TimeStep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, AccretedMass, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WritePlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i]*MassTaper;
    AccretedMass = sys->accreted_mass[i];
    WritePlanetFile (t, i);
  }
}
   

void WriteBigPlanetFile (TimeStep, n)
     int TimeStep;
     int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  sprintf (name, "%sbigplanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'bigplanet.dat' file. Aborting.\n");
    prs_exit (1);
  }
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", \ 
    TimeStep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, AccretedMass, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp);
  fclose (output);
}

void WriteBigPlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i]*MassTaper;
    AccretedMass = sys->accreted_mass[i]; 
    WriteBigPlanetFile (t, i);
  }
}

real GetfromPlanetFile (TimeStep, column, n)
     int TimeStep, column, n;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'planet%d.dat' file. Aborting restart.\n",n);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'planet%d.dat'. Aborting restart.\n",n);
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != TimeStep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'planet%d.dat' file. Aborting restart.\n", TimeStep,n);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}

void RestartPlanetarySystem (timestep, sys)
     PlanetarySystem *sys;
     int timestep;
{
  int k;
  for (k = 0; k < sys->nb; k++) {
    sys->x[k] = GetfromPlanetFile (timestep, 2, k);
    sys->y[k] = GetfromPlanetFile (timestep, 3, k);
    sys->vx[k] = GetfromPlanetFile (timestep, 4, k);
    sys->vy[k] = GetfromPlanetFile (timestep, 5, k);
    sys->mass[k] = GetfromPlanetFile (timestep, 6, k);
    sys->accreted_mass[k] = GetfromPlanetFile (timestep, 7, k);
  }
}

void WriteDiskPolar(array, number)
     PolarGrid 	*array;
     int 	number;
{
  int           Nr, Ns;
  FILE          *dump;
  char 		name[80];
  real 		*ptr;
  ptr = array->Field;
  if (CPU_Master)
    sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Nr = array->Nrad;
  Ns = array->Nsec;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s%d.dat'...", array->Name, number);
  fflush (stdout);
  MPI_Barrier (MPI_COMM_WORLD);
/* We strip the first CPUOVERLAP rings if the current CPU is not the 
   innermost one */
  if (CPU_Rank > 0) {
    ptr += CPUOVERLAP*Ns;
    Nr -=CPUOVERLAP ;
  }
/* We strip the last CPUOVERLAP rings if the current CPU is not the outermost
   one, equal to CPU_Highest in all cases */
  if (CPU_Rank != CPU_Highest) {
    Nr -=CPUOVERLAP;
  }
  fwrite (ptr, sizeof(real), Nr*Ns,dump);
  fclose(dump);
  fprintf(stdout, "%d/", CPU_Rank);  
  fflush(stdout);
  MPI_Barrier (MPI_COMM_WORLD);
  masterprint("done\n");
}

void WriteDim () {	  
  char filename[200];
  FILE 	*dim;
  if (!CPU_Master) return;
  sprintf (filename, "%sdims.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    prs_exit (1);
  }
  fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n",0,0,0,0,RMAX, NTOT/NINTERM, GLOBALNRAD, NSEC);
  fclose (dim);
}

void SendOutput (index, dens, gasvr, gasvt, gasenerg, label)
     int          index;
     PolarGrid   *dens, *gasvr, *gasvt, *label, *gasenerg;
{
  if (CPU_Master)
    printf ("\n*** OUTPUT %d ***\n", index);
  if (IsDisk == YES) {
      if (Write_Density == YES) WriteDiskPolar (dens, index);
      if (Write_Velocity == YES) {
	WriteDiskPolar (gasvr, index);
	WriteDiskPolar (gasvt, index);
      }
      if (Write_Energy == YES) WriteDiskPolar (gasenerg, index);
      if (Write_Temperature == YES) WriteDiskPolar (Temperature, index);
      if (Write_DivV == YES) WriteDiskPolar (DivergenceVelocity, index);
      if (Write_Qplus == YES)  WriteDiskPolar (Qplus, index);
      if (AdvecteLabel == YES) WriteDiskPolar (label, index);
      MPI_Barrier (MPI_COMM_WORLD);
      if (Merge && (CPU_Number > 1)) merge (index);
  }
}

 
