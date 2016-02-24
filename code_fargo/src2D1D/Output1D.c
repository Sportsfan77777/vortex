#include "mp.h"

//extern boolean  Write_Density, Write_Velocity, IsDisk;

static real writebuf[MAX1D];


void WriteDiskPolar1D(array, number)
PolarGrid1D 	*array;
int 	         number;
{
  int             Nr, Increment;
  FILE           *dump;
  char 		name[80];
  MPI_Request req;
  real 		*ptr;
  ptr = array->Field;
  if (!(IAmTheFirst || IAmTheLast)) return;
  if (CPU_Master)
    sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  Nr = array->Nrad;
  memcpy (writebuf, ptr, sizeof(real)*Nr);
  Increment = IINNER+CPUOVERLAP;
  if (RMIN1D > RMIN) Increment = 0;
  if (IAmTheLast) {
    MPI_Isend (writebuf+Increment,NRAD1D-Increment,MPI_DOUBLE,\
	       0,24,MPI_COMM_WORLD,&req);
    MPI_Wait (&req, &stat);
  }
  if (IAmTheFirst) {
    MPI_Irecv (writebuf+Increment,NRAD1D-Increment,MPI_DOUBLE,\
	       CPU_Number-1,24,MPI_COMM_WORLD,&req);
    MPI_Wait (&req, &stat);
  }
  if (!(IAmTheFirst)) return;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s%d.dat'...", array->Name, number);
  fflush (stdout);
  fwrite (writebuf, sizeof(real), Nr,dump);
  fclose(dump);
  fprintf(stdout, "%d/", CPU_Rank);  
  fflush(stdout);
  masterprint("done\n");
  fflush(stdout);
}

void WriteDim1D () {	  
  char filename[200];
  FILE 	*dim;
  if (!CPU_Master) return;
  sprintf (filename, "%sdims1D.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    prs_exit (1);
  }
  fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\n",0,0,0,0,RMAX1D, NTOT/NINTERM, NRAD1D);
  fclose (dim);
}


void Write1DSummary(array, array1D, number)
     PolarGrid 	*array;
     PolarGrid1D 	*array1D;
     int 	         number;
{
  int             i,i0,i1,j,Ns;
  FILE           *dump;
  char 		  name[80];
  char            command[1024];
  real            average;
  real 		 *field, *field1D;
  field = array->Field;
  field1D = array1D->Field;
  if (CPU_Master)
    sprintf (name, "%s%s.ascii_rad.%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s.ascii_rad.%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Ns = array->Nsec;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf (stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s.ascii_rad.%d.dat'...", array->Name, number);
  fflush (stdout);
  fprintf(stdout, "%d/", CPU_Rank);  
  if ((RMIN1D<RMIN) && IAmTheFirst) {
    for (i = 0; i < IINNER+CPUOVERLAP; i++) {
      fprintf(dump, "%.8g %.18g\n",Rmed1D[i],field1D[i]);
    }
  }
  i0=CPUOVERLAP;
  i1=NRAD-CPUOVERLAP;
  if ((IAmTheLast) && (RMAX1D < RMAX)) i1 = NRAD;
  if ((IAmTheFirst) && (RMIN1D > RMIN)) i0 = 0;
  for (i = i0; i < i1; i++) {
    average=0.;
    for (j=0; j<Ns; j++) average+=field[j+i*Ns];
    average/=Ns;
    fprintf(dump, "%.8g %.18g\n",Rmed[i],average);
  }
  if ((RMAX1D>RMAX) && (IAmTheLast)) {
    for (i = IINNER+GLOBALNRAD-CPUOVERLAP; i < NRAD1D; i++)
      fprintf(dump, "%.8g %.18g\n",Rmed1D[i],field1D[i]);
  }
  fclose (dump);
  MPI_Barrier (MPI_COMM_WORLD);
  if (Merge) {
    if (!CPU_Master) return;
    message ("Merging output files...");
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; cat %s.ascii_rad.%d.dat.%05d >> %s.ascii_rad.%d.dat",\
	       OUTPUTDIR, array->Name, number, i, array->Name, number);
      system (command);
    }
    sprintf (command, "cd %s; rm -f %s.ascii_rad.%d.dat.0*",\
	     OUTPUTDIR, array->Name, number);
    system (command);
  }
  masterprint("done\n");
}

void Write1DWeightedSummary(array, weight, array1D, number)
PolarGrid 	*array, *weight;
PolarGrid1D 	*array1D;
int 	         number;
{
  int             i,i0,i1,j,l,Ns;
  FILE           *dump;
  char 		  name[80];
  char            command[1024];
  real            sum, w_sum;
  real 		 *field, *field1D, *w;
  field = array->Field;
  field1D = array1D->Field;
  w = weight->Field;
  if (CPU_Master)
    sprintf (name, "%s%s.ascii_rad.%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s.ascii_rad.%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Ns = array->Nsec;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s.ascii_rad.%d.dat'...", array->Name, number);
  fflush (stdout);
  fprintf(stdout, "%d/", CPU_Rank);  
  if ((RMIN1D<RMIN) && IAmTheFirst) {
    for (i = 0; i < IINNER+CPUOVERLAP; i++)
      fprintf(dump, "%g %g\n",Rmed1D[i],field1D[i]);
  }
  i0=CPUOVERLAP;
  i1=NRAD-CPUOVERLAP;
  if ((IAmTheLast) && (RMAX1D < RMAX)) i1 = NRAD;
  if ((IAmTheFirst) && (RMIN1D > RMIN)) i0 = 0;
  for (i = i0; i < i1; i++) {
    sum=0.;
    w_sum=0.;
    for (j=0; j<Ns; j++) {
      l = j+i*Ns;
      sum+=field[l]*w[l];
      w_sum+=w[l];
    }
    if(w_sum==0.) sum=0.;
    else sum/=w_sum;
    fprintf(dump, "%g %g\n",Rmed[i],sum);
  }
  if ((RMAX1D>RMAX) && (IAmTheLast)) {
    for (i = IINNER+GLOBALNRAD-CPUOVERLAP; i < NRAD1D; i++)
      fprintf(dump, "%g %g\n",Rmed1D[i],field1D[i]);
  }
  fclose (dump);
  MPI_Barrier (MPI_COMM_WORLD);
  if (Merge) {
    if (!CPU_Master) return;
    message ("Merging output files...");
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; cat %s.ascii_rad.%d.dat.%05d >> %s.ascii_rad.%d.dat",\
	       OUTPUTDIR, array->Name, number, i, array->Name, number);
      system (command);
    }
    sprintf (command, "cd %s; rm -f %s.ascii_rad.%d.dat.0*",\
	     OUTPUTDIR, array->Name, number);
    system (command);
  }
  masterprint("done\n");
}
