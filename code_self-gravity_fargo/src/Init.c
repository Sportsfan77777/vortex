#include "mp.h"

extern boolean Restart;
extern int     NbRestart;

void ReadfromFile (array, fileprefix, filenumber)
     PolarGrid *array;
     char *fileprefix;
     int filenumber;
{
  int nr, ns, c, foo=0;
  real *field;
  char name[256];
  FILE *input;
  /* Simultaneous read access to the same file have been observed to
     give wrong results. */
  /* A sequential reading is imposed below. */
  /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Prev, 10, MPI_COMM_WORLD, &stat);
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
    return;
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); 
    /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message
     that it expects */
  if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything
				   meanwhile */
}

void InitLabel (array, sys)
     PolarGrid *array;
     PlanetarySystem *sys;
{
  int nr, ns, i, j, l;
  real xp, yp, rp;
  real x, y, angle, distance, rhill;
  real *field;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  xp = sys->x[0];
  yp = sys->y[0];
  rp = sqrt ( xp*xp + yp*yp );
  rhill = rp * pow( sys->mass[0]/3., 1./3 );
   /* Initialize label as you wish. In this example, label only takes
      into account fluid elements inside the planet's Hill Sphere */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      angle = (real)j/(real)ns*2.0*PI;
      x = Rmed[i] * cos(angle);
      y = Rmed[i] * sin(angle);
      distance = sqrt( (x - xp)*(x - xp) + (y - yp)*(y - yp) );
      if ( distance < rhill )
	field[l] = 1.0;
      else
	field[l] = 0.0;
    }
  }
}

void Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, pla_sys)
     PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_energy, *gas_label;
     PlanetarySystem *pla_sys;
{
  extern boolean Adiabatic;
  real *energ, *dens;
  int i, j, l, nr, ns;
  energ = gas_energy->Field;
  nr = gas_energy->Nrad;
  ns = gas_energy->Nsec;
  ReadPrevDim ();
  InitEuler (gas_v_rad, gas_v_theta, gas_density, gas_energy);
  InitLabel (gas_label, pla_sys);
  if (Restart == YES) {
    CheckRebin (NbRestart);
    MPI_Barrier (MPI_COMM_WORLD);
    /* Don't start reading before master has finished rebining... */
    /* It shouldn't be a problem though since a sequential read is */
    /* imposed in the ReadfromFile function below */
    mastererr ("Reading restart files...");
    fflush (stderr);
    ReadfromFile (gas_density, "gasdens", NbRestart);
    ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
    ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
    if (Adiabatic) {
      ReadfromFile (gas_energy, "gasTemperature", NbRestart);
      /* ! gas_energy accounts for the gas temperature... */
      dens = gas_density->Field;
      for (i=0; i<nr; i++) {
	for (j=0; j<ns; j++) {
	  l = i*ns + j;
	  energ[l] = dens[l]*energ[l]/(ADIABATICINDEX-1.0);
	  /* this is e = dens*temp / (gamma-1) */
	}
      }
    }
    ReadfromFile (gas_label, "gaslabel", NbRestart);
    if (StoreSigma) RefillSigma (gas_density);
    /* StoreEnergy = NO if Adiabatic = NO */
    if (StoreEnergy) RefillEnergy (gas_energy);
    fprintf (stderr, "done\n");
    fflush (stderr);
  }
  WriteDim (); 
}
