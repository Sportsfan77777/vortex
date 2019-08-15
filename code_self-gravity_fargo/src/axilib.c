#include "mp.h"

void mpi_make1Dprofile (gridfield, axifield)
     real* gridfield;
     real* axifield;
{
  MPI_Request req1, req2;
  int i, j, l;
  real *localaxifield;
  localaxifield = (real*) malloc(sizeof(real) * NRAD);
  if ( localaxifield == NULL ) 
    erreur ("Not enough memory in axisglib.c ; suspect...");
  
  for ( i = 0; i < NRAD; i++ )
    localaxifield[i] = 0.;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      localaxifield[i] += gridfield[l];
    }
    localaxifield[i] /= (real)NSEC;
  }
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD; i++ )
      axifield[i] = localaxifield[i];
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for ( i = 0; i < GLOBALNRAD; i++ ) {
	if ( i < Max_or_active )
	  axifield[i] = localaxifield[i];
	else
	  axifield[i] = 0.;
      }
      MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &stat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &stat);
      for (i = Zero_or_active; i < Max_or_active; i++)
	axifield[i+IMIN] = localaxifield[i];
      if ( CPU_Rank != CPU_Highest ) {
	MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait (&req1, &stat);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  }
  free (localaxifield);
}
