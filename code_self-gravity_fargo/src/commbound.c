#include "mp.h"

static real *SendInnerBoundary;
static real *SendOuterBoundary;
static real *RecvInnerBoundary;
static real *RecvOuterBoundary;
static int allocated_com = 0;
static int size_com;
extern boolean Adiabatic;

void AllocateComm () {
  size_com = 3;
  if (Adiabatic) size_com += 1;
  if (AdvecteLabel == YES) size_com += 1;
  size_com *= NSEC * CPUOVERLAP;
  SendInnerBoundary = malloc (size_com * sizeof(real));
  SendOuterBoundary = malloc (size_com * sizeof(real));
  RecvInnerBoundary = malloc (size_com * sizeof(real));
  RecvOuterBoundary = malloc (size_com * sizeof(real));
  if ((SendInnerBoundary == NULL) ||\
      (SendOuterBoundary == NULL) ||\
      (RecvInnerBoundary == NULL) ||\
      (RecvOuterBoundary == NULL)) {
    fprintf (stderr, "CPU %d had not enough memory to allocate communicators.\n", CPU_Rank);
    prs_exit(0);
  }
  allocated_com = 1;
}

void CommunicateBoundaries (Density, Vrad, Vtheta, Energy, Label)
     PolarGrid *Density, *Vrad, *Vtheta, *Energy, *Label;
{
  MPI_Request req1, req2, req3, req4;
  int l, oo, o, nr;
  if (!allocated_com) AllocateComm ();
  l = CPUOVERLAP*NSEC;
  nr = Density->Nrad;
  oo = (nr-CPUOVERLAP)*NSEC;
  o = (nr-2*CPUOVERLAP)*NSEC;
  memcpy (SendInnerBoundary, Density->Field+l, l*sizeof(real));
  memcpy (SendInnerBoundary+l, Vrad->Field+l, l*sizeof(real));
  memcpy (SendInnerBoundary+2*l, Vtheta->Field+l, l*sizeof(real));
  memcpy (SendOuterBoundary, Density->Field+o, l*sizeof(real));
  memcpy (SendOuterBoundary+l, Vrad->Field+o, l*sizeof(real));
  memcpy (SendOuterBoundary+2*l, Vtheta->Field+o, l*sizeof(real));
  if (Adiabatic) {
    memcpy (SendInnerBoundary+3*l, Energy->Field+l, l*sizeof(real));
    memcpy (SendOuterBoundary+3*l, Energy->Field+o, l*sizeof(real));
  }
  if (AdvecteLabel == YES) {
    memcpy (SendInnerBoundary+(3+(Adiabatic == YES ? 1:0))*l, Label->Field+l, l*sizeof(real));
    memcpy (SendOuterBoundary+(3+(Adiabatic == YES ? 1:0))*l, Label->Field+o, l*sizeof(real));
  }
  /* ------------------------------------------ */
  /* Note that boundary exchange is independant */
  /* from chosen domain decomposition           */
  /* ------------------------------------------ */
  if (CPU_Rank%2 == 0) {
    if (CPU_Rank > 0) {
      MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
    }
    if (CPU_Rank != CPU_Highest) {
      MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
      MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
    }
  } else {
    if (CPU_Rank != CPU_Highest) {
      MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
      MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
    }
    if (CPU_Rank > 0) {
      MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
    }
  }
  if (CPU_Rank > 0) {
    MPI_Wait (&req1, &stat);
    MPI_Wait (&req2, &stat);
    memcpy (Density->Field, RecvInnerBoundary, l*sizeof(real));
    memcpy (Vrad->Field, RecvInnerBoundary+l, l*sizeof(real));
    memcpy (Vtheta->Field, RecvInnerBoundary+2*l, l*sizeof(real));
    if (Adiabatic)
      memcpy (Energy->Field, RecvInnerBoundary+3*l, l*sizeof(real));
    if (AdvecteLabel == YES)
      memcpy (Label->Field, RecvInnerBoundary+(3+(Adiabatic == YES ? 1:0))*l, l*sizeof(real));
  }
  if (CPU_Rank != CPU_Highest) {
    MPI_Wait (&req3, &stat);
    MPI_Wait (&req4, &stat);
    memcpy (Density->Field+oo, RecvOuterBoundary, l*sizeof(real));
    memcpy (Vrad->Field+oo, RecvOuterBoundary+l, l*sizeof(real));
    memcpy (Vtheta->Field+oo, RecvOuterBoundary+2*l, l*sizeof(real));
    if (Adiabatic)
      memcpy (Energy->Field+oo, RecvOuterBoundary+3*l, l*sizeof(real));
    if (AdvecteLabel == YES)
      memcpy (Label->Field+oo, RecvOuterBoundary+(3+(Adiabatic == YES ? 1:0))*l, l*sizeof(real));
  }
}
