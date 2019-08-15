/***********************************************/
/*                                             */
/*                                             */
/*  Fake functions library for non fftw built  */
/*                                             */
/*                                             */
/***********************************************/

#define fftw_real 2
#define fftw_complex 3

typedef int rfftwnd_mpi_plan;

void fftw_malloc();
void rfftw2d_mpi_create_plan();
void rfftwnd_mpi();
