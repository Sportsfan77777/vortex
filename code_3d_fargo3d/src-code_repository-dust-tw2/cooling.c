//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>
//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>
void Edamp_cpu(real dt) {
//<USER_DEFINED>
  INPUT(Energy);
  OUTPUT(Energy);
  //<\USER_DEFINED>
  //<INTERNAL>
  int i;
  int j;
  int k;
  //<\INTERNAL>

//<EXTERNAL>
real* e = Energy->field_cpu;
real edamp = EDAMP;
int pitch = Pitch_cpu;
int stride = Stride_cpu;
int size_x = Nx;
int size_y = Ny+2*NGHY;
int size_z = Nz+2*NGHZ;
//<\EXTERNAL>
//<MAIN_LOOP>
for (k=0; k<size_z; k++) {
for (j=0; j<size_y; j++) {
for (i=0; i<size_x; i++) {
//<#>
e[l] *= 1.0/(1.0+edamp*dt);
//<\#>
}
}
}
//<\MAIN_LOOP>
}