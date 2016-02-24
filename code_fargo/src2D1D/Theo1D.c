#include "mp.h"

void FillSigma1D() {
  int i;
  for (i = 0; i < NRAD1D; i++) {
    SigmaMed1D[i] = Sigma(Rmed1D[i]);
    SigmaInf1D[i] = Sigma(Rinf1D[i]);
  }
}

void RefillSigma1D (Surfdens)
     PolarGrid1D *Surfdens;
{
  int i, nr;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = field[i];
    SigmaMed1D[i] = moy;
  }
  SigmaInf1D[0] = SigmaMed1D[0];
  for (i = 1; i < nr; i++) {
    SigmaInf1D[i] = (SigmaMed1D[i-1]*(Rmed1D[i]-Rinf1D[i])+\
		    SigmaMed1D[i]*(Rinf1D[i]-Rmed1D[i-1]))/\
                    (Rmed1D[i]-Rmed1D[i-1]);
  }
}
