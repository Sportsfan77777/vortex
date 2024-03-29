/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
  #if GEOMETRY == CYLINDRICAL || POLAR
      R = x1;
  #elif GEOMETRY == SPHERICAL
      R = x1 * sin(x2);
      z = x1 * cos(x2);
  #endif

  *nu1 = viscosityNu(R, z) * v[RHO]; // nu = mu * rho
  *nu2 = 0.0; // bulk viscosity
}
