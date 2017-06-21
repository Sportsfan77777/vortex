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
      R = x1 * sin(x2)
      z = x1 * cos(x2)
  #endif

  *nu1 = kinematic_viscosity(R, z);
  *nu2 = 0.0; // bulk viscosity
}

double kinematic_viscosity(double R, double z) {
  /* viscosity profile */
  double visc_lower_amplitude, visc_upper_amplitude;
  double lower_alpha, upper_alpha, lower_accretion_rate, upper_accretion_rate;
  double ramp, ramp_amplitude, ramp_center, ramp_width, negative_z_angle, positive_z_angle;
  double viscosity;

  if (g_inputParam[P_ViscosityType] == 1) {
    // alpha viscosity
    lower_alpha = g_inputParam[P_BaseViscosity];
    upper_alpha = g_inputParam[P_MaxViscosity];

    visc_lower_amplitude = lower_alpha * soundSpeed(r0) * scaleHeight(r0);
    visc_upper_amplitude = upper_alpha * soundSpeed(r0) * scaleHeight(r0);
  }
  else if (g_inputParam[P_ViscosityType] == 2) {
    // mass accretion rate
    lower_accretion_rate = g_inputParam[P_BaseViscosity];
    upper_accretion_rate = g_inputParam[P_MaxViscosity];

    visc_lower_amplitude = lower_accretion_rate / density2D(r0);
    visc_upper_amplitude = upper_accretion_rate / density2D(r0);
  }
  else if (g_inputParam[P_ViscosityType] == 3) {
    // constant
    visc_lower_amplitude = g_inputParam[P_BaseViscosity];
    visc_upper_amplitude = g_inputParam[P_MaxViscosity];
  }
  else {
    // zero
    return 0.0; // exit!
  }

  // The rest of the amplitude

  // Z-profile
  ramp_amplitude = visc_upper_amplitude - visc_lower_amplitude;

  positive_z_angle = (z - ramp_center) / ramp_width;
  negative_z_angle = (z + ramp_center) / ramp_width;
  ramp = 0.5 * (2.0 + tanh(positive_z_angle) - tanh(negative_z_angle))

  z_profile = 1.0 + (ramp_amplitude - 1.0) * ramp;
  
  viscosity = visc_amplitude * z_profile;
  
  return viscosity;
}