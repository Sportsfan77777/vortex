/* ///////////////////////////////////////////////////////////////////// */
/* 

Initializes variables with functions (to be used in init.c) 


*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/// Density ///

double density2D(double R) {
  // 2-D surface density (sigma) --- Power Law
  return g_inputParam[Sigma0_Param] * pow(R / r0, -g_inputParam[DensityPower]);
}

double density3D(double R, double z) {
  // 3-D density (rho) --- Exponential set by Hydrostatic Equilibrium
  double R_sq, z_sq;
  double coeff, exp_coeff, term_a, term_b;

  R_sq = pow(R, 2);
  z_sq = pow(z, 2);

  coeff = density2D(R) / (pow(2.0 * CONST_PI, 0.5) * scaleHeight(R));
  exp_coeff = 1.0 / pow(aspectRatio(R), 2);
  term_a = R / pow(R_sq + z_sq, 0.5);
  term_b = 1.0;

  return coeff * exp(exp_coeff * (term_a - term_b));
}

/// Pressure + Temperature ///
double flaringIndex() {
   // Flaring Index
   return (-g_inputParam[TemperaturePower] + 1) / 2.0;
}

double scaleHeight(double R) {
   // Scale Height (H) --- Set by temperature profile (flaring index)
   double h, f;
   h = g_inputParam[AspectRatio_Param];
   f = flaringIndex();

   return (h * R) * pow((R / r0), f); 
}

double aspectRatio(double R) {
   // Aspect Ratio (H / R) --- Set by temperature profile
   return scaleHeight(R) / R;
}

double soundSpeed(double R) {
    // Sound Speed (c_s) --- Set by temperature profile (aspect ratio)
    return scaleHeight(R) * omegaK(R);
}

double pressure(double R, double z) {
    // Pressure --- Set by ideal gas law (?)
    return pow(soundSpeed(R), 2) * density3D(R, z);
}


/// Velocity ///

double omegaK(double R) {
   // Keplerian Angular Velocity --- Set by Kepler's 3rd Law
   double omega_sq;
   omega_sq = bigG * g_inputParam[Mstar] / pow(R, 3);
   return sqrt(omega_sq);
}


double vtheta2D(double R) {
   // Azimuthal Velocity (v_theta = v_Keplerian)
   return R * omegaK(R);
}

double vtheta3D(double R, double z) {
   // Azimuthal Velocity (v_theta)
   double q, p;
   double R_sq, z_sq;
   double coeff_a, term_a, term_b, term_c;
   
   q = g_inputParam[TemperaturePower];
   p = g_inputParam[DensityPower];

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   coeff_a = vtheta2D(R);
   term_a = q * (R / pow(R_sq + z_sq, 0.5));
   term_b = 1 - q;
   term_c = (q + p) * pow(aspectRatio(R), 2);

   return coeff_a * sqrt(term_a + term_b - term_c);
}
