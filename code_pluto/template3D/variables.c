/* ///////////////////////////////////////////////////////////////////// */
/* 

Initializes variables with functions (to be used in init.c) 


*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/// Density ///

double density2D(double R) {
  // 2-D surface density (sigma) --- Power Law
  return g_inputParam[P_Sigma0] * pow(R / r0, -g_inputParam[P_DensityPower]);
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
   return (-g_inputParam[P_TemperaturePower] + 1) / 2.0;
}

double scaleHeight(double R) {
   // Scale Height (H) --- Set by temperature profile (flaring index)
   double h, f;
   h = g_inputParam[P_AspectRatio];
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
   //return 2.0*CONST_PI/(R*sqrt(R)); // old method
   double omega_sq;
   omega_sq = bigG * g_inputParam[P_Mstar] / pow(R, 3);
   return sqrt(omega_sq);
}

double rotating_omegaK() {
   // Angular Velocity of Rotating Frame
   if (ROTATING_FRAME) return omegaK(r0);
   else return 0;
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
   
   q = g_inputParam[P_TemperaturePower];
   p = g_inputParam[P_DensityPower];

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   coeff_a = vtheta2D(R);
   term_a = q * (R / pow(R_sq + z_sq, 0.5));
   term_b = 1 - q;
   term_c = (q + p) * pow(aspectRatio(R), 2);

   return coeff_a * sqrt(term_a + term_b - term_c) - R * rotating_omegaK();
}

/// Potential ///

double planetMass() {
   // Planet's Mass in units of Jupiter-to-Sun mass ratio (M_p)
   return (g_inputParam[P_Mplanet]) * (CONST_Mjup / CONST_Msun); // in units of the Jupiter-to-Sun mass ratio
}

double distanceToPlanet(double r, double R, double angle) {
   // Radial Distance to Planet w/ Smoothing Length (d)
   double d_sq, rs_sq;

   d_sq = pow(r, 2) + pow(r0, 2) - (2 * R * r0) * cos(angle);
   rs_sq = pow(smoothingLength(), 2);

   return sqrt(d_sq + rs_sq);
}

double smoothingLength() {
   // Smoothing Length (either H, r_h, or 0.0)
   double reference_radius;

   if (g_inputParam[P_SMOOTH_SCALE_HEIGHT]) reference_radius = g_inputParam[P_AspectRatio];
   else if (g_inputParam[P_SMOOTH_HILL_RADIUS]) reference_radius = r0 * pow(planetMass() / g_inputParam[P_Mstar], 1.0 / 3.0);
   else reference_radius = 0.0;

   return 0.6 * reference_radius;
}

double stellarPotential(double r) {
   // Gravitational Potential of Star (\phi_*)
   return -bigG * g_inputParam[P_Mstar] / r;
}

double planetPotential(double r, double R, double angle) {
   // Gravitational Potential of Planet (\phi_p)
   // Note: r is spherical radius, R is cylindrical radius
   double d, planet_mass, indirect_term, phi;

   d = distanceToPlanet(r, R, angle); // including smoothing length
   planet_mass = planetMass();

   phi = -bigG * planet_mass / d;

   if (g_inputParam[P_INDIRECT_TERM]) {
       indirect_term = -planet_mass * R * cos(angle) / pow(r0, 2);
       phi += indirect_term;
   }

   return phi;
}

