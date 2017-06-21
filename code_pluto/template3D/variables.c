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
   return (2.0*CONST_PI) * sqrt(omega_sq);
}

double rotating_omegaK() {
   // Angular Velocity of Rotating Frame
   if (ROTATING_FRAME) return omegaK(r0);
   else return 0;
}

double omega3D(double R, double z) {
   // Angular Velocity (Omega) --- does not include rotating frame
   double q, p;
   double R_sq, z_sq;
   double coeff_a, term_a, term_b, term_c;
   
   q = g_inputParam[P_TemperaturePower];
   p = g_inputParam[P_DensityPower];

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   coeff_a = omegaK(R)
   term_a = q * (R / pow(R_sq + z_sq, 0.5));
   term_b = 1 - q;
   term_c = (q + p) * pow(aspectRatio(R), 2);

   return coeff_a * (term_a + term_b - term_c);
}

double omegaPower(double R, double z) {
   // radial exponential dependence of omega (d ln Omega / d ln r)
   double p, q, f;
   double R_sq, z_sq;
   double zeroth_order_term;
   double coeff, second_order_term;

   zeroth_order_term = -1.5;

   f = flaringIndex(); 
   q = g_inputParam[P_TemperaturePower];
   p = g_inputParam[P_DensityPower];

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   term_a = f * R * (z_sq / pow(R_sq + z_sq, 1.5));
   term_b = -f * (p + q) * pow(R, f);
   
   coeff = omegaK(R) / omega3D(R, z);
   second_order_term = term_a + term_b;

   return zeroth_order_term + 0.5 * pow(coeff, 2) * second_order_term;
}

double vtheta2D(double R) {
   // Azimuthal Velocity (v_theta = v_Keplerian)
   return R * omegaK(R);
}

double vtheta3D(double R, double z) {
   // Azimuthal Velocity (v_theta) --- includes rotating frame
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

double radialVelocity(double R, double z) {
   // Radial Velocity (v_rad)
   return omegaPower(R, z) * viscosityNu(R, z) / r
}

/// Viscosity ///

double viscosityNu(double R, double z) {
  // "viscosity component" of kinematic viscosity (nu = mu / rho)
  // viscosity profile: See MKL 2014, Section 4.1.1
  // Parameters: Cylindrical R and z

  double visc_lower_amplitude, visc_upper_amplitude;
  double lower_alpha, upper_alpha, lower_accretion_rate, upper_accretion_rate;
  double z_coor, ramp;
  double ramp_amplitude, ramp_center, ramp_width, negative_z_angle, positive_z_angle;
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
  density_factor = density3D(r0, z) / density3D(R, z);
  visc_lower_amplitude *= (density_factor * omegaPower(R, z));
  visc_upper_amplitude *= (density_factor * omegaPower(R, z));

  // Z-profile
  z_coor = z / scaleHeight(r0); // z in number of scaleights
  ramp_amplitude = visc_upper_amplitude - visc_lower_amplitude;
  ramp_center = g_inputParam[P_ViscosityRampCenter]; // in number of scale heights
  ramp_width = g_inputParam[P_ViscosityRampWidth]; // in number of scale heights

  positive_z_angle = (z_coor - ramp_center) / ramp_width;
  negative_z_angle = (z_coor + ramp_center) / ramp_width;
  ramp = 0.5 * (2.0 + tanh(positive_z_angle) - tanh(negative_z_angle))

  z_profile = 1.0 + (ramp_amplitude - 1.0) * ramp;
  
  // Full Expression
  viscosity = visc_amplitude * z_profile;
  return viscosity;
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
   return (4.0 * CONST_PI*CONST_PI) * -bigG * g_inputParam[P_Mstar] / r;
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

