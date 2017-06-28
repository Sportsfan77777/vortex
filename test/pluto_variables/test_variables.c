/* ///////////////////////////////////////////////////////////////////// */
/* 

Initializes variables with functions (to be used in init.c) 


*/
/* ///////////////////////////////////////////////////////////////////// */

#include "stdio.h"
#include "math.h"
#include "test_variables.h"

/// Density ///

double density2D(double R) {
  // 2-D surface density (sigma) --- Power Law
  return 5.787e-5 * pow(R / r0, -densityPower);
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
   return (-temperaturePower + 1) / 2.0;
}

double scaleHeight(double R) {
   // Scale Height (H) --- Set by temperature profile (flaring index)
   double h, f;
   h = 0.06;
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


/// Azimuthal Velocity ///

double omegaK(double R) {
   // Keplerian Angular Velocity --- Set by Kepler's 3rd Law
   //return 2.0*CONST_PI/(R*sqrt(R)); // old method
   double omega_sq;
   omega_sq = bigG * 1.0 / pow(R, 3);
   return (2.0*CONST_PI) * sqrt(omega_sq);
}

double rotating_omegaK() {
   // Angular Velocity of Rotating Frame
   if (0.0) return omegaK(r0);
   else return 0;
}

double omega3D(double R, double z) {
   // Angular Velocity (Omega) --- does not include rotating frame
   double q, p;
   double R_sq, z_sq;
   double coeff_a, term_a, term_b, term_c;
   
   q = temperaturePower;
   p = densityPower;

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   coeff_a = omegaK(R);
   term_a = q * (R / pow(R_sq + z_sq, 0.5));
   term_b = 1 - q;
   term_c = (q + p) * pow(aspectRatio(R), 2);

   return coeff_a * sqrt(term_a + term_b - term_c);
}

double omegaPower(double R, double z) {
   // radial exponential dependence of omega (d ln Omega / d ln r)
   double p, q, f;
   double R_sq, z_sq;
   double zeroth_order_term;
   double coeff, second_order_term, term_a, term_b;

   zeroth_order_term = -1.5;

   f = flaringIndex(); 
   q = temperaturePower;
   p = densityPower;

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   term_a = q * R * (z_sq / pow(R_sq + z_sq, 1.5));
   term_b = -2.0 * f * (p + q) * (aspectRatio(r0)) * pow(R, 2.0 * f);
   
   coeff = omegaK(R) / omega3D(R, z);
   second_order_term = term_a + term_b;

   return zeroth_order_term + 0.5 * pow(coeff, 2) * second_order_term;
}

double azimuthalVelocity2D(double R) {
   // Azimuthal Velocity (v_phi = v_Keplerian)
   return R * omegaK(R);
}

double azimuthalVelocity3D(double R, double z) {
   // Azimuthal Velocity (v_phi) --- includes rotating frame
   double q, p;
   double R_sq, z_sq;
   double coeff_a, term_a, term_b, term_c;
   
   q = temperaturePower;
   p = densityPower;

   R_sq = pow(R, 2);
   z_sq = pow(z, 2);

   coeff_a = azimuthalVelocity2D(R);
   term_a = q * (R / pow(R_sq + z_sq, 0.5));
   term_b = 1 - q;
   term_c = (q + p) * pow(aspectRatio(R), 2);

   return coeff_a * sqrt(term_a + term_b - term_c) - R * rotating_omegaK();
}

/// Radial Velocity ///

double cylindricalRadialVelocity(double R, double z) {
   // Radial Velocity (v_R) --- set by viscous torque (v_R = nu / R * (d ln Omega / d ln R))
   return omegaPower(R, z) * viscosityNu(R, z) / R;
}

double radialVelocity_rComponent(double R, double theta, double z) {
   // Spherical Radial Velocity (v_r) --- set by coordinate change (v_r = V_R * sin theta)
   return cylindricalRadialVelocity(R, z) * sin(theta);
}

double radialVelocity_thetaComponent(double R, double theta, double z) {
   // Spherical Radial Velocity: Small Theta Component (v_theta) --- set by coordinate change (v_r = V_R * cos theta)
   return cylindricalRadialVelocity(R, z) * cos(theta);
}

/// Viscosity ///

double viscosityNu(double R, double z) {
  // "viscosity component" of kinematic viscosity (nu = mu / rho)
  // viscosity profile: See MKL 2014, Section 4.1.1
  // Parameters: Cylindrical R and z

  double visc_lower_amplitude, visc_upper_amplitude, density_factor, omega_factor;
  double lower_alpha, upper_alpha, lower_accretion_rate, upper_accretion_rate;
  double z_coor, ramp;
  double ramp_amplitude, ramp_center, ramp_width, negative_z_angle, positive_z_angle;
  double viscosity, z_profile;

  // constant
  visc_lower_amplitude = 1e-5;
  visc_upper_amplitude = 1e-4;

  //printf("Lower Viscosity: %f\n", visc_lower_amplitude);
  //printf("Upper Viscosity: %f\n", visc_upper_amplitude);

  // The rest of the amplitude
  density_factor = density3D(r0, z) / density3D(R, z);
  omega_factor = omegaPower(r0, z) / omegaPower(R, z);
  //visc_lower_amplitude *= (density_factor * omega_factor);
  //visc_upper_amplitude *= (density_factor * omega_factor);

  if (visc_lower_amplitude >= visc_upper_amplitude) {
     // No Ramp!
     //printf("No Ramp!\n");
     return visc_lower_amplitude;
  }

  ///// With Ramp /////

  // Z-profile
  z_coor = z / scaleHeight(r0); // z in number of scaleights
  ramp_amplitude = (visc_upper_amplitude - visc_lower_amplitude) / visc_lower_amplitude;
  ramp_center = 1.0; // in number of scale heights
  ramp_width = 0.3; // in number of scale heights

  positive_z_angle = (z_coor - ramp_center) / ramp_width;
  negative_z_angle = (z_coor + ramp_center) / ramp_width;
  ramp = 0.5 * (2.0 + tanh(positive_z_angle) - tanh(negative_z_angle));

  z_profile = 1.0 + (ramp_amplitude) * ramp;
/*
  printf("Postive Angle: %f, %f\n", positive_z_angle, tanh(positive_z_angle));
  printf("Negative Angle: %f, %f\n", negative_z_angle, tanh(negative_z_angle));

  printf("Ramp: %f\n", ramp);
  printf("Ramp Amplitude: %f\n", ramp_amplitude);

  printf("Z-Profile: %f\n", z_profile); */
  
  // Full Expression
  viscosity = visc_lower_amplitude * z_profile;
  return viscosity;
}

int main() {
  // do stuff
  double r, angle, z;

  r = 1.0;
  angle = (CONST_PI / 2.0) - 0.0;
  z = cos(angle) + 0.0;

  printf("z: %f\n", z);
  printf("z / H: %f\n", z / 0.06);

  printf("Azimuthal Velocity 2D: %f\n", azimuthalVelocity2D(r, z));
  printf("Omega3D: %f\n", omega3D(r, z));
  printf("Azimuthal Velocity 3D: %f\n", azimuthalVelocity3D(r, z));
  printf("Omega Power: %f\n", omegaPower(r, z));

  //printf("Viscosity: %f\n", viscosityNu(1.0, z));
  printf("Velocity: %f\n", cylindricalRadialVelocity(r, z));

  printf("Radial Velocity: %f\n", radialVelocity_rComponent(r, angle, z));
  printf("Theta  Velocity: %f\n", radialVelocity_thetaComponent(r, angle, z));
  return 0;
}
