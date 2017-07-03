#ifndef VARIABLES_H_   /* Include guard */
#define VARIABLES_H_

/// Constants ///
#define bigG 1.0
#define r0 1.0

/// Density ///
double density2D(double R);
double density3D(double R, double z);

/// Pressure + Temperature ///
double flaringIndex();
double scaleHeight(double R);
double aspectRatio(double R);
double soundSpeed(double R);
double pressure(double R, double z);

/// Azimuthal Velocity ///
double omegaK(double R);
double rotating_omegaK();
double omega3D(double R, double z);
double omegaPower(double R, double z);
double azimuthalVelocity2D(double R);
double azimuthalVelocity3D(double R, double z);

/// Radial Velocity ///
double cylindricalRadialVelocity(double R, double z);
double radialVelocity_rComponent(double R, double theta, double z);
double radialVelocity_thetaComponent(double R, double theta, double z);

/// Viscosity ///
double simpleViscosityNu(double input, double R, double z);
double zProfileViscosity(double z, double visc_lower_amplitude, double visc_upper_amplitude);
double simpleViscosityRadialOffset(double visc, double R, double z);
double combinedViscosityRadialOffset(double visc, double R, double z);
double viscosityNu(double R, double z);

/// External Torque ///
double zProfile_MagneticWind(double z);
double zProfile_HallEffect(double z);
double integratedMagneticTorqueTerm(double R, double z);
double magneticAccretionRate(double z);
double externalTorque(double R, double z);
double externalMagneticForce(double R, double z, double density);

/// Potential ///
double planetMass();
double distanceToPlanet(double r, double R, double angle);
double smoothingLength();
double stellarPotential(double r);
double planetPotential(double r, double R, double angle);

#endif // VARIABLES_H_
