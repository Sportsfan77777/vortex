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
double viscosityNu(double R, double z);

/// Potential ///
double planetMass();
double distanceToPlanet(double r, double R, double angle);
double smoothingLength();
double stellarPotential(double r);
double planetPotential(double r, double R, double angle);

/// External Torque ///
double externalTorque(double R, double z);
double externalTorque_WindOnly(double R, double z);
double externalTorque_WindWithViscosity(double R, double z);

#endif // VARIABLES_H_
