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

/// Velocity ///
double omegaK(double R);
double vtheta2D(double R);
double vtheta3D(double R, double z);

#endif // VARIABLES_H_
