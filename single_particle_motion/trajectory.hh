#include <cmath>

// Particle and equation constants
const double q = 1.0;         // charge
const double m = 1.0;         // mass
const double c = 1.0;         // speed of light
const double mc = m * c;      // mass times speed of light
const double q_c = q / c;     // charge divided by speed of light
const double s = 2.0;         // mystery constant

// Magnetic field
const double Bx = 0.0;        // magnetic field in the x direction
const double By = 0.0;        // magnetic field in the y direction
const double Bz = 1.0;        // magnetic field in the z direction

// Electric field
const double Ex = 0.0;        // electric field in the x direction
const double Ey = 0.0;        // electric field in the y direction
const double Ez = 0.0;        // electric field in the z direction

// Position
const double Rx = -10.0;      // initial x position
const double Ry = 0.0;        // initial y position
const double Rz = 0.0;        // initial z position

// Momentum
const double Px = 0.0;        // initial x momentum
const double Py = 1.0;        // initial y momentum
const double Pz = 0.0;        // initial z momentum

// Numerical integration parameters
const double dmax = 5;      // maximum spatial displacement per step
const double tmax = 1000;    // maximum time to integrate
const double CFL = 0.5;       // CFL condition

// Square function
inline double Sqr(double x)
{
   return (x * x);
};

// Norm function
inline double Norm(double* v)
{
   return sqrt(Sqr(v[0])+Sqr(v[1])+Sqr(v[2]));
};

// Cross product a x b = c
inline void Cross(double* a, double* b, double* c)
{
   c[0] = a[1]*b[2]-a[2]*b[1];
   c[1] = a[2]*b[0]-a[0]*b[2];
   c[2] = a[0]*b[1]-a[1]*b[0];
};

// Function to convert momentum to velocity in magnitude
inline double Vel(double mom)
{
   return mom / (m * sqrt(1.0 + Sqr(mom / mc)));
};

// Function to convert momentum to velocity in vector form
inline void Vel(double* mom, double* vel)
{
   double mmag = Norm(mom);
   double ratio = Vel(mmag) / mmag;
 
   vel[0] = mom[0] * ratio;
   vel[1] = mom[1] * ratio;
   vel[2] = mom[2] * ratio;
};

// Function to compute gyroradius
inline double LarmorRadius(double mmag, double Bmag)
{
   return mmag * c / (q * Bmag);
};
