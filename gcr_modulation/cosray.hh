#include <cmath>
#include <cstring>
#include <string>

// Physical scales
#define unit_length_fluid 1.496e+13
#define unit_velocity_fluid 1.000e+07
#define unit_number_density_fluid 1.0
#define unit_energy_particle 1.602e-12
#define unit_charge_particle 4.803e-10
#define unit_frequency_fluid (unit_velocity_fluid / unit_length_fluid)
#define unit_diffusion_fluid (unit_velocity_fluid * unit_length_fluid)
#define unit_mass_particle (unit_energy_particle / unit_velocity_fluid / unit_velocity_fluid)
#define unit_density_fluid (unit_mass_particle * unit_number_density_fluid)
#define unit_magnetic_fluid (unit_velocity_fluid * sqrt(unit_density_fluid))
#define charge_mass_particle (unit_charge_particle * unit_magnetic_fluid / unit_mass_particle / unit_velocity_fluid / unit_frequency_fluid)

// Properties of the simulation
const int nbins = 100;
const int ntraj = 8000;
const double c_code = 2.998e+10 / unit_velocity_fluid;
const double c2_code = c_code * c_code;
const double m_p = 1.673e-24 / unit_mass_particle;
const double q_p = 4.803e-10 * charge_mass_particle / unit_charge_particle;
const double B_0 = 10.0e-05 / unit_magnetic_fluid;
const double u_0 = 9.5e+07 / unit_velocity_fluid;
const double omega = 2.0 * M_PI / (25.0 * 24.0 * 3600.0) / unit_frequency_fluid;
const double kappa_0 = 1.5e22 / unit_diffusion_fluid;
const double eta = 0.5;
const double f_0 = 1.0;
const double J_0 = 100;
const double AU_cgs = 1.496e+13 / unit_length_fluid;
const double r_in = 1.0 * AU_cgs;
const double lambda = 0.05 * AU_cgs;
const double r_0 = 1.0 * AU_cgs;
const double r_min = 0.05 * AU_cgs;
const double r_max = 80.0 * AU_cgs;
const double MeV_cgs = 1.602e-06 / unit_energy_particle;
const double T_b = 150.0 * MeV_cgs;
const double T_0 = 1000.0 * MeV_cgs;
const double T_min = 90.0 * MeV_cgs;
const double T_max = 50000.0 * MeV_cgs;
const double cfl = 0.5;
const double delta = 1.0e-04 * AU_cgs;
const double lc_0 = 0.001 * AU_cgs;
const double dB2_0 = 5.0e-11 / unit_magnetic_fluid / unit_magnetic_fluid;

const std::string trajectory_fname = "trajectory.dat";
const std::string intensity_fname = "intensity_1au.dat";

// Square function
inline double Sqr(double x) {return x * x;};

// Cube function
inline double Cube(double x) {return x * x * x;};

// Quintic function
inline double Quint(double x) {return x * x * x * x * x;};

// Kinetic energy from momentum
inline double EnrKin(double mom) {return c_code * sqrt(Sqr(mom) + Sqr(m_p * c_code)) - m_p * c2_code;};

// Momentum from kinetic energy
inline double Mom(double T) {return sqrt(T * (T + 2.0 * m_p * c2_code)) / c_code;};

const double p_0 = Mom(T_0);

// Velocity from momentum
inline double Vel(double mom) {return mom / (m_p * sqrt(1.0 + Sqr(mom / (m_p * c_code))));};

// Relativistic (Lorentz) factor
inline double RelFactor(double vel) {return 1.0 / sqrt(1.0 - Sqr(vel / c_code));};

// Cyclotron frequency
inline double CyclotronFreq(double vel, double Bmag) {return q_p * Bmag / (RelFactor(vel) * m_p * c_code);};

const double lnp_min = log(Mom(T_min));
const double lnp_max = log(Mom(T_max));
const double dlnp = (lnp_max - lnp_min) / nbins;

// Norm function
inline double Norm(double* v) {return sqrt(Sqr(v[0]) + Sqr(v[1]) + Sqr(v[2]));};

// Cross product a x b = c
inline void Cross(double* a, double* b, double* c)
{
   c[0] = a[1]*b[2]-a[2]*b[1];
   c[1] = a[2]*b[0]-a[0]*b[2];
   c[2] = a[0]*b[1]-a[1]*b[0];
};

// Unit vector function
inline void UnitVec(double* v_in, double* v_out)
{
   double mag = Norm(v_in);
   v_out[0] = v_in[0] / mag;
   v_out[1] = v_in[1] / mag;
   v_out[2] = v_in[2] / mag;
};

// Get perpendicular vector function
inline void GetPerpVec(double* v_in, double* v_out)
{
   int max_ang;
   double unit_vec[3] = {0.0};

// Find the unit vector making the largest angle with "first" (smallest component of "first")
   if(fabs(v_in[0]) < fabs(v_in[1])) {
      if(fabs(v_in[0]) < fabs(v_in[2])) max_ang = 0;
      else max_ang = 2;
   }
   else {
      if(fabs(v_in[1]) < fabs(v_in[2])) max_ang = 1;
      else max_ang = 2;
   };

   unit_vec[max_ang] = 1.0;
   Cross(v_in, unit_vec, v_out);
};

// Change vector from basis to Cartesian
inline void ChangeFromBasis(double basis[][3], double* v)
{
   double v_new[3];
   v_new[0] = v[0] * basis[0][0] + v[1] * basis[1][0] + v[2] * basis[2][0];
   v_new[1] = v[0] * basis[0][1] + v[1] * basis[1][1] + v[2] * basis[2][1];
   v_new[2] = v[0] * basis[0][2] + v[1] * basis[1][2] + v[2] * basis[2][2];
   std::memcpy(v, v_new, 3 * sizeof(double));
};