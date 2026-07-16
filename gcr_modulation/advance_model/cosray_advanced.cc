#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>

#include "cosray.hh"

#define PRINT_TRAJECTORY
#define PRINT_INTENSITY
#define PARKER_FIELD
#define CURRENT_SHEET
#define MAGNETIC_DRIFT
#define SMITH_BIEBER

// Compute the background flow and field
void GetFields(double t, double* pos, double* u, double* B)
{
   double r, s, costheta, sintheta, sinphi, cosphi;

// Polar (xy) and spherical (xyz) radii
   s = sqrt(Sqr(pos[0]) + Sqr(pos[1]));
   r = Norm(pos);

// Angles for conversion to Cartesian
   costheta = pos[2] / r;
   sintheta = sqrt(fmax(1.0 - Sqr(costheta), 0.0));
   if (s < 1.0e-12) {
      cosphi = 0.0;
      sinphi = 0.0;
   } else {
      cosphi = pos[0] / s;
      sinphi = pos[1] / s;
   };
   
// Radial flow with constant speed
   u[0] = u_0 * sintheta * cosphi;
   u[1] = u_0 * sintheta * sinphi;
   u[2] = u_0 * costheta;

#ifdef PARKER_FIELD
// Parker spiral magnetic field
   double Br = B_0 * Sqr(r_0 / r);
   double Bp = -Br * r * omega / u_0 * sintheta;
#ifdef SMITH_BIEBER
   Bp -= 0.1 * Br * r * omega / u_0;
#endif
   B[0] = Br * sintheta * cosphi - Bp * sinphi;
   B[1] = Br * sintheta * sinphi + Bp * cosphi;
   B[2] = Br * costheta;
#else
// Radial magnetic field with constant strength
   B[0] = B_0 * sintheta * cosphi;
   B[1] = B_0 * sintheta * sinphi;
   B[2] = B_0 * costheta;
#endif
#ifdef CURRENT_SHEET
// Flat equatorial current sheet
// if z < 0: positive cycle (A > 0)
// if z > 0: negative cycle (A < 0)
   if (pos[2] < 0.0) {
      B[0] *= -1.0;
      B[1] *= -1.0;
      B[2] *= -1.0;
   };
#endif
};

// Compute parallel component of diffusion
inline double GetKappaPara(double t, double* pos, double mom, double Bmag)
{
   double r = Norm(pos);
   double v = Vel(mom);
   double Omega = CyclotronFreq(v, Bmag);
   double l_b = 1.2618 * lc_0 * sqrt(r / r_0);
   double A2_sl = 0.2 * dB2_0 * Sqr(r_0 / r) / Sqr(Bmag);
   return 3.0 * Cube(v) / (20.0 * l_b * Sqr(Omega) * A2_sl * sin(3.0 * M_PI / 5.0)) * (1.0 + (72.0 / 7.0) * cbrt(Quint(l_b * Omega / v)));
};

// Compute perpendicular component of diffusion
inline double GetKappaPerp(double t, double* pos, double mom, double Bmag)
{
   double A2_2D = 0.8 * dB2_0 * Sqr(r_0 / Norm(pos)) / Sqr(Bmag);
   double kappa_para = GetKappaPara(t, pos, mom, Bmag);
   return eta * A2_2D * kappa_para;
};

// Compute all 9 components of the diffusion tensor
void GetKappaTensor(double t, double* pos, double mom, double* B, double Kappa[][3])
{
   int i, j, delta_ij;
   double k_para, k_perp, Bmag, Bmag2;

   Bmag = Norm(B);
   Bmag2 = Sqr(Bmag);
   k_para = GetKappaPara(t, pos, mom, Bmag);
   k_perp = GetKappaPerp(t, pos, mom, Bmag);

   for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
         delta_ij = (i == j);
         Kappa[i][j] = (k_perp * delta_ij + (k_para - k_perp) * B[i] * B[j] / Bmag2);
      };
   };
};

// Boundary condition
inline double OuterBoundaryCondition(double* pos, double mom)
{
   double T = EnrKin(mom) / T_0;
   return (12.0 / Sqr(mom)) * pow(T, -2.6) / (1.0 + 5.3 * pow(T, -1.22) + 1.3 * pow(T, -2.8) + 0.0087 * pow(T, -4.32));
};

// Set initial conditions
void InitializeTrajectory(double* pos_init, double& lnp_init)
{
// Initial position
   pos_init[0] = r_in;
   pos_init[1] = 0.0;
   pos_init[2] = 0.0;
// Initial momentum
   lnp_init = lnp_min + drand48() * (lnp_max - lnp_min);
};

// Integrate a single trajectory
void IntegrateTrajectory(double* pos_in, double lnp_in, double* pos_out, double& lnp_out, bool print)
{
// Declare variables and arrays
   int i, j, step = 0;
   double t, r, lnp, lnp_rate, mom, k_para, k_perp, dl, dt;
   static double rn[4];
   double pos[3], pos1[3], pos_rate[3], pos_rate2[3], basis[3][3];
   double u[3], B[3], u1[3], B1[3], Kappa[3][3], Kappa1[3][3], divK[3];
   double magB, vrL_3, curlB[3], gradB[3], gradBxB[3];
   std::ofstream trajectory_file;

// Initialize local variables
   t = 0.0;
   pos[0] = pos_in[0];
   pos[1] = pos_in[1];
   pos[2] = pos_in[2];
   lnp = lnp_in;
   r = Norm(pos);
   mom = exp(lnp);

// Open trajectory file if trajectory path is being printed
   if(print) {
      trajectory_file.open(trajectory_fname.c_str(), std::ofstream::out);
      trajectory_file << std::setprecision(6);
   };

// Time loop
   while(r > r_min && r < r_max) {
// Get fields and diffusion tensor
      GetFields(t, pos, u, B);
      magB = Norm(B);
      k_para = GetKappaPara(t, pos, mom, magB);
      k_perp = GetKappaPerp(t, pos, mom, magB);
      GetKappaTensor(t, pos, mom, B, Kappa);

// Field-aligned (FA) frame
      UnitVec(B, basis[2]);
      GetPerpVec(basis[2], basis[1]);
      Cross(basis[1], basis[2], basis[0]);
      
// Compute the divergence of K (numerical derivative)
      divK[0] = 0.0;
      divK[1] = 0.0;
      divK[2] = 0.0;
      for(j = 0; j < 3; j++) {
         for(i = 0; i < 3; i++) pos1[i] = pos[i];
         pos1[j] += delta;
         GetFields(t, pos1, u1, B1);
         GetKappaTensor(t, pos1, mom, B1, Kappa1);
         for(i = 0; i < 3; i++) divK[i] += (Kappa1[j][i] - Kappa[j][i]) / delta;
      };

// Advective terms
      pos_rate[0] = u[0] - divK[0];
      pos_rate[1] = u[1] - divK[1];
      pos_rate[2] = u[2] - divK[2];
      lnp_rate = -2.0 * u_0 / (3.0 * r);

#ifdef MAGNETIC_DRIFT
// Compute and add magnetic drifts (numerical derivative)
      curlB[0] = 0.0;
      curlB[1] = 0.0;
      curlB[2] = 0.0;
      gradB[0] = 0.0;
      gradB[1] = 0.0;
      gradB[2] = 0.0;
      for(j = 0; j < 3; j++) {
         for(i = 0; i < 3; i++) pos1[i] = pos[i];
         pos1[j] += delta;
         GetFields(t, pos1, u1, B1);
         gradB[j] = (Norm(B1) - magB) / delta;
         if (j == 0) {
            curlB[1] -= (B1[2] - B[2]) / delta;
            curlB[2] += (B1[1] - B[1]) / delta;
         } else if (j == 1) {
            curlB[0] += (B1[2] - B[2]) / delta;
            curlB[2] -= (B1[0] - B[0]) / delta;
         } else {
            curlB[0] -= (B1[1] - B[1]) / delta;
            curlB[1] += (B1[0] - B[0]) / delta;
         };
      };
      Cross(gradB, B, gradBxB);
      vrL_3 = mom * Vel(mom) * c_code / (3.0 * q_p * magB);

      pos_rate[0] += vrL_3 * (curlB[0] - 2.0 * gradBxB[0] / magB);
      pos_rate[1] += vrL_3 * (curlB[1] - 2.0 * gradBxB[1] / magB);
      pos_rate[2] += vrL_3 * (curlB[2] - 2.0 * gradBxB[2] / magB);
#endif

// Diffusion in the x' and y' direction (FA frame)
      rn[0] = drand48();
      rn[1] = drand48();
      pos_rate2[0] = pos_rate2[1] = 2.0 * sqrt(-k_perp * log(rn[0]));
      pos_rate2[0] *= cos(2.0 * M_PI * rn[1]);
      pos_rate2[1] *= sin(2.0 * M_PI * rn[1]);

// Diffusion in the z' direction (FA frame)
      if(!(step % 2)) {
         rn[2] = drand48();
         rn[3] = drand48();
         pos_rate2[2] = 2.0 * sqrt(-k_para * log(rn[2])) * cos(2.0 * M_PI * rn[3]);
      }
// Use two random numbers from the previous step for efficiency
      else {
         pos_rate2[2] = 2.0 * sqrt(-k_para * log(rn[2])) * sin(2.0 * M_PI * rn[3]);
      };
      ChangeFromBasis(basis, pos_rate2);

// Estimate the time step: advective space, diffusive (space), and advective momentum
      dl = 0.05 * r;
      dt = dl / Norm(pos_rate);
      dt = fmin(dt, Sqr(dl) / k_para);
      dt = fmin(dt, 0.05 / fabs(lnp_rate));
      dt *= cfl;

// Update the position and momentum (backward in time)
      t += dt;
      pos[0] -= pos_rate[0] * dt + pos_rate2[0] * sqrt(dt);
      pos[1] -= pos_rate[1] * dt + pos_rate2[1] * sqrt(dt);
      pos[2] -= pos_rate[2] * dt + pos_rate2[2] * sqrt(dt);
      lnp -= lnp_rate * dt;
      r = Norm(pos);
      mom = exp(lnp);
      step++;

// Print the projection on the XY plane
      if(print) {
         trajectory_file << std::setw(15) << pos[0]
                         << std::setw(15) << pos[1]
                         << std::setw(15) << pos[2]
                         << std::setw(15) << mom
                         << std::endl;
      };
   };

// Close trajectory file if necessary and save final position and momentum for binning
   if(print) trajectory_file.close();
   pos_out[0] = pos[0];
   pos_out[1] = pos[1];
   pos_out[2] = pos[2];
   lnp_out = lnp;
};

// Bin trajectory function
void BinTrajectory(double lnp_init, double* pos_final, double lnp_final, double* counts, double* distro)
{
   int bin;
   double mom_final = exp(lnp_final);

// Find momentum bin based on starting momentum and add one to the counts array
   bin = (lnp_init - lnp_min) / dlnp;
   counts[bin] += 1.0;
// Find momentum weight based on the final momentum and add that value to the distribution array
   distro[bin] += (Norm(pos_final) >= r_max ? OuterBoundaryCondition(pos_final, mom_final) : 0.0);
};

int main(void)
{
// Initialize the MPI environment
   MPI_Init(NULL, NULL);
   
// Declare local variables and arrays for computations
   int part, bin;
   double mom, lnp_in, lnp_out;
   double pos_in[3], pos_out[3];
   time_t time_start, time_end;

// Find MPI communicator quantities
   int comm_rank, comm_size;
   MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

// Initialize distribution arrays
   double* counts_partial;
   double* distro_partial;
   double* counts = new double[nbins];
   double* distro = new double[nbins];
   memset(counts, 0, nbins * sizeof(double));
   memset(distro, 0, nbins * sizeof(double));
   srand48(time(NULL)+comm_rank);

// Assign tasks to workers
   int ntraj_proc;
   if(comm_rank == 0) {
      time_start = time(NULL);
      ntraj_proc = ntraj / comm_size;
   };
   MPI_Bcast(&ntraj_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

// Loop over all trajectories
   for(part = 0; part < ntraj_proc; part++) {
// Initialize trajectory
      InitializeTrajectory(pos_in, lnp_in);
// Integrate trajectory
      IntegrateTrajectory(pos_in, lnp_in, pos_out, lnp_out, false);
// Bin trajectory
      BinTrajectory(lnp_in, pos_out, lnp_out, counts, distro);
   };

// Gather results from all workers
   if(comm_rank == 0) {
      int cpu;
      counts_partial = new double[nbins];
      distro_partial = new double[nbins];
      
      for(cpu = 1; cpu < comm_size; cpu++) {
         MPI_Recv(counts_partial, nbins, MPI_DOUBLE, cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         MPI_Recv(distro_partial, nbins, MPI_DOUBLE, cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         for(bin = 0; bin < nbins; bin++) {
            counts[bin] += counts_partial[bin];
            distro[bin] += distro_partial[bin];
         };
      };

      time_end = time(NULL);
// Output total number of trajectories simulated as a sanity check
      double total_counts = 0.0;
      for(bin = 0; bin < nbins; bin++) total_counts += counts[bin];
      std::cerr << (int)(total_counts) << " trajectories completed\n";
      std::cerr << "Runtime was " << (int)(time_end - time_start) << "s\n";
   }
   else {
      MPI_Send(counts, nbins, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(distro, nbins, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
   };


   if(comm_rank == 0) {
// Print one trajectory
#ifdef PRINT_TRAJECTORY
      const double T_traj = 500.0 * MeV_cgs;
      IntegrateTrajectory(pos_in, log(Mom(T_traj)), pos_out, lnp_out, true);
#endif

// Print the spectrum
#ifdef PRINT_INTENSITY
      pos_out[0] = r_max;
      pos_out[1] = 0.0;
      pos_out[2] = 0.0;
      std::ofstream intensity_file(intensity_fname.c_str(), std::ofstream::out);
      intensity_file << std::setprecision(6);
      for(bin = 0; bin < nbins; bin++) {
         mom = exp(lnp_min + (bin + 0.5) * dlnp);
         intensity_file << std::setw(15) << EnrKin(mom) / MeV_cgs
                        << std::setw(15) << Sqr(mom) * (counts[bin] >= 0.999 ? distro[bin] / counts[bin] : 0.0)
                        << std::setw(15) << Sqr(mom) * OuterBoundaryCondition(pos_out, mom)
                        << std::endl;
      };
      intensity_file.close();
#endif
   };

// Clean up
   delete[] counts;
   delete[] distro;
   if(comm_rank == 0) {
      delete[] counts_partial;
      delete[] distro_partial;
   };

// Finalize the MPI environment.
   MPI_Finalize();
};

