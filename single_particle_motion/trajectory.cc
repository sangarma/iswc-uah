// Include packages/dependencies here
#include <iostream>
#include <iomanip>

#include "trajectory.hh"

// Function for magnetic field
inline void Bfield(double t, double* pos, double* B)
{
// Constant field
   if (pos[0] <= 0){
   B[0] = 0.0;
   B[1] = 0.0;
   B[2] = 0.4;
   } else{
      B[0] = 0;
      B[1] = 0;
      B[2] = 0.8;
   }
}

// Function for electric field
inline void Efield(double t, double* pos, double* E)
{
// Constant field
	E[0] = 0;
	E[1] = 0.01;
	E[2] = 0;
};

// Lorentz force function
inline void Lorentz(double t, double* pos, double* vel, double* force)
{
   double B[3], E[3];

// Magnetic (v x B) force
   Bfield(t, pos, B);
   Cross(vel, B, force); 
   force[0] *= q_c;
   force[1] *= q_c;
   force[2] *= q_c;

// Electric force
   Efield(t, pos, E);
   force[0] += q*E[0];
   force[1] += q*E[1];
   force[2] += q*E[2];
};

// Initial condions function
inline void SetInitialConditions(double* pos, double* mom)
{
// Position
   pos[0] = Rx;
   pos[1] = Ry;
   pos[2] = Rz;
// Momentum
   mom[0] = Px;
   mom[1] = Py;
   mom[2] = Pz;
};

// Timestep function
inline double TimeStep(double t, double* pos, double* mom)
{
   double B[3];
   Bfield(t, pos, B);
   double rg = LarmorRadius(Norm(mom),Norm(B));
   double dr = 2.0 * rg * M_PI / 20.0;
// Two constraints: physical gyro-radius and absolute maximum distance traveled
   return CFL * fmin(dr, dmax) / Vel(Norm(mom));
};

// Advance position and momentum using Euler's method
void AdvanceEuler(double t, double* pos, double* mom, double dt)
{
   double vel[3], F[3];

// Find velocity (position slope) and force (momentum slope)
   Vel(mom, vel);
   Lorentz(t, pos, vel, F);
   
// Advance position and momentum to next step
   pos[0] += vel[0] * dt;
   pos[1] += vel[1] * dt;
   pos[2] += vel[2] * dt;
   mom[0] += F[0] * dt;
   mom[1] += F[1] * dt;
   mom[2] += F[2] * dt;
};

// Advance position and momentum using Heun's method
void AdvanceHeun(double t, double* pos, double* mom, double dt)
{
// TODO
   double vel1[3], vel2[3], F1[3], F2[3], pos_temp[3], mom_temp[3];
   Vel(mom,vel1);
   Lorentz(t,pos,vel1,F1);
   pos_temp[0] = pos[0] + vel1[0] * dt;
   pos_temp[1] = pos[1] + vel1[1] * dt;
   pos_temp[2] = pos[2] + vel1[2] * dt;
   mom_temp[0] = mom[0] + F1[0] * dt;
   mom_temp[1] = mom[1] + F1[1] * dt;
   mom_temp[2] = mom[2] + F1[1] * dt;


   Vel(mom_temp,vel2);
   Lorentz(t+dt,pos_temp,vel2,F2);
   pos[0] += (vel1[0] + vel2[0]) * dt/2;
   pos[1] += (vel1[1] + vel2[1]) * dt/2;
   pos[2] += (vel1[2] + vel2[2]) * dt/2;
   mom[0] += (F1[0] + F2[0]) * dt/2;
   mom[1] += (F1[1] + F2[1]) * dt/2;
   mom[2] += (F1[2] + F2[2]) * dt/2;


};

// Advance position and momentum using RK4 method
void AdvanceRK4(double t, double* pos, double* mom, double dt)
{
// TODO
   double vel1[3], vel2[3], vel3[3], vel4[3], F1[3], F2[3], F3[3], F4[3], pos_temp[3], mom_temp[3];
   Vel(mom,vel1);
   Lorentz(t,pos,vel1,F1);
   for(int i = 0; i<3; i++){
      pos_temp[i] = pos[i] + vel1[i]*dt/2.0;
      mom_temp[i] = mom[i] + F1[i] * dt/2.0;
   }

   Vel(mom_temp,vel2);
   Lorentz(t+dt,pos_temp,vel2,F2);

   for(int i = 0; i<3; i++){
   pos_temp[i] = pos[i] + vel2[i]*dt/2.0;
   mom_temp[i] = mom[i] + F2[i] * dt/2.0;
   } 

   Vel(mom_temp,vel3);
   Lorentz(t+dt,pos_temp,vel3,F3);

   for(int i = 0; i<3; i++){
   pos_temp[i] = pos[i] + vel3[i]*dt;
   mom_temp[i] = mom[i] + F3[i] * dt;
   }
   Vel(mom_temp,vel4);
   Lorentz(t+dt,pos_temp,vel4,F4);
   for(int i = 0; i<3; i++){
      pos[i] += (vel1[i] + 2*vel2[i] + 2* vel3[i] + vel4[i])*dt/6.0;
      mom[i] += (F1[i] + 2*F2[i] + 2* F3[i] + F4[i])*dt/6.0;
   }




};

int main(int argc, char** argv)
{
   double t, dt, pos[3], mom[3];

// Set your initial position and momentum
   t = 0.0;
   SetInitialConditions(pos, mom);

// Iterate over time
   std::setprecision(8);
   while(t < tmax) {
// Output trajectory
      std::cerr << "\r" << std::setw(16) << t;
      std::cout << std::setw(16) << t
                << std::setw(16) << pos[0]
                << std::setw(16) << pos[1]
                << std::setw(16) << pos[2]
                << std::setw(16) << mom[0]
                << std::setw(16) << mom[1]
                << std::setw(16) << mom[2]
                << std::endl;

// Figure out a suitable timestep
      dt = TimeStep(t, pos, mom);

// Advance one step using the uncommented method
      AdvanceRK4(t, pos, mom, dt);
      t += dt;
   };
   std::cerr << "\r" << std::setw(16) << t << std::endl;

   return 0;
};
