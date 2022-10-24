/*

  Functions for the calculation of the kinematic shift of the FP
  for the SESPS @ FSU.

>>>  Delta_Z(int...) returns the shift of the FP in the z-direction in
     cm. A negative (<0) delta-z is defined as a shift towards the
     magnet.

     Arguments: Delta_Z(int ZT, int AT, int ZP, int AP, int ZE, int AE,
	                double EP, double angle, double B),
        where Z,A are atomic number and mass number for each particle,
        EP is the KE of the projectile (i.e. beam energy in MeV),
        angle is the spectrograph angle (in degrees),
        B is the field in Gauss.


>>>  Wire_Dist() returns the distance (in cm) between the wires for
     the calculation of relative weights between FP1 and FP2.

  //format: T(P,E)R
  //   i.e., T = target,
  //         P = projectile,
  //         E = ejectile,
  //         R = residual;
  //expects angle in degs, B in G, masses and KE in MeV

  KGH -- Jul19

  More clean up. Restructure a bit. GWM -- Oct. 2022

*/

#ifndef FP_KINEMATICS
#define FP_KINEMATICS

#include <cstdint>

//requires (Z,A) for T, P, and E, as well as energy of P,
// spectrograph angle of interest, and field value

struct FPParameters
{
    uint64_t ZT = 0; //Target nucleus Z
    uint64_t AT = 0; //Target nucleus A
    uint64_t ZP = 0; //Projectile nucleus Z
    uint64_t AP = 0; //Projectile nucleus A
    uint64_t ZE = 0; //Ejectile nucleus Z
    uint64_t AE = 0; //Ejectile nucleus A
    double energyP = 0.0; //Projectile Kinetic Energy (MeV)
    double spsAngle = 0.0; //SPS Lab angle (degrees)
    double bField = 0.0; //SPS Magnetic Field (G)
};

struct XWeights
{
    double x1Weight = 0.0;
    double x2Weight = 0.0;
};

double Delta_Z(const FPParameters& params); 

static constexpr double Wire_Dist() { return 4.28625; }

XWeights GetXWeights(const FPParameters& params);

#endif
