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

  Small modifications for use with the MassLookup class GWM -- Jan 2021

  More clean up. Restructure a bit. GWM -- Oct. 2022

*/
#include <iostream>
#include <cmath>
#include "MassLookup.h"
#include "FP_kinematics.h"

// requires (Z,A) for T, P, and E, as well as energy of P,
//  spectrograph angle of interest, and field value
double Delta_Z(const FPParameters& params)
{

    /* CONSTANTS */
    static constexpr double s_MeV2J = 1.60218E-13;               // J per MeV
    static constexpr double s_unitCharge = 1.602E-19;            // Coulombs
    static constexpr double s_c = 2.9979E8;                       // m/s
    /* SESPS-SPECIFIC */
    static constexpr double s_dispersion = 1.96; // dispersion (x/rho)
    static constexpr double s_magnification = 0.39;  // magnification in x
    static constexpr double s_deg2rad = M_PI / 180.;

    uint64_t ZR = params.ZT + params.ZP - params.ZE;
    uint64_t AR = params.AT + params.AP - params.AE;

    double bField_tesla = params.bField/10000; // convert to tesla
    double spsAngle_rads = params.spsAngle * s_deg2rad;

    const MassLookup& masses = MassLookup::GetInstance();
    double massT = masses.FindMass(params.ZT, params.AT);
    double massP = masses.FindMass(params.ZP, params.AP);
    double massE = masses.FindMass(params.ZE, params.AE);
    double massR = masses.FindMass(ZR, AR);

    if ((massT * massP * massE * massR) == 0)
    {
        std::cerr << "***WARNING: error loading one or more masses; returning 0\n";
        return 0;
    }

    double Q = massT + massP - massE - massR; // Q-value

    // kinematics a la Iliadis p.590
    double term1 = sqrt(massP * massE * params.energyP) / (massE + massR) * std::cos(spsAngle_rads);
    double term2 = (params.energyP * (massR - massP) + massR * Q) / (massE + massR);

    double energyE = term1 + sqrt(term1 * term1 + term2);
    energyE *= energyE;

    // ejectile momentum
    double pE = std::sqrt(energyE * (energyE + 2.0 * massE));

    // calculate rho from B a la B*rho = (proj. momentum)/(proj. charge)
    double rho = (pE * s_MeV2J) / (params.ZE * s_unitCharge * s_c * bField_tesla) * 100; // in cm

    double K = std::sqrt(massP * massE * params.energyP / energyE) * std::sin(spsAngle_rads);

    double denom = massE + massR - sqrt(massP * massE * params.energyP / energyE) * std::cos(spsAngle_rads);

    K /= denom;
    return -1.0 * rho * s_dispersion * s_magnification * K; // delta-Z in cm
}

XWeights GetXWeights(const FPParameters &params)
{
    XWeights weights;
    double dz = Delta_Z(params);
    weights.x1Weight = ((Wire_Dist() / 2.0) - dz) / Wire_Dist();
    weights.x2Weight= 1.0 - weights.x1Weight;
    return weights;
}