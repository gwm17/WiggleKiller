/*
WiggleKiller.h
Class for correcting focal plane spectra for differential-nonlinearity using a new solution or optimizing a solution. Can be used to simply test a solution
(guess values of lambda) or optimize a guess using a simluated annealing method. Can generate calibration files for use with the WiggleCorrector class. Note:
requires particle ID (scint/cathode) and x1-x2 correlation cuts.

Gordon M. May 2021
*/

#ifndef WIGGLEKILLER_H
#define WIGGLEKILLER_H

#include <string>
#include <vector>
#include <iostream>
#include <TF1.h>
#include <TTree.h>
#include "CutHandler.h"
#include "DataStructs.h"
#include "CubicSpline.h"
#include "RandomGenerator.h"
#include "FP_kinematics.h"

// parsed down information for optimization
struct CorrectionPoint
{
    double leftPSD;
    double rightPSD;
    double leftTime;
    double rightTime;
};

// simple structure for holding calibration data
struct Calibration
{
    double time_i;
    std::vector<double> deltas;

    double GetDelta()
    {
        if (deltas.size() == 0)
            return 0.0;

        double delta = 0.0;
        for (auto d : deltas)
        {
            delta += d;
        }
        return delta / deltas.size();
    }
};

struct OptimizationParamters
{
    double lambdaFront = 0.0;
    double lambdaBack = 0.0;
};



class WiggleKiller
{
public:
    WiggleKiller(double lambdaFront, double lambdaBack);
    ~WiggleKiller();
    void SetCuts(const std::string& edename, const std::string& xxname);
    void OptimizeParameters(const std::string& inputFileName, const std::string& outputFileName);
    void ApplyLambdas(const std::string& inputFileName, const std::string& outputFileName);
    void GenerateCalibrationFiles(const std::string& path, const std::string& calFileNameTag);

    void SetFPParameters(const FPParameters& params) { m_fpParams = params; }

private:

    double GetX1Peak2Trough(double lambda);
    double GetX2Peak2Trough(double lambda);

    OptimizationParamters MinimizeSimAnneal();
    double ApplyPSDCorrection(double time, double psd, int delay, double lambda);

    OptimizationParamters m_optParams;
    XWeights m_xWeights;
    FPParameters m_fpParams;

    std::normal_distribution<double> m_gaussian{0.0, 1.0};
    std::uniform_real_distribution<double> m_uniform{0.0, 1.0};

    std::vector<std::vector<double>> delayBinnedAvgPSD;
    std::vector<std::vector<double>> delayBinnedCounts;

    std::vector<CubicSpline> PSDsplines;

    std::vector<CorrectionPoint> frontPoints;
    std::vector<CorrectionPoint> backPoints;
    std::vector<std::vector<Calibration>> calibData;

    std::vector<double> widths1;
    std::vector<double> widths2;
    std::vector<double> times;

    // intrinsic binning information
    static constexpr double BINWIDTH = 0.5;
    static constexpr double BINMAX = 1500.0;
    static constexpr double BINMIN = 0.0;
    int nbins;

    CutHandler cuts;

    // Masks for quick indexing
    enum Delays
    {
        FrontRight,
        FrontLeft,
        BackRight,
        BackLeft,
        AnodeFront,
        AnodeBack
    };
};

#endif