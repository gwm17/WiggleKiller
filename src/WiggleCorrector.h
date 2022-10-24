/*
WiggleCorrector.h
Class for correcting focal plane spectra for differential-nonlinearity using an existing solution. Reads in calibration from file
and generates the cubic spline for interpolation. Appies to the data set. Correction can generates histograms as well as an optional
tree of corrected data. Note: requires particle ID (scint/cathode) and x1-x2 correltaion cuts. Resulting data has cuts applied.

Gordon M. May 2021
*/
#ifndef WIGGLECORRECTOR_H
#define WIGGLECORRECTOR_H

#include <string>
#include <THashTable.h>
#include "CubicSpline.h"
#include "DataStructs.h"
#include "CutHandler.h"
#include "FP_kinematics.h"

class WiggleCorrector
{
public:
    WiggleCorrector();
    WiggleCorrector(const std::string& path, const std::string& tag, double lambdaFront, double lambdaBack);
    void ReadCorrectionData(const std::string& path, const std::string& tag, double lambdaFront, double lambdaBack);
    void SetCuts(const std::string& edename, const std::string& xxname) { cuts.SetCuts(edename, xxname); }
    void ApplyCorrection(const std::string& inname, const std::string& outname, bool writeTree=true);

    void SetFPParameters(const FPParameters& params) { m_fpParams = params; }

private:
    void CopyBulkEvent(ProcessedEvent *inevent, ProcessedEvent &outevent);
    bool CheckSplines()
    { 
        return (frontLeftSpline.IsValid() && frontRightSpline.IsValid() && backLeftSpline.IsValid() && backRightSpline.IsValid()); 
    }

    FPParameters m_fpParams;

    CubicSpline frontLeftSpline, frontRightSpline, backLeftSpline, backRightSpline;
    CutHandler cuts;
};

#endif