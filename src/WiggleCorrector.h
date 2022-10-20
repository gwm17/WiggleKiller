/*
WiggleCorrector.h
Class for correcting focal plane spectra for differential-nonlinearity using an existing solution. Reads in calibration from file and generates the cubic spline
for interpolation. Appies to the data set. Correction can generates histograms as well as an optional tree of corrected data. Note: requires particle ID 
(scint/cathode) and x1-x2 correltaion cuts. Resulting data has cuts applied. 

Gordon M. May 2021
*/
#ifndef WIGGLECORRECTOR_H
#define WIGGLECORRECTOR_H

#include <string>
#include <THashTable.h>
#include "CubicSpline.h"
#include "DataStructs.h"
#include "CutHandler.h"

class WiggleCorrector {
public:
	WiggleCorrector();
	WiggleCorrector(std::string& flname, std::string& frname, std::string& blname, std::string& brname);
	void ReadCorrectionData(std::string& flname, std::string& frname, std::string& blname, std::string& brname);
	void SetCuts(std::string& edename, std::string& xxname) { cuts.SetCuts(edename, xxname); };
	void ApplyCorrection(std::string& inname, std::string& outname);
	void ApplyCorrection_NoTree(std::string& inname, std::string& outname);

private:
	void MyFill(THashTable* table, std::string name, std::string titlex, int binsx, double minx, double maxx, double valuex);
	void MyFill(THashTable* table, std::string name, std::string titlex, std::string titley, int binsx, double minx, double maxx, double valuex, int binsy, double miny, double maxy, double valuey);
	void CopyBulkEvent(ProcessedEvent* inevent, ProcessedEvent& outevent);
	bool CheckSplines() {return (frontLeftSpline.IsValid() && frontRightSpline.IsValid() && backLeftSpline.IsValid() && backRightSpline.IsValid()); };
	std::pair<double, double> GetXavgWeights(int zt, int at, int zp, int ap, int ze, int ae, double bke, double theta, double b);

	CubicSpline frontLeftSpline, frontRightSpline, backLeftSpline, backRightSpline;
	CutHandler cuts;

};

#endif