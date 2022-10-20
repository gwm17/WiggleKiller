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
#include <THashTable.h>
#include <TTree.h>
#include <TRandom3.h>
#include "CutHandler.h"
#include "DataStructs.h"
#include "CubicSpline.h"

//parsed down information for optimization
struct CorrectionPoint {
	double leftPSD;
	double rightPSD;
	double leftTime;
	double rightTime;
};

//simple structure for holding calibration data
struct Calibration {
	double time_i;
	std::vector<double> deltas;

	double GetDelta() {
		if(deltas.size() == 0) return 0.0;

		double delta = 0.0;
		for(auto d : deltas) {
			delta += d; 
		}
		return delta/deltas.size();
	}
};

class WiggleKiller {
public:
	WiggleKiller();
	WiggleKiller(std::string& edename, std::string& xxname, double lf, double lb);
	~WiggleKiller();
	void SetCuts(std::string& edename, std::string& xxname);
	void SetLambdas(double lf, double lb);
	void OptimizeParameters(std::string& inname, std::string& outname);
	void ApplyLambdas(std::string& inname, std::string& outname);
	void GenerateCalibrationFiles(std::string& dflname, std::string& dfrname, std::string& dblname, std::string& dbrname);

private:
	void MyFill(THashTable* table, std::string name, std::string titlex, int binsx, double minx, double maxx, double valuex);
	void MyFill(THashTable* table, std::string name, std::string titlex, std::string titley, int binsx, double minx, double maxx, double valuex, int binsy, double miny, double maxy, double valuey);
	void MyClear(THashTable* table, std::string name);
	std::pair<double, double> GetXavgWeights(int zt, int at, int zp, int ap, int ze, int ae, double bke, double theta, double b);
	double GetTCheck1Width(double lambda);
	double GetTCheck2Width(double lambda);
	double GetX1Peak2Trough(double lambda);
	double GetX2Peak2Trough(double lambda);
	double GetThetaWidth(double lambda1, double lambda2);
	std::pair<double,double> MinimizeSimAnneal();
	double ApplyPSDCorrection(double time, double psd, int delay, double lambda);

	double lambdaFront;
	double lambdaBack;
	std::pair<double, double> weights;

	TRandom3* generator;

	std::vector<std::vector<double>> delayBinnedAvgPSD;
	std::vector<std::vector<double>> delayBinnedCounts;

	std::vector<CubicSpline> PSDsplines;

	std::vector<CorrectionPoint> frontPoints;
	std::vector<CorrectionPoint> backPoints;
	std::vector<std::vector<Calibration>> calibData;

	std::vector<double> widths1;
	std::vector<double> widths2;
	std::vector<double> times;

	static constexpr double GOLDEN_RATIO = (1.0 + std::sqrt(5))/2.0; //UNUSED
	//intrinsic binning information
	static constexpr double BINWIDTH = 0.5;
	static constexpr double BINMAX = 1500.0;
	static constexpr double BINMIN = 0.0;
	int nbins;

	CutHandler cuts;
	
	//Masks for quick indexing
	enum Delays {
		FrontRight,
		FrontLeft,
		BackRight,
		BackLeft,
		AnodeFront,
		AnodeBack
	};
	
};

#endif