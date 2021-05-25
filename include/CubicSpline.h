/*
CubicSpline.h 
Class for generating cubic splines of data in tables or held in memory. Cubic splines are a form of interpolation,
somewhat more advanced than linear interpolation, but not so complicated that it significantly slows down a calculation.
For more information see Wikipedia or Num. Rec. in C.

Gordon M. May 2021
*/
#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include <vector>
#include <string>

//struct holding final spline info
struct Spline {
	double y1=0, y2=0;
	double x1=0, x2=0;
	double k1=0, k2=0;
};

class CubicSpline {
public:
	CubicSpline();
	CubicSpline(std::string& filename);
	~CubicSpline();
	void ReadFile(std::string& filename);
	void ReadData(std::vector<double>& x, std::vector<double>& y) { 
		data_x = x; 
		data_y = y; 
		validFlag=true;
		MakeSplines();
	};
	bool IsValid() { return validFlag; };
	double Evaluate(double x);
	double EvaluateROOT(double* x, double* p); //for plotting as ROOT function

private:
	void MakeSplines();
	std::vector<double> data_x;
	std::vector<double> data_y;
	std::vector<Spline> splines;

	bool validFlag;
};

#endif