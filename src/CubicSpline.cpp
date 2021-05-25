/*
CubicSpline.cpp
Class for generating cubic splines of data in tables or held in memory. Cubic splines are a form of interpolation,
somewhat more advanced than linear interpolation, but not so complicated that it significantly slows down a calculation.
For more information see Wikipedia or Num. Rec. in C.

Gordon M. May 2021
*/
#include "CubicSpline.h"
#include <fstream>
#include <iostream>
#include <cmath>

CubicSpline::CubicSpline() :
	validFlag(false)
{
}

CubicSpline::CubicSpline(std::string& filename) :
	validFlag(false)
{
	ReadFile(filename);
}

CubicSpline::~CubicSpline() {}

/*
Expected file format is a naked (no header) single space separated table:
x y\n
*/
void CubicSpline::ReadFile(std::string& filename) {
	std::ifstream input(filename);
	if(!input.is_open()) {
		std::cerr<<"Unable to open input data at CubicSpline::ReadFile from filename: "<<filename<<std::endl;
		validFlag = false;
		return;
	}

	std::string junk;
	//std::getline(input, junk);

	double x, y;
	while(input>>x) {
		input>>y;
		data_x.push_back(x);
		data_y.push_back(y);
	}

	if(data_x.size() != data_y.size()) {
		std::cerr<<"Error in CubicSpline::ReadFile! Number of x points not equal to number of y points!"<<std::endl;
		validFlag = false;
		return;
	}

	validFlag = true;

	input.close();

	MakeSplines();
}

/*
After data is read in splines can be solved. Each data point is referred to as a knot. Endpoint conditions are
derivatives of neighbors must be equal and second derivatives must be zero (natural cubic splines). Solved using gaussian elimination.
*/
void CubicSpline::MakeSplines() {
	if(!validFlag) {
		std::cerr<<"Error at CubicSpline::MakeSplines! Unable to generate splines without first initializing data."<<std::endl;
		return;
	}

	int knots = data_x.size();

	//Matrix and vector data init
	double** a = new double*[knots];
	for(int i=0; i<knots; i++) {
		a[i] = new double[knots];
	}
	double* b = new double[knots];
	double* k = new double[knots];

	Spline s;
	splines.clear();
	double x0, x1, x2;
	double y0, y1, y2;

	//Setup matrix eqn.
	for(int i=0; i<knots; i++) {
		if(i == 0) {
			x1 = data_x[i];
			x2 = data_x[i+1];
			y1 = data_y[i];
			y2 = data_y[i+1];
			a[i][i] = 2.0/(x2-x1);
			a[i][i+1] = 1.0/(x2-x1);
			b[i] = 3.0*(y2-y1)/(std::pow((x2-x1), 2.0));
      		s.x1 = x1; s.x2 = x2;
      		s.y1 = y1; s.y2 = y2;
      		splines.push_back(s);
		} else if(i == (knots-1)) {
			x0 = data_x[i-1];
			x1 = data_x[i];
      		y0 = data_y[i-1];
      		y1 = data_y[i];
      		a[i][i-1] = 1.0/(x1-x0);
      		a[i][i] = 2.0/(x1-x0);
      		b[i] = 3.0*(y1-y0)/(std::pow((x1-x0), 2.0));
		} else {
	    	x0 = data_x[i-1]; 
	    	x1 = data_x[i]; 
	    	x2 = data_x[i+1];
      		y0 = data_y[i-1];
      		y1 = data_y[i];
      		y2 = data_y[i+1];
      		a[i][i-1] = 1.0/(x1-x0);
      		a[i][i] = 2.0/(x1-x0)+2.0/(x2-x1);
      		a[i][i+1] = 1.0/(x2-x1);
      		b[i] = 3.0*(y1-y0)/(std::pow((x1-x0), 2.0))+3.0*(y2-y1)/(std::pow((x2-x1), 2.0));
      		s.x1 = x1; s.x2 = x2;
      		s.y1 = y1; s.y2 = y2;
      		splines.push_back(s);
		}
	}

	//solve for curvature vector k using gaussian elimination
	a[0][1] /= a[0][0];
	b[0] /= a[0][0];
	for(int i=1; i<(knots-1); i++) {
		a[i][i+1] /= a[i][i]-a[i][i-1]*a[i-1][i];
   	 	b[i] = (b[i] - a[i][i-1]*b[i-1])/(a[i][i]-a[i][i-1]*a[i-1][i]);
	}
	int g1 = knots-1;
  	int g2 = knots-2;
  	b[g1] = (b[g1]-a[g1][g2]*b[g2])/(a[g1][g1]-a[g1][g2]*a[g2][g1]);

  	k[g1] = b[g1];
  	for(int i=(knots-2); i>=0; i--) {
    	k[i] = b[i] - a[i][i+1]*k[i+1];
  	}

  	//Fill the spline data
  	for(unsigned int i=0; i<splines.size(); i++) {
  		splines[i].k1 = k[i];
  		splines[i].k2 = k[i+1];
  	}

  	//deallocate
  	delete[] b;
  	delete[] k;
  	for(int i=0; i<knots; i++) {
  		delete[] a[i];
  	}
  	delete[] a;
}

double CubicSpline::Evaluate(double x) {
	if(!validFlag) {
		std::cerr<<"Error at CubicSpline::Evaluate! Unable to evaluate without first generating splines."<<std::endl;
		return 0.0;
	}

	Spline s;
	for(unsigned int i=0; i<splines.size(); i++) {
		auto& spline = splines[i];
		if(x >= spline.x1 && x < spline.x2) {
			s = spline;
			break;
		} else if (i == (splines.size() -1)) {
			//std::cerr<<"Error at CubicSpline::Evaluate! Input x value: "<<x<<" is not within the spline range min: "<<splines[0].x1<<" max: "<<splines[splines.size()-1].x2<<std::endl;
			return 0.0;
		}
	}

	double t = (x-s.x1)/(s.x2-s.x1);
    double a = s.k1*(s.x2-s.x1)-(s.y2-s.y1);
    double b = -s.k2*(s.x2-s.x1)+(s.y2-s.y1);
    return (1.0-t)*s.y1+t*s.y2+t*(1.0-t)*((1.0-t)*a+t*b);
}

//Purely for plotting in ROOT, do not use for caluculations. 
double CubicSpline::EvaluateROOT(double* x, double* p) {
	if(!validFlag) {
		std::cerr<<"Error at CubicSpline::Evaluate! Unable to evaluate without first generating splines."<<std::endl;
		return 0.0;
	}

	double xval = x[0];
	Spline s;
	for(unsigned int i=0; i<splines.size(); i++) {
		auto& spline = splines[i];
		if(xval >= spline.x1 && xval < spline.x2) {
			s = spline;
			break;
		} else if (i == (splines.size() -1)) {
			//std::cerr<<"Error at CubicSpline::Evaluate! Input x value: "<<x<<" is not within the spline range min: "<<splines[0].x1<<" max: "<<splines[splines.size()-1].x2<<std::endl;
			return 0.0;
		}
	}

	double t = (xval-s.x1)/(s.x2-s.x1);
    double a = s.k1*(s.x2-s.x1)-(s.y2-s.y1);
    double b = -s.k2*(s.x2-s.x1)+(s.y2-s.y1);
    return (1.0-t)*s.y1+t*s.y2+t*(1.0-t)*((1.0-t)*a+t*b);
}

