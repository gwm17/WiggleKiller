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



class CubicSpline
{
public:

    // struct holding spline info
    struct Spline
    {
        double y1 = 0, y2 = 0;
        double x1 = 0, x2 = 0;
        double k1 = 0, k2 = 0;

        //Evaluate function within this spline
        const double Evaluate(double x) const
        {
            double t = (x - x1) / (x2 - x1);
            double a = k1 * (x2 - x1) - (y2 - y1);
            double b = -k2 * (x2 - x1) + (y2 - y1);
            return (1.0 - t) * y1 + t * y2 + t * (1.0 - t) * ((1.0 - t) * a + t * b);
        }
    };

    CubicSpline();
    CubicSpline(const std::string &filename);
    ~CubicSpline();

    void ReadFile(const std::string &filename);
    void ReadData(const std::vector<double> &x, const std::vector<double>& y)
    {
        m_isValid = true;
        MakeSplines(x, y);
    };

    double Evaluate(double x);
    double EvaluateROOT(double *x, double *p); // for plotting as ROOT function

    const bool IsValid() const { return m_isValid; };

private:
    void MakeSplines(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<Spline> m_splineList;

    bool m_isValid;
};

#endif