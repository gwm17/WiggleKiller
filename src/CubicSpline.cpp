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
    m_isValid(false)
{
}

CubicSpline::CubicSpline(const std::string &filename) :
    m_isValid(false)
{
    ReadFile(filename);
}

CubicSpline::~CubicSpline() {}

/*
Expected file format is a naked (no header) single space separated table:
x y\n
*/
void CubicSpline::ReadFile(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input.is_open())
    {
        std::cerr << "Unable to open input data at CubicSpline::ReadFile from filename: " << filename << std::endl;
        m_isValid = false;
        return;
    }

    std::string junk;
    // std::getline(input, junk);

    double x, y;
    std::vector<double> data_x, data_y;
    while (input >> x)
    {
        input >> y;
        data_x.push_back(x);
        data_y.push_back(y);
    }

    if (data_x.size() != data_y.size())
    {
        std::cerr << "Error in CubicSpline::ReadFile! Number of x points not equal to number of y points!" << std::endl;
        m_isValid = false;
        return;
    }

    m_isValid = true;

    input.close();

    MakeSplines(data_x, data_y);
}

/*
After data is read in splines can be solved. Each data point is referred to as a knot. Endpoint conditions are
derivatives of neighbors must be equal and second derivatives must be zero (natural cubic splines). Solved using gaussian elimination.
*/
void CubicSpline::MakeSplines(const std::vector<double>& x, const std::vector<double>& y)
{
    if (!m_isValid)
    {
        std::cerr << "Error at CubicSpline::MakeSplines! Unable to generate splines without first initializing data." << std::endl;
        return;
    }

    std::size_t knots = x.size();

    // Matrix and vector data init
    double **a = new double *[knots];
    for (int i = 0; i < knots; i++)
    {
        a[i] = new double[knots];
    }
    double *b = new double[knots];
    double *k = new double[knots];

    Spline spline;
    m_splineList.clear();
    double x0, x1, x2;
    double y0, y1, y2;

    //Setup matrix equations

    //for first step in spline
    x1 = x[0];
    x2 = x[1];
    y1 = y[0];
    y2 = y[1];
    a[0][0] = 2.0 / (x2 - x1);
    a[0][1] = 1.0 / (x2 - x1);
    b[0] = 3.0 * (y2 - y1) / (std::pow((x2 - x1), 2.0));
    spline.x1 = x1;
    spline.x2 = x2;
    spline.y1 = y1;
    spline.y2 = y2;
    m_splineList.push_back(spline);

    // Setup non-edge matrix elements
    std::size_t lastIndex = knots - 1;
    for (std::size_t i = 1; i < lastIndex; i++)
    {
        x0 = x[i - 1];
        x1 = x[i];
        x2 = x[i + 1];
        y0 = y[i - 1];
        y1 = y[i];
        y2 = y[i + 1];
        a[i][i - 1] = 1.0 / (x1 - x0);
        a[i][i] = 2.0 / (x1 - x0) + 2.0 / (x2 - x1);
        a[i][i + 1] = 1.0 / (x2 - x1);
        b[i] = 3.0 * (y1 - y0) / (std::pow((x1 - x0), 2.0)) + 3.0 * (y2 - y1) / (std::pow((x2 - x1), 2.0));
        spline.x1 = x1;
        spline.x2 = x2;
        spline.y1 = y1;
        spline.y2 = y2;
        m_splineList.push_back(spline);
    }

    //For last step in spline
    x0 = x[lastIndex - 1];
    x1 = x[lastIndex];
    y0 = y[lastIndex - 1];
    y1 = y[lastIndex];
    a[lastIndex][lastIndex - 1] = 1.0 / (x1 - x0);
    a[lastIndex][lastIndex] = 2.0 / (x1 - x0);
    b[lastIndex] = 3.0 * (y1 - y0) / (std::pow((x1 - x0), 2.0));

    // solve for curvature vector k using gaussian elimination
    a[0][1] /= a[0][0];
    b[0] /= a[0][0];
    for (int i = 1; i < (knots - 1); i++)
    {
        a[i][i + 1] /= a[i][i] - a[i][i - 1] * a[i - 1][i];
        b[i] = (b[i] - a[i][i - 1] * b[i - 1]) / (a[i][i] - a[i][i - 1] * a[i - 1][i]);
    }
    int g1 = knots - 1;
    int g2 = knots - 2;
    b[g1] = (b[g1] - a[g1][g2] * b[g2]) / (a[g1][g1] - a[g1][g2] * a[g2][g1]);

    k[g1] = b[g1];
    for (int i = (knots - 2); i >= 0; i--)
    {
        k[i] = b[i] - a[i][i + 1] * k[i + 1];
    }

    // Fill the spline data
    for (unsigned int i = 0; i < m_splineList.size(); i++)
    {
        m_splineList[i].k1 = k[i];
        m_splineList[i].k2 = k[i + 1];
    }

    // deallocate
    delete[] b;
    delete[] k;
    for (int i = 0; i < knots; i++)
    {
        delete[] a[i];
    }
    delete[] a;
}

double CubicSpline::Evaluate(double x)
{
    if (!m_isValid)
    {
        std::cerr << "Error at CubicSpline::Evaluate! Unable to evaluate without first generating splines." << std::endl;
        return 0.0;
    }

    for (auto& spline : m_splineList)
    {
        if (x >= spline.x1 && x <= spline.x2)
            return spline.Evaluate(x);
    }

    // std::cerr<<"Error at CubicSpline::Evaluate! Input x value: "<<x<<" is not within the spline range min: "<<m_splineList[0].x1<<" max: "<<m_splineList[m_splineList.size()-1].x2<<std::endl;
    return 0.0;
}

// Purely for plotting in ROOT, do not use for caluculations.
// Note that for splines the parameter pointer p is unused
double CubicSpline::EvaluateROOT(double *x, double *p)
{
    if (!m_isValid)
    {
        std::cerr << "Error at CubicSpline::Evaluate! Unable to evaluate without first generating splines." << std::endl;
        return 0.0;
    }

    double xval = x[0]; //Ugh ROOT
    for (auto& spline : m_splineList)
    {
        if (xval >= spline.x1 && xval <= spline.x2)
            return spline.Evaluate(xval);
    }

    // std::cerr<<"Error at CubicSpline::EvaluateRoot! Input x value: "<<x<<" is not within the spline range min: "<<m_splineList[0].x1<<" max: "<<m_splineList[m_splineList.size()-1].x2<<std::endl;
    return 0.0;
}
