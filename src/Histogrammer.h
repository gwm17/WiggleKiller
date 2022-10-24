#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include <string>
#include <unordered_map>
#include <memory>
#include "TH1.h"
#include "TH2.h"


struct Histo1DParameters
{
    std::string name = "";
    std::string title = "";
    int nBins = 0;
    double minVal = 0.0;
    double maxVal = 0.0;
};

struct Histo2DParameters
{
    std::string name = "";
    std::string title = "";
    int nBinsX = 0;
    double minValX = 0.0;
    double maxValX = 0.0;
    int nBinsY = 0;
    double minValY = 0.0;
    double maxValY = 0.0;
};

void FillHistogram1D(std::unordered_map<std::string, std::shared_ptr<TObject>>& map, const Histo1DParameters& params, double value);

void FillHistogram2D(std::unordered_map<std::string, std::shared_ptr<TObject>>& map, const Histo2DParameters& params,
                     double valueX, double valueY);

#endif