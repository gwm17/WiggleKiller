#include "Histogrammer.h"

void FillHistogram1D(std::unordered_map<std::string, std::shared_ptr<TObject>>& map, const Histo1DParameters& params, double value)
{
    auto iter = map.find(params.name.c_str());
    if(iter == map.end())
    {
        std::shared_ptr<TH1F> hRef = std::make_shared<TH1F>(params.name.c_str(), params.title.c_str(), params.nBins,
                                                            params.minVal, params.maxVal);
        hRef->Fill(value);
        map[params.name] = hRef;
    }
    else
    {
        std::static_pointer_cast<TH1F>(iter->second)->Fill(value);
    }
}

void FillHistogram2D(std::unordered_map<std::string, std::shared_ptr<TObject>>& map, const Histo2DParameters& params,
                     double valueX, double valueY)
{
    auto iter = map.find(params.name.c_str());
    if(iter == map.end())
    {
        std::shared_ptr<TH2F> hRef = std::make_shared<TH2F>(params.name.c_str(), params.title.c_str(), params.nBinsX,
                                                            params.minValX, params.maxValX, params.nBinsY, params.minValY, params.maxValY);
        hRef->Fill(valueX, valueY);
        map[params.name] = hRef;
    }
    else
    {
        std::static_pointer_cast<TH2F>(iter->second)->Fill(valueX, valueY);
    }
}