/*
WiggleCorrector.cpp
Class for correcting focal plane spectra for differential-nonlinearity using an existing solution. Reads in calibration from file and generates the cubic spline
for interpolation. Appies to the data set. Correction can generates histograms as well as an optional tree of corrected data. Note: requires particle ID
(scint/cathode) and x1-x2 correltaion cuts. Resulting data has cuts applied.

Gordon M. May 2021
*/
#include "WiggleCorrector.h"
#include "Histogrammer.h"
#include <TROOT.h>
#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <cmath>

WiggleCorrector::WiggleCorrector() {}

WiggleCorrector::WiggleCorrector(const std::string& path, const std::string& tag, double lambdaFront, double lambdaBack)
{
    ReadCorrectionData(path, tag, lambdaFront, lambdaBack);
}

void WiggleCorrector::ReadCorrectionData(const std::string& path, const std::string& tag, double lambdaFront, double lambdaBack)
{
    std::string delayFrontLeftName = path + "delayFrontLeft_lambda" + std::to_string(lambdaFront) + "_" + tag + ".txt";
    std::string delayFrontRightName = path + "delayFrontRight_lambda" + std::to_string(lambdaFront) + "_" + tag + ".txt";
    std::string delayBackLeftName = path + "delayBackLeft_lambda" + std::to_string(lambdaBack) + "_" + tag + ".txt";
    std::string delayBackRightName = path + "delayBackRight_lambda" + std::to_string(lambdaBack) + "_" + tag + ".txt";
    frontLeftSpline.ReadFile(delayFrontLeftName);
    frontRightSpline.ReadFile(delayFrontRightName);
    backLeftSpline.ReadFile(delayBackLeftName);
    backRightSpline.ReadFile(delayBackRightName);

    if (!CheckSplines())
    {
        std::cerr << "Error at WiggleCorrector::ReadCorrectionData! Splines not properly initialized!" << std::endl;
    }
}

void WiggleCorrector::CopyBulkEvent(ProcessedEvent *inevent, ProcessedEvent &outevent)
{
    outevent = *inevent;
}

void WiggleCorrector::ApplyCorrection(const std::string& inputFileName, const std::string& outputFileName, bool writeTree)
{
    if (!CheckSplines() || !cuts.IsValid())
    {
        std::cerr << "Error at WiggleCorrector::ApplyCorrection! Cannot apply without splines!" << std::endl;
        return;
    }

    TFile *input = TFile::Open(inputFileName.c_str(), "READ");
    if (!input->IsOpen())
    {
        std::cerr << "Unable to open input data at WiggleCorrector::ApplyCorrection! Filename: " << inputFileName << std::endl;
        return;
    }
    TTree *intree = (TTree *)input->Get("SPSTree");

    ProcessedEvent *inevent = new ProcessedEvent();
    intree->SetBranchAddress("event", &inevent);

    XWeights xWeights = GetXWeights(m_fpParams);

    TFile *output = TFile::Open(outputFileName.c_str(), "RECREATE");
    if (!output->IsOpen())
    {
        std::cerr << "Unable to open output data at WiggleCorrector::ApplyCorrection! Filename: " << outputFileName << std::endl;
        return;
    }

    TTree* outtree = nullptr;
    if(writeTree)
        outtree = new TTree("SPSTree", "SPSTree");

    ProcessedEvent outevent;
    if(outtree)
        outtree->Branch("event", &outevent);

    std::unordered_map<std::string, std::shared_ptr<TObject>> histoMap;

    long blentries = intree->GetEntries();
    long count = 0, flush_count = 0, flush_val = 0.1 * blentries;
    double dfr_ti, dfl_ti, dbl_ti, dbr_ti;
    bool front, back;

    for (long i = 0; i < blentries; i++)
    {
        intree->GetEntry(i);
        count++;
        if (count == flush_val)
        {
            flush_count++;
            count = 0;
            std::cout << "\rPercent of file processed: " << flush_count * 10 << "%" << std::flush;
        }
        front = false;
        back = false;
        if (!cuts.IsInside(inevent->cathode, inevent->scintLeft, inevent->x1, inevent->x2, inevent->anodeFront, inevent->anodeBack))
            continue;

        CopyBulkEvent(inevent, outevent);
        if (inevent->delayFrontLeftE != -1 && inevent->delayFrontRightE != -1 && inevent->anodeFront != -1)
        {
            dfr_ti = inevent->delayFrontRightTime - inevent->anodeFrontTime;
            dfl_ti = inevent->delayFrontLeftTime - inevent->anodeFrontTime;
            outevent.delayFrontLeftTime = (dfl_ti + frontLeftSpline.Evaluate(dfl_ti)) + inevent->anodeFrontTime;
            outevent.delayFrontRightTime = (dfr_ti + frontRightSpline.Evaluate(dfr_ti)) + inevent->anodeFrontTime;
            outevent.fp1_tdiff = (outevent.delayFrontLeftTime - outevent.delayFrontRightTime) * 0.5;
            outevent.fp1_tsum = (outevent.delayFrontLeftTime + outevent.delayFrontRightTime) - 2.0 * inevent->scintLeftTime;
            outevent.fp1_tcheck = (outevent.delayFrontLeftTime + outevent.delayFrontRightTime) * 0.5 - inevent->anodeFrontTime;
            outevent.x1 = outevent.fp1_tdiff * 0.5050;
            FillHistogram1D(histoMap, {"x1_before", "x1", 600, -300, 300}, inevent->x1);
            FillHistogram1D(histoMap, {"x1_after", "x1", 600, -300, 300}, outevent.x1);
            front = true;
        }

        if (inevent->delayBackLeftE != -1 && inevent->delayBackRightE != -1 && inevent->anodeBack != -1)
        {
            dbr_ti = inevent->delayBackRightTime - inevent->anodeBackTime;
            dbl_ti = inevent->delayBackLeftTime - inevent->anodeBackTime;
            outevent.delayBackLeftTime = (dbl_ti + backLeftSpline.Evaluate(dbl_ti)) + inevent->anodeBackTime;
            outevent.delayBackRightTime = (dbr_ti + backRightSpline.Evaluate(dbr_ti)) + inevent->anodeBackTime;
            outevent.fp2_tdiff = (outevent.delayBackLeftTime - outevent.delayBackRightTime) * 0.5;
            outevent.fp2_tsum = (outevent.delayBackLeftTime + outevent.delayBackRightTime) - 2.0 * inevent->scintLeftTime;
            outevent.fp2_tcheck = (outevent.delayBackLeftTime + outevent.delayBackRightTime) * 0.5 - inevent->anodeBackTime;
            outevent.x2 = outevent.fp2_tdiff * 0.476;
            FillHistogram1D(histoMap, {"x2_before", "x2", 600, -300, 300}, inevent->x2);
            FillHistogram1D(histoMap, {"x2_after", "x2", 600, -300, 300}, outevent.x2);
            back = true;
        }

        if (front && back)
        {
            outevent.xavg = outevent.x1 * xWeights.x1Weight + outevent.x2 * xWeights.x2Weight;
            outevent.theta = std::atan((outevent.x2 - outevent.x1) / (Wire_Dist() * 10));
            FillHistogram1D(histoMap, {"xavg_before", "xavg", 600, -300, 300}, inevent->xavg);
            FillHistogram1D(histoMap, {"xavg_after", "xavg", 600, -300, 300}, outevent.xavg);
            FillHistogram2D(histoMap, {"xavg_theta_before", ";xavg;theta", 600, -300, 300, 180, 0, 90}, inevent->xavg,  inevent->theta * 180.0 / M_PI);
            FillHistogram2D(histoMap, {"xavg_theta_after", ";xavg;theta", 600, -300, 300, 180, 0, 90}, outevent.xavg, outevent.theta * 180.0 / M_PI);
        }

        if(outtree)
            outtree->Fill();
    }
    std::cout << std::endl;

    input->Close();
    output->cd();
    if(outtree)
        outtree->Write(outtree->GetName(), TObject::kOverwrite);
    for(auto& gram : histoMap)
        gram.second->Write(gram.second->GetName(), TObject::kOverwrite);
    output->Close();

    delete inevent;
}