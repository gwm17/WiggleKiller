/*
WiggleCorrector.cpp
Class for correcting focal plane spectra for differential-nonlinearity using an existing solution. Reads in calibration from file and generates the cubic spline
for interpolation. Appies to the data set. Correction can generates histograms as well as an optional tree of corrected data. Note: requires particle ID 
(scint/cathode) and x1-x2 correltaion cuts. Resulting data has cuts applied. 

Gordon M. May 2021
*/
#include "WiggleCorrector.h"
#include "FP_kinematics.h"
#include <TROOT.h>
#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <cmath>

WiggleCorrector::WiggleCorrector() {}

WiggleCorrector::WiggleCorrector(std::string& flname, std::string& frname, std::string& blname, std::string& brname) {
	ReadCorrectionData(flname, frname, blname, brname);
}

void WiggleCorrector::MyFill(THashTable* table, std::string name, std::string titlex, int binsx, double minx, double maxx, double valuex)
{
	TH1* h = (TH1*) table->FindObject(name.c_str());
	if(h == nullptr) {
		std::string title = name + ";" + titlex;
		TH1F* histo = new TH1F(name.c_str(), title.c_str(), binsx, minx, maxx);
		histo->Fill(valuex);
		table->Add(histo);
	} else {
		h->Fill(valuex);
	}
}

void WiggleCorrector::MyFill(THashTable* table, std::string name, std::string titlex, std::string titley, int binsx, double minx, double maxx, double valuex, int binsy, double miny, double maxy, double valuey)
{
	TH2* h = (TH2*) table->FindObject(name.c_str());
	if(h == nullptr) {
		std::string title = name + ";" + titlex + ";" + titley;
		TH2F* histo = new TH2F(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
		histo->Fill(valuex, valuey);
		table->Add(histo);
	} else {
		h->Fill(valuex, valuey);
	}
}

void WiggleCorrector::ReadCorrectionData(std::string& flname, std::string& frname, std::string& blname, std::string& brname) {
	frontLeftSpline.ReadFile(flname);
	frontRightSpline.ReadFile(frname);
	backLeftSpline.ReadFile(blname);
	backRightSpline.ReadFile(brname);

	if(!CheckSplines()) {
		std::cerr<<"Error at WiggleCorrector::ReadCorrectionData! Splines not properly initialized!"<<std::endl;
	}
}

std::pair<double, double> WiggleCorrector::GetXavgWeights(int zt, int at, int zp, int ap, int ze, int ae, double bke, double theta, double b) {
	double dz = Delta_Z(zt, at, zp, ap, ze, ae, bke, theta, b);
	double w1 = ((Wire_Dist()/2.0) - dz)/Wire_Dist();
	double w2 = 1.0 - w1;
	return std::make_pair(w1, w2);
}

void WiggleCorrector::CopyBulkEvent(ProcessedEvent* inevent, ProcessedEvent& outevent) {
	outevent.anodeFront = inevent->anodeFront;
	outevent.anodeBack = inevent->anodeBack;
	outevent.scintLeft = inevent->scintLeft;
	outevent.scintRight = inevent->scintRight;
	outevent.cathode = inevent->cathode;
	outevent.delayFrontRightE = inevent->delayFrontRightE;
	outevent.delayFrontLeftE = inevent->delayFrontLeftE;
	outevent.delayBackRightE = inevent->delayBackRightE;
	outevent.delayBackLeftE = inevent->delayBackLeftE;
	outevent.monitorE = inevent->monitorE;
	outevent.scintLeftShort = inevent->scintLeftShort;
	outevent.scintRightShort = inevent->scintRightShort;
	outevent.delayFrontRightShort = inevent->delayFrontRightShort;
	outevent.delayFrontLeftShort = inevent->delayFrontLeftShort;
	outevent.delayBackRightShort = inevent->delayBackRightShort;
	outevent.delayBackLeftShort = inevent->delayBackLeftShort;
	outevent.monitorShort = inevent->monitorShort;
	outevent.scintLeftTime = inevent->scintLeftTime;
	outevent.scintRightTime = inevent->scintRightTime;
	outevent.anodeFrontTime = inevent->anodeFrontTime;
	outevent.cathodeTime = inevent->cathodeTime;
	outevent.anodeBackTime = inevent->anodeBackTime;
	outevent.monitorTime = inevent->monitorTime;
	for(int i=0; i<5; i++) {
		outevent.sabreRingE[i] = inevent->sabreRingE[i];
		outevent.sabreWedgeE[i] = inevent->sabreWedgeE[i];
		outevent.sabreRingChannel[i] = inevent->sabreRingChannel[i];
		outevent.sabreWedgeChannel[i] = inevent->sabreWedgeChannel[i];
		outevent.sabreRingTime[i] = inevent->sabreRingTime[i];
		outevent.sabreWedgeTime[i] = inevent->sabreWedgeTime[i];
		outevent.sabreArray[i] = inevent->sabreArray[i];
	}
}

void WiggleCorrector::ApplyCorrection(std::string& inname, std::string& outname) {
	if(!CheckSplines() || !cuts.IsValid()) {
		std::cerr<<"Error at WiggleCorrector::ApplyCorrection! Cannot apply without splines!"<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inname.c_str(), "READ");
	if(!input->IsOpen()) {
		std::cerr<<"Unable to open input data at WiggleCorrector::ApplyCorrection! Filename: "<<inname<<std::endl;
		return;
	}
	TTree* intree = (TTree*) input->Get("SPSTree");
	ProcessedEvent* inevent = new ProcessedEvent();
	intree->SetBranchAddress("event", &inevent);
	TParameter<double>* at = (TParameter<double>*) input->Get("AT");
	TParameter<double>* zt = (TParameter<double>*) input->Get("ZT");
	TParameter<double>* ap = (TParameter<double>*) input->Get("AP");
	TParameter<double>* zp = (TParameter<double>*) input->Get("ZP");
	TParameter<double>* ae = (TParameter<double>*) input->Get("AE");
	TParameter<double>* ze = (TParameter<double>*) input->Get("ZE");
	TParameter<double>* b = (TParameter<double>*) input->Get("Bfield");
	TParameter<double>* bke = (TParameter<double>*) input->Get("BeamKE");
	TParameter<double>* theta = (TParameter<double>*) input->Get("Theta");
	auto weights = GetXavgWeights(zt->GetVal(), at->GetVal(), zp->GetVal(), ap->GetVal(), ze->GetVal(), ae->GetVal(), bke->GetVal(), theta->GetVal(), b->GetVal());

	TFile* output = TFile::Open(outname.c_str(), "RECREATE");
	if(!output->IsOpen()) {
		std::cerr<<"Unable to open output data at WiggleCorrector::ApplyCorrection! Filename: "<<outname<<std::endl;
		return;
	}
	TTree* outtree = new TTree("SPSTree", "SPSTree");
	ProcessedEvent outevent;
	ProcessedEvent blank;
	outtree->Branch("event", &outevent);
	THashTable* table = new THashTable();

	long blentries = intree->GetEntries();
	long count=0, flush_count=0, flush_val = 0.1*blentries;
	double dfr_ti, dfl_ti, dbl_ti, dbr_ti;
	bool front, back;

	for(long i=0; i<blentries; i++) {
		intree->GetEntry(i);
		count++;
		outevent = blank;
		if(count == flush_val) {
			flush_count++;
			count = 0;
			std::cout<<"\rPercent of file processed: "<<flush_count*10<<"%"<<std::flush;
		}
		front = false;
		back = false;
		if(!cuts.IsInside(inevent->cathode, inevent->scintLeft, inevent->x1, inevent->x2, inevent->anodeFront, inevent->anodeBack)) continue;

		CopyBulkEvent(inevent, outevent);
		if(inevent->delayFrontLeftE != -1 && inevent->delayFrontRightE != -1 && inevent->anodeFront != -1) {
			dfr_ti = inevent->delayFrontRightTime - inevent->anodeFrontTime;
			dfl_ti = inevent->delayFrontLeftTime - inevent->anodeFrontTime;
			outevent.delayFrontLeftTime = (dfl_ti + frontLeftSpline.Evaluate(dfl_ti))+inevent->anodeFrontTime;
			outevent.delayFrontRightTime = (dfr_ti + frontRightSpline.Evaluate(dfr_ti))+inevent->anodeFrontTime;
			outevent.fp1_tdiff = (outevent.delayFrontLeftTime - outevent.delayFrontRightTime)*0.5;
			outevent.fp1_tsum = (outevent.delayFrontLeftTime + outevent.delayFrontRightTime) - 2.0*inevent->scintLeftTime;
			outevent.fp1_tcheck = (outevent.delayFrontLeftTime + outevent.delayFrontRightTime)*0.5 - inevent->anodeFrontTime;
			outevent.x1 = outevent.fp1_tdiff*0.5050;
			MyFill(table,"x1_before","x1",600,-300,300,inevent->x1); 
			MyFill(table,"x1_after","x1",600,-300,300,outevent.x1); 
			front = true;
		}

		if(inevent->delayBackLeftE != -1 && inevent->delayBackRightE != -1 && inevent->anodeBack != -1) {
			dbr_ti = inevent->delayBackRightTime - inevent->anodeBackTime;
			dbl_ti = inevent->delayBackLeftTime - inevent->anodeBackTime;
			outevent.delayBackLeftTime = (dbl_ti + backLeftSpline.Evaluate(dbl_ti))+inevent->anodeBackTime;
			outevent.delayBackRightTime = (dbr_ti + backRightSpline.Evaluate(dbr_ti))+inevent->anodeBackTime;
			outevent.fp2_tdiff = (outevent.delayBackLeftTime - outevent.delayBackRightTime)*0.5;
			outevent.fp2_tsum = (outevent.delayBackLeftTime + outevent.delayBackRightTime) - 2.0*inevent->scintLeftTime;
			outevent.fp2_tcheck = (outevent.delayBackLeftTime + outevent.delayBackRightTime)*0.5 - inevent->anodeBackTime;
			outevent.x2 = outevent.fp2_tdiff*0.476;
			MyFill(table,"x2_before","x2",600,-300,300,inevent->x2); 
			MyFill(table,"x2_after","x2",600,-300,300,outevent.x2); 
			back = true;
		}

		if(front && back) {
			outevent.xavg = outevent.x1*weights.first + outevent.x2*weights.second;
			outevent.theta = std::atan((outevent.x2 - outevent.x1)/(Wire_Dist()*10));
			MyFill(table,"xavg_before","xavg",600,-300,300,inevent->xavg); 
			MyFill(table,"xavg_after","xavg",600,-300,300,outevent.xavg);
			MyFill(table,"xavg_theta_before","xavg","theta",600,-300,300,inevent->xavg,180,0,90,inevent->theta*180.0/M_PI);
			MyFill(table,"xavg_theta_after","xavg","theta",600,-300,300,outevent.xavg,180,0,90,outevent.theta*180.0/M_PI);
		}
		outtree->Fill();
	}
	std::cout<<std::endl;

	input->Close();
	output->cd();
	outtree->Write(outtree->GetName(), TObject::kOverwrite);
	table->Write();
	output->Close();

	delete inevent;
}

void WiggleCorrector::ApplyCorrection_NoTree(std::string& inname, std::string& outname) {
	if(!CheckSplines() || !cuts.IsValid()) {
		std::cerr<<"Error at WiggleCorrector::ApplyCorrection! Cannot apply without splines!"<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inname.c_str(), "READ");
	if(!input->IsOpen()) {
		std::cerr<<"Unable to open input data at WiggleCorrector::ApplyCorrection! Filename: "<<inname<<std::endl;
		return;
	}
	TTree* intree = (TTree*) input->Get("SPSTree");
	ProcessedEvent* inevent = new ProcessedEvent();
	intree->SetBranchAddress("event", &inevent);
	TParameter<double>* at = (TParameter<double>*) input->Get("AT");
	TParameter<double>* zt = (TParameter<double>*) input->Get("ZT");
	TParameter<double>* ap = (TParameter<double>*) input->Get("AP");
	TParameter<double>* zp = (TParameter<double>*) input->Get("ZP");
	TParameter<double>* ae = (TParameter<double>*) input->Get("AE");
	TParameter<double>* ze = (TParameter<double>*) input->Get("ZE");
	TParameter<double>* b = (TParameter<double>*) input->Get("Bfield");
	TParameter<double>* bke = (TParameter<double>*) input->Get("BeamKE");
	TParameter<double>* theta = (TParameter<double>*) input->Get("Theta");
	auto weights = GetXavgWeights(zt->GetVal(), at->GetVal(), zp->GetVal(), ap->GetVal(), ze->GetVal(), ae->GetVal(), bke->GetVal(), theta->GetVal(), b->GetVal());

	TFile* output = TFile::Open(outname.c_str(), "RECREATE");
	if(!output->IsOpen()) {
		std::cerr<<"Unable to open output data at WiggleCorrector::ApplyCorrection! Filename: "<<outname<<std::endl;
		return;
	}
	THashTable* table = new THashTable();

	long blentries = intree->GetEntries();
	long count=0, flush_count=0, flush_val = 0.1*blentries;
	double dfr_ti, dfl_ti, dbl_ti, dbr_ti;
	double dfr_tc, dfl_tc, dbl_tc, dbr_tc;
	double x1c, x2c, xac, thetac;
	bool front, back;

	for(long i=0; i<blentries; i++) {
		intree->GetEntry(i);
		count++;
		if(count == flush_val) {
			flush_count++;
			count = 0;
			std::cout<<"\rPercent of file processed: "<<flush_count*10<<"%"<<std::flush;
		}
		front = false;
		back = false;
		if(!cuts.IsInside(inevent->cathode, inevent->scintLeft, inevent->x1, inevent->x2, inevent->anodeFront, inevent->anodeBack)) continue;

		if(inevent->delayFrontLeftE != -1 && inevent->delayFrontRightE != -1 && inevent->anodeFront != -1) {
			dfr_ti = inevent->delayFrontRightTime - inevent->anodeFrontTime;
			dfl_ti = inevent->delayFrontLeftTime - inevent->anodeFrontTime;
			dfl_tc= (dfl_ti + frontLeftSpline.Evaluate(dfl_ti))+inevent->anodeFrontTime;
			dfr_tc = (dfr_ti + frontRightSpline.Evaluate(dfr_ti))+inevent->anodeFrontTime;
			x1c = (dfl_tc - dfr_tc)*0.5*0.5050;
			MyFill(table,"x1_before","x1",600,-300,300,inevent->x1); 
			MyFill(table,"x1_after","x1",600,-300,300,x1c); 
			front = true;
		}

		if(inevent->delayBackLeftE != -1 && inevent->delayBackRightE != -1 && inevent->anodeBack != -1) {
			dbr_ti = inevent->delayBackRightTime - inevent->anodeBackTime;
			dbl_ti = inevent->delayBackLeftTime - inevent->anodeBackTime;
			dbl_tc = (dbl_ti + backLeftSpline.Evaluate(dbl_ti))+inevent->anodeBackTime;
			dbr_tc = (dbr_ti + backRightSpline.Evaluate(dbr_ti))+inevent->anodeBackTime;
			x2c = (dbl_tc - dbr_tc)*0.5*0.476;
			MyFill(table,"x2_before","x2",600,-300,300,inevent->x2); 
			MyFill(table,"x2_after","x2",600,-300,300,x2c); 
			back = true;
		}

		if(front && back) {
			xac = x1c*weights.first + x2c*weights.second;
			thetac = std::atan((x2c - x1c)/(Wire_Dist()*10));
			MyFill(table,"xavg_before","xavg",600,-300,300,inevent->xavg); 
			MyFill(table,"xavg_after","xavg",600,-300,300,xac);
			MyFill(table,"xavg_theta_before","xavg","theta",600,-300,300,inevent->xavg,180,0,90,inevent->theta*180.0/M_PI);
			MyFill(table,"xavg_theta_after","xavg","theta",600,-300,300,xac,180,0,90,thetac*180.0/M_PI);
		}
	}
	std::cout<<std::endl;

	input->Close();
	output->cd();
	table->Write();
	output->Close();

	delete inevent;
}