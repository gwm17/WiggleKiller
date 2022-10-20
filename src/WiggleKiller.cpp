/*
WiggleKiller.cpp
Class for correcting focal plane spectra for differential-nonlinearity using a new solution or optimizing a solution. Can be used to simply test a solution
(guess values of lambda) or optimize a guess using a simluated annealing method. Can generate calibration files for use with the WiggleCorrector class. Note:
requires particle ID (scint/cathode) and x1-x2 correlation cuts.

Gordon M. May 2021
*/
#include "WiggleKiller.h"
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TParameter.h>
#include "FP_kinematics.h"
#include <fstream>
#include <cmath>

WiggleKiller::WiggleKiller() {
	nbins = BINMAX/BINWIDTH;
	delayBinnedAvgPSD.resize(4);
	delayBinnedCounts.resize(4);
	PSDsplines.resize(4);
	calibData.resize(4);
	for(int i=0; i<4; i++) {
		delayBinnedAvgPSD[i].resize(nbins);
		delayBinnedCounts[i].resize(nbins);
		calibData[i].resize(nbins);
		for(int j=0; j<nbins; j++) {
			calibData[i][j].time_i = j*BINWIDTH;
		}
	}
	generator = new TRandom3(0);
}

WiggleKiller::WiggleKiller(std::string& edename, std::string& xxname, double lf, double lb) :
	cuts(edename, xxname)
{
	SetLambdas(lf, lb);
	nbins = BINMAX/BINWIDTH;
	delayBinnedAvgPSD.resize(4);
	delayBinnedCounts.resize(4);
	PSDsplines.resize(4);
	calibData.resize(4);
	for(int i=0; i<4; i++) {
		delayBinnedAvgPSD[i].resize(nbins);
		delayBinnedCounts[i].resize(nbins);
		calibData[i].resize(nbins);
		for(int j=0; j<nbins; j++) {
			calibData[i][j].time_i = j*BINWIDTH;
		}
	}
	generator = new TRandom3(0);
}

WiggleKiller::~WiggleKiller() {
	delete generator;
}

void WiggleKiller::MyFill(THashTable* table, std::string name, std::string titlex, int binsx, double minx, double maxx, double valuex)
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

void WiggleKiller::MyFill(THashTable* table, std::string name, std::string titlex, std::string titley, int binsx, double minx, double maxx, double valuex, int binsy, double miny, double maxy, double valuey)
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

void WiggleKiller::MyClear(THashTable* table, std::string name) {
	TH1* h = (TH1*) table->FindObject(name.c_str());
	if(h) {
		h->Reset();
	} else {
		TH2* h2 = (TH2*) table->FindObject(name.c_str());
		if(h2) h2->Reset();
	}
}

std::pair<double, double> WiggleKiller::GetXavgWeights(int zt, int at, int zp, int ap, int ze, int ae, double bke, double theta, double b) {
	double dz = Delta_Z(zt, at, zp, ap, ze, ae, bke, theta, b);
	double w1 = ((Wire_Dist()/2.0) - dz)/Wire_Dist();
	double w2 = 1.0 - w1;
	return std::make_pair(w1, w2);
}

void WiggleKiller::SetLambdas(double lf, double lb) {
	lambdaFront = lf;
	lambdaBack = lb;
}

void WiggleKiller::SetCuts(std::string& edename, std::string& xxname) {
	cuts.SetCuts(edename, xxname);
}

/*
Core concept of the correction. Wiggles (differential non-linearity) is a result of varying pulse shape causing distortions in the PSD analysis
of the delay line. To correct we apply:

t' = t + lambda*psd

lambda is a free parameter deterimined/optimized. In application of a calibration, the assumption is that lambda*psd is consistent between data sets. In other
words, the global offset of time is irrelevant, only the slight deviations between left and right matter.
*/
double WiggleKiller::ApplyPSDCorrection(double time, double psd, int delay, double lambda) {
	int this_bin = floor(time/BINWIDTH);
	if(this_bin <0 || this_bin > nbins) {
		return time;
	}
	
	double tprime = time + lambda*(PSDsplines[delay].Evaluate(time));
	return tprime;
}

//Old optimization parameter
double WiggleKiller::GetTCheck1Width(double lambda) {
	double dfr_tc, dfl_tc;
	double tcheck1c;

	std::vector<double> vals;
	long double mean = 0;

	for(auto& point : frontPoints) {
		dfr_tc = ApplyPSDCorrection(point.rightTime, point.rightPSD, FrontRight, lambda);
		dfl_tc = ApplyPSDCorrection(point.leftTime, point.leftPSD, FrontLeft, lambda);
		tcheck1c = (dfr_tc + dfl_tc)/2.0;
		vals.push_back(tcheck1c);
		mean += tcheck1c;
	}
	mean /= vals.size();
	double variance = 0;
	for(auto& value : vals) {
		variance += std::pow(value-mean, 2.0);
	}

	double sigma1 = std::sqrt(variance/vals.size());
	return sigma1;
}

//Old optimization parameter
double WiggleKiller::GetTCheck2Width(double lambda) {
	double dbr_tc, dbl_tc;
	double tcheck2c;

	std::vector<double> vals;
	double mean = 0;

	for(auto& point : backPoints) {
		dbr_tc = ApplyPSDCorrection(point.rightTime, point.rightPSD, BackRight, lambda);
		dbl_tc = ApplyPSDCorrection(point.leftTime, point.leftPSD, BackLeft, lambda);
		tcheck2c = (dbr_tc + dbl_tc)/2.0;
		vals.push_back(tcheck2c);
		mean += tcheck2c;
	}
	mean /= vals.size();
	double variance = 0;
	for(auto& value : vals) {
		variance += std::pow(value-mean, 2.0);
	}

	double sigma2 = std::sqrt(variance/vals.size());
	return sigma2;
}

/*
Current optimization is to reduce peak-to-trough distance of wiggles in a region of clear DNL. Squicky as it uses the spectrum of interest as well
as a strong bias towards intrinsic smoothness. But it works far better than any other metric. 
*/
double WiggleKiller::GetX1Peak2Trough(double lambda) {
	double peak_location = -127.0;
	double trough_location = -117.0;

	double dfr_tc, dfl_tc, x1c;
	double bin_min = -150;
	double bin_max = -110;
	double bin_width = 1.0;
	int nbins_local =  ceil((bin_max - bin_min)/bin_width);
	int peak_bin = floor((peak_location-bin_min)/bin_width);
	int trough_bin = floor((trough_location-bin_min)/bin_width);

	std::vector<double> counts;
	counts.resize(nbins_local);
	int this_bin;
	for(auto& point : frontPoints) {
		dfr_tc = ApplyPSDCorrection(point.rightTime, point.rightPSD, FrontRight, lambda);
		dfl_tc = ApplyPSDCorrection(point.leftTime, point.leftPSD, FrontLeft, lambda);
		x1c = (dfl_tc - dfr_tc)*0.5*0.5050;
		this_bin = floor((x1c - bin_min)/bin_width);
		if(this_bin < 0 || this_bin > nbins_local-1) continue;
		counts[this_bin]++;
	}
	double p2t = std::fabs(counts[peak_bin]-counts[trough_bin]);
	return p2t;

}

double WiggleKiller::GetX2Peak2Trough(double lambda) {
	double peak_location = -105.0;
	double trough_location = -90.0;

	double dbr_tc, dbl_tc, x2c;
	double bin_min = -120;
	double bin_max = -70;
	double bin_width = 1.0;
	int nbins_local =  ceil((bin_max - bin_min)/bin_width);
	int peak_bin = floor((peak_location-bin_min)/bin_width);
	int trough_bin = floor((trough_location-bin_min)/bin_width);

	std::vector<double> counts;
	counts.resize(nbins_local);
	int this_bin;
	for(auto& point : backPoints) {
		dbr_tc = ApplyPSDCorrection(point.rightTime, point.rightPSD, FrontRight, lambda);
		dbl_tc = ApplyPSDCorrection(point.leftTime, point.leftPSD, FrontLeft, lambda);
		x2c = (dbl_tc - dbr_tc)*0.5*0.5050;
		this_bin = floor((x2c - bin_min)/bin_width);
		if(this_bin < 0 || this_bin > nbins_local-1) continue;
		counts[this_bin]++;
	}
	double p2t = std::fabs(counts[peak_bin]-counts[trough_bin]);
	return p2t;

}

//Experimetnal optimization parameter; not to be used
double WiggleKiller::GetThetaWidth(double lambda1, double lambda2) {
	double dfr_tc, dfl_tc, dbr_tc, dbl_tc;
	double x1c, x2c;
	double thetac;

	std::vector<double> vals;
	double mean = 0;

	for(unsigned int i=0; i<frontPoints.size(); i++) {
		CorrectionPoint& front = frontPoints[i];
		CorrectionPoint& back = backPoints[i];

		dfr_tc = ApplyPSDCorrection(front.rightTime, front.rightPSD, FrontRight, lambda1);
		dfl_tc = ApplyPSDCorrection(front.leftTime, front.leftPSD, FrontLeft, lambda1);
		dbr_tc = ApplyPSDCorrection(back.rightTime, back.rightPSD, BackRight, lambda1);
		dbl_tc = ApplyPSDCorrection(back.leftTime, back.leftPSD, BackLeft, lambda1);

		x1c = (dfl_tc - dfr_tc)*0.5*0.505;
		x2c = (dbl_tc - dbr_tc)*0.5*0.476;
		thetac = std::atan((x2c - x1c)/(Wire_Dist()*10.0));
		vals.push_back(thetac);
		mean += thetac;
	}
	mean /= vals.size();
	double variance = 0;
	for(auto& value : vals) {
		variance += std::pow(value-mean, 2.0);
	}
	double sigma = std::sqrt(variance/vals.size());
	return sigma;
}

std::pair<double,double> WiggleKiller::MinimizeSimAnneal() {
	std::cout<<"Simulated Annealing method... reliability over speed"<<std::endl;
	int time;
	double maxTemp = 100000.0;
	double minTemp = 0.1;
	double curTemp;
	double tau = 1e4;

	
	double cur_lambda1, cur_lambda2;
	double new_lambda1, new_lambda2;
	double cur_sigma1, cur_sigma2;
	double new_sigma1, new_sigma2;
	
	//double cur_sigma, new_sigma;
	cur_lambda1 = lambdaFront;
	cur_sigma1 = GetX1Peak2Trough(cur_lambda1);
	cur_lambda2 = lambdaBack;
	cur_sigma2 = GetX2Peak2Trough(cur_lambda2);

	time=0;
	int counter = 0;
	int flush_val = 1000;
	int flush_count = 0;
	curTemp = maxTemp;
	std::cout<<"Simultaneous minimiziation starting... Initial widths -> Sigma1: "<<cur_sigma1<<" Sigma2: "<<cur_sigma2<<std::endl;
	while(curTemp > minTemp) {
		time++;
		counter++;
		if(counter == flush_val) {
			counter = 0;
			flush_count++;
			std::cout<<"\rNumber of time steps: "<<flush_count*1000<<" Current temperature: "<<curTemp<<std::flush;
		}

		curTemp = maxTemp*std::exp(-time/tau); //slowly cool the system

		
		new_lambda1 = cur_lambda1 + generator->Gaus(); //sample moves from gaussian
		new_lambda2 = cur_lambda2 + generator->Gaus();
		new_sigma1 = GetX1Peak2Trough(new_lambda1);
		new_sigma2 = GetX2Peak2Trough(new_lambda2);
		
		if(generator->Rndm() < std::exp(-1.0*(new_sigma1 - cur_sigma1)/curTemp)) { //Metropolis method
			cur_lambda1 = new_lambda1;
			cur_sigma1 = new_sigma1;
		}

		if(generator->Rndm() < std::exp(-1.0*(new_sigma2 - cur_sigma2)/curTemp)) {
			cur_lambda2 = new_lambda2;
			cur_sigma2 = new_sigma2;
		}
		
		widths1.push_back(cur_sigma1);
		widths2.push_back(cur_sigma2);
		times.push_back(time);
	}
	std::cout<<std::endl;

	return std::make_pair(cur_lambda1, cur_lambda2);
}

void WiggleKiller::OptimizeParameters(std::string& inname, std::string& outname) {

	if(!cuts.IsValid()) {
		std::cerr<<"Invalid cuts! Unable to run."<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inname.c_str(), "READ");
	TTree* intree = (TTree*) input->Get("SPSTree");

	ProcessedEvent* event_add = new ProcessedEvent();

	intree->SetBranchAddress("event", &event_add);

	TParameter<double>* at = (TParameter<double>*) input->Get("AT");
	TParameter<double>* zt = (TParameter<double>*) input->Get("ZT");
	TParameter<double>* ap = (TParameter<double>*) input->Get("AP");
	TParameter<double>* zp = (TParameter<double>*) input->Get("ZP");
	TParameter<double>* ae = (TParameter<double>*) input->Get("AE");
	TParameter<double>* ze = (TParameter<double>*) input->Get("ZE");
	TParameter<double>* b = (TParameter<double>*) input->Get("Bfield");
	TParameter<double>* bke = (TParameter<double>*) input->Get("BeamKE");
	TParameter<double>* theta = (TParameter<double>*) input->Get("Theta");
	weights = GetXavgWeights(zt->GetVal(), at->GetVal(), zp->GetVal(), ap->GetVal(), ze->GetVal(), ae->GetVal(), bke->GetVal(), theta->GetVal(), b->GetVal());

	TFile* output = TFile::Open(outname.c_str(), "RECREATE");

	THashTable* table = new THashTable();

	long blentries = intree->GetEntries();
	long flush = 0.1*blentries;
	long count=0, flush_count=0;

	double dfr_ti, dfl_ti, dbr_ti, dbl_ti;
	double dfr_psd, dfl_psd, dbr_psd, dbl_psd;
	double dfr_tc, dfl_tc, dbr_tc, dbl_tc;
	double x1c, x2c, xac, thetac, thetai;
	double tcheck1, tcheck2, tcheck1c, tcheck2c;
	int this_bin;

	CorrectionPoint this_point;
	for(int i=0; i<4; i++) {
		for(int j=0; j<nbins; j++) {
			delayBinnedAvgPSD[i][j] = 0.0;
			delayBinnedCounts[i][j] = 0.0;
			calibData[i][j].deltas.clear();
		}
	}

	std::vector<double> delaytimes;
	for(int i=0; i<nbins; i++) {
		delaytimes.push_back(i*BINWIDTH);
	}


	std::cout<<"Generating PSD correction data..."<<std::endl;

	for(long i=0; i<blentries; i++) {
		intree->GetEntry(i);
		count++;
		if(count == flush) {
			count =0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*10<<"%"<<std::flush;
		}

		if(event_add->xavg == -1e6) continue;

		if(!cuts.IsInside(event_add->cathode, event_add->scintLeft, event_add->x1, event_add->x2, event_add->anodeFrontTime, event_add->anodeBackTime)) continue;



		dfr_ti = event_add->delayFrontRightTime - event_add->anodeFrontTime;
		dfl_ti = event_add->delayFrontLeftTime - event_add->anodeFrontTime;
		dbr_ti = event_add->delayBackRightTime - event_add->anodeBackTime;
		dbl_ti = event_add->delayBackLeftTime - event_add->anodeBackTime;

		if(dfr_ti <0 || dfl_ti<0 || dbr_ti<0 || dbl_ti<0) {
			continue;
		}

		tcheck1 = (dfr_ti + dfl_ti)/2.0;
		tcheck2 = (dbr_ti + dbl_ti)/2.0;

		dfr_psd = (1.0 - event_add->delayFrontRightShort/event_add->delayFrontRightE);
		dfl_psd = (1.0 - event_add->delayFrontLeftShort/event_add->delayFrontLeftE);
		dbr_psd = (1.0 - event_add->delayBackRightShort/event_add->delayBackRightE);
		dbl_psd = (1.0 - event_add->delayBackLeftShort/event_add->delayBackLeftE);

		thetai = std::atan((event_add->x2 - event_add->x1)/(Wire_Dist()*10));

		MyFill(table,"x1psdl","x1","PSD Left",600,-300,300, event_add->x1,1000,0,1.0, dfl_psd);
		MyFill(table,"x1psdr","x1","PSD Right",600,-300,300, event_add->x1,1000,0,1.0, dfr_psd);
		MyFill(table,"x2psdl","x2","PSD Left",600,-300,300, event_add->x2,1000,0,1.0, dbl_psd);
		MyFill(table,"x2psdr","x2","PSD Right",600,-300,300, event_add->x2,1000,0,1.0, dbr_psd);
		MyFill(table,"x1","x1",600,-300,300,event_add->x1);
		MyFill(table,"x2","x2",600,-300,300,event_add->x2);
		MyFill(table,"xavg","xavg",600,-300,300,event_add->xavg);
		MyFill(table,"dfr_ti","delayFR Time",2000,-500,1500,dfr_ti);
		MyFill(table,"dfl_ti","delayFL Time",2000,-500,1500,dfl_ti);
		MyFill(table,"dbr_ti","delayBR Time",2000,-500,1500,dbr_ti);
		MyFill(table,"dbl_ti","delayBL Time",2000,-500,1500,dbl_ti);
		MyFill(table,"tcheck1","tcheck1",500,600,650,tcheck1);
		MyFill(table,"tcheck2","tcheck2",500,630,680,tcheck2);
		MyFill(table,"tc1_x1","x1","tcheck1",600,-300,300,event_add->x1,500,600,650, tcheck1);
		MyFill(table,"tc1_dfr","delayFR Time","tcheck1",1500,0,1500,dfr_ti,500,600,650, tcheck1);
		MyFill(table,"tc1_dfl","delayFL Time","tcheck1",1500,0,1500,dfl_ti,500,600,650, tcheck1);
		MyFill(table,"tc2_dbr","delayBR Time","tcheck2",1500,0,1500,dbr_ti,500,630,680, tcheck2);
		MyFill(table,"tc2_dbl","delayBL Time","tcheck2",1500,0,1500,dbl_ti,500,630,680, tcheck2);
		MyFill(table,"tc2_x2","x2","tcheck2",600,-300,300,event_add->x2,500,630,680, tcheck2);
		MyFill(table,"dfr_ti_psd","delayFR Time","PSD",2000,-500,1500,dfr_ti,1000,0,1.0,dfr_psd);
		MyFill(table,"dfl_ti_psd","delayFL Time","PSD",2000,-500,1500,dfl_ti,1000,0,1.0,dfl_psd);
		MyFill(table,"dbr_ti_psd","delayBR Time","PSD",2000,-500,1500,dbr_ti,1000,0,1.0,dbr_psd);
		MyFill(table,"dbl_ti_psd","delayBL Time","PSD",2000,-500,1500,dbl_ti,1000,0,1.0,dbl_psd);
		MyFill(table,"xavg_theta","xavg","theta",600,-300,300,event_add->xavg,180,0,90, thetai*180.0/M_PI);

		this_bin = floor(dfr_ti/BINWIDTH);
		delayBinnedAvgPSD[FrontRight][this_bin] += dfr_psd;
		delayBinnedCounts[FrontRight][this_bin]++;


		this_bin = floor(dfl_ti/BINWIDTH);
		delayBinnedAvgPSD[FrontLeft][this_bin] += dfl_psd;
		delayBinnedCounts[FrontLeft][this_bin]++;
		
		this_bin = floor(dbr_ti/BINWIDTH);
		delayBinnedAvgPSD[BackRight][this_bin] += dbr_psd;
		delayBinnedCounts[BackRight][this_bin]++;

		this_bin = floor(dbl_ti/BINWIDTH);
		delayBinnedAvgPSD[BackLeft][this_bin] += dbl_psd;
		delayBinnedCounts[BackLeft][this_bin]++;

		//if((dfl_ti>470.0 && dfl_ti<550.0)) {
		//if(generator->Rndm() <0.1) {
		//if(event_add->x1 > 100.0 && event_add->x1 < 160.0) {
		if(event_add->x1 > -150.0 && event_add->x1 < -110.0) {
		//if(event_add->xavg > 80.0 && event_add->xavg < 120.0) {
		//if(event_add->xavg > -200.0 && event_add->xavg < -160.0) {
			this_point.leftTime = dfl_ti;
			this_point.rightTime = dfr_ti;
			this_point.leftPSD = dfl_psd;
			this_point.rightPSD = dfr_psd;
			frontPoints.push_back(this_point);
			
		}
		//if((dbl_ti>370.0 && dbl_ti<450.0)) {
		//if((event_add->x2>150.0 && event_add->x2<210.0)) {
		if((event_add->x2>-120.0 && event_add->x2<-70.0)) {
			this_point.leftTime = dbl_ti;
			this_point.rightTime = dbr_ti;
			this_point.leftPSD = dbl_psd;
			this_point.rightPSD = dbr_psd;
			backPoints.push_back(this_point);
			
		}
	}
	std::cout<<std::endl;

	for(int i=0; i<4; i++) {
		for(int j=0; j<nbins; j++) {
			if(delayBinnedCounts[i][j] == 0) continue;
			delayBinnedAvgPSD[i][j] /= delayBinnedCounts[i][j];
		}
	}

	PSDsplines[FrontLeft].ReadData(delaytimes, delayBinnedAvgPSD[FrontLeft]);
	PSDsplines[FrontRight].ReadData(delaytimes, delayBinnedAvgPSD[FrontRight]);
	PSDsplines[BackLeft].ReadData(delaytimes, delayBinnedAvgPSD[BackLeft]);
	PSDsplines[BackRight].ReadData(delaytimes, delayBinnedAvgPSD[BackRight]);
	
	std::cout<<"Sample size front: "<<frontPoints.size()<<" back: "<<backPoints.size()<<std::endl;
	std::cout<<"Finished. Applying PSD correction..."<<std::endl;
	std::cout<<"Minimizing tcheck to find optimal PSD parameters... this make take some time so put on your big person pants"<<std::endl;
	auto min_vals = MinimizeSimAnneal();
	lambdaFront = min_vals.first;
	lambdaBack = min_vals.second;
	std::cout<<"Finished. Lambda1: "<<lambdaFront<<" Lambda2: "<<lambdaBack<<std::endl;
	std::cout<<"Generating final histograms..."<<std::endl;
	
	std::cout<<std::endl;
	flush_count=0;
	count=0;
	MyClear(table, "tcheck1_cor");
	MyClear(table, "tcheck2_cor");
	for(long i=0; i<blentries; i++) {
		intree->GetEntry(i);
		count++;
		if(count == flush) {
			count =0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*10<<"%"<<std::flush;
		}

		if(event_add->xavg == -1e6) continue;
		if(!cuts.IsInside(event_add->cathode, event_add->scintLeft, event_add->x1, event_add->x2, event_add->anodeFrontTime, event_add->anodeBackTime)) continue;

		
		dfr_ti = event_add->delayFrontRightTime - event_add->anodeFrontTime;
		dfl_ti = event_add->delayFrontLeftTime - event_add->anodeFrontTime;
		dbr_ti = event_add->delayBackRightTime - event_add->anodeBackTime;
		dbl_ti = event_add->delayBackLeftTime - event_add->anodeBackTime;

		dfr_psd = (1.0 - event_add->delayFrontRightShort/event_add->delayFrontRightE);
		dfl_psd = (1.0 - event_add->delayFrontLeftShort/event_add->delayFrontLeftE);
		dbr_psd = (1.0 - event_add->delayBackRightShort/event_add->delayBackRightE);
		dbl_psd = (1.0 - event_add->delayBackLeftShort/event_add->delayBackLeftE);

		dfr_tc = ApplyPSDCorrection(dfr_ti, dfr_psd, FrontRight, min_vals.first);
		dfl_tc = ApplyPSDCorrection(dfl_ti, dfl_psd, FrontLeft, min_vals.first);
		dbr_tc = ApplyPSDCorrection(dbr_ti, dbr_psd, BackRight, min_vals.second);
		dbl_tc = ApplyPSDCorrection(dbl_ti, dbl_psd, BackLeft, min_vals.second);

		this_bin = floor(dfr_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[FrontRight][this_bin].deltas.push_back(dfr_tc - dfr_ti);

		this_bin = floor(dfl_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[FrontLeft][this_bin].deltas.push_back(dfl_tc - dfl_ti);

		this_bin = floor(dbl_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[BackLeft][this_bin].deltas.push_back(dbl_tc - dbl_ti);

		this_bin = floor(dbr_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[BackRight][this_bin].deltas.push_back(dbr_tc - dbr_ti);

		x1c = (dfl_tc - dfr_tc)*0.5*0.505;
		x2c = (dbl_tc - dbr_tc)*0.5*0.476;
		xac = x1c*weights.first + x2c*weights.second;
		thetac = std::atan((x2c - x1c)/(Wire_Dist()*10.0));
		
		tcheck1c = (dfr_tc + dfl_tc)/2.0;
		tcheck2c = (dbr_tc + dbl_tc)/2.0;


		MyFill(table,"x1psdl_cor","x1","PSD Left",600,-300,300,x1c,1000,0.0,1.0, dfl_psd);
		MyFill(table,"x1psdr_cor","x1","PSD Right",600,-300,300,x1c,1000,0.0,1.0, dfr_psd);
		MyFill(table,"x2psdl_cor","x2","PSD Left",600,-300,300,x2c,1000,0.0,1.0, dbl_psd);
		MyFill(table,"x2psdr_cor","x2","PSD Right",600,-300,300,x1c,1000,0.0,1.0, dbr_psd);
		MyFill(table,"x1_cor","x1",600,-300,300,x1c);
		MyFill(table,"x2_cor","x2",600,-300,300,x2c);
		MyFill(table,"xavg_cor","xavg",600,-300,300,xac);
		MyFill(table,"xavg_theta_cor","xavg","theta",600,-300,300,xac,180,0,90,thetac*180.0/M_PI);
		MyFill(table,"tcheck1c","tcheck1",500,600,650,tcheck1c);
		MyFill(table,"tcheck2c","tcheck2",500,630,680,tcheck2c);
		MyFill(table,"tc1_dfrc","delayFR Time","tcheck1",1500,0,1500,dfr_tc,500,600,650, tcheck1c);
		MyFill(table,"tc1_dflc","delayFL Time","tcheck1",1500,0,1500,dfl_tc,500,600,650, tcheck1c);
		MyFill(table,"tc2_dbrc","delayBR Time","tcheck2",1500,0,1500,dbr_tc,500,630,680, tcheck2c);
		MyFill(table,"tc2_dblc","delayBL Time","tcheck2",1500,0,1500,dbl_tc,500,630,680, tcheck2c);
		MyFill(table,"tc1_x1c","x1","tcheck1",600,-300,300,x1c,500,600,650, tcheck1c);
		MyFill(table,"tc2_x2c","x2","tcheck2",600,-300,300,x2c,500,630,680,tcheck2c);
		MyFill(table,"dfr_tc_psd","delayFR Time","PSD",2000,-500,1500,dfr_tc,1000,0,1.0,dfr_psd);
		MyFill(table,"dfl_tc_psd","delayFL Time","PSD",2000,-500,1500,dfl_tc,1000,0,1.0,dfl_psd);
		MyFill(table,"dbr_tc_psd","delayBR Time","PSD",2000,-500,1500,dbr_tc,1000,0,1.0,dbr_psd);
		MyFill(table,"dbl_tc_psd","delayBL Time","PSD",2000,-500,1500,dbl_tc,1000,0,1.0,dbl_psd);
	}
	std::cout<<std::endl;
	std::cout<<"Finished. Writing to disk."<<std::endl;

	TGraph* graph1 = new TGraph(times.size(), &(times[0]), &(widths1[0]));
	graph1->SetName("anneal1_graph");
	graph1->SetTitle("Front Annealing Evo;time;#sigma_{tc1}");
	TGraph* graph2 = new TGraph(times.size(), &(times[0]), &(widths2[0]));
	graph2->SetName("anneal2_graph");
	graph2->SetTitle("Back Annealing Evo;time;#sigma_{tc2}");

	TGraph** graphs = new TGraph*[4];
	for(int i=0; i<4; i++) {
		auto& avgs = delayBinnedAvgPSD[i];
		std::string name = "delay"+std::to_string(i)+"_avgPSDgraph";
		graphs[i] = new TGraph(nbins, &(delaytimes[0]), &(avgs[0]));
		graphs[i]->SetName(name.c_str());
		graphs[i]->SetTitle(name.c_str());
		graphs[i]->SetMarkerColor(2);
	}

	input->Close();
	output->cd();
	table->Write();
	graph1->Write();
	graph2->Write();
	for(int i=0; i<4; i++) {
		graphs[i]->Write();
	}
	output->Close();
	delete event_add;
	delete graph1;
	delete graph2;
	for(int i=0; i<4; i++) {
		delete graphs[i];
	}
	delete[] graphs;
}

void WiggleKiller::ApplyLambdas(std::string& inname, std::string& outname) {
	if(!cuts.IsValid()) {
		std::cerr<<"Invalid cuts! Unable to run."<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inname.c_str(), "READ");
	TTree* intree = (TTree*) input->Get("SPSTree");

	ProcessedEvent* event_add = new ProcessedEvent();

	intree->SetBranchAddress("event", &event_add);

	TParameter<double>* at = (TParameter<double>*) input->Get("AT");
	TParameter<double>* zt = (TParameter<double>*) input->Get("ZT");
	TParameter<double>* ap = (TParameter<double>*) input->Get("AP");
	TParameter<double>* zp = (TParameter<double>*) input->Get("ZP");
	TParameter<double>* ae = (TParameter<double>*) input->Get("AE");
	TParameter<double>* ze = (TParameter<double>*) input->Get("ZE");
	TParameter<double>* b = (TParameter<double>*) input->Get("Bfield");
	TParameter<double>* bke = (TParameter<double>*) input->Get("BeamKE");
	TParameter<double>* theta = (TParameter<double>*) input->Get("Theta");
	weights = GetXavgWeights(zt->GetVal(), at->GetVal(), zp->GetVal(), ap->GetVal(), ze->GetVal(), ae->GetVal(), bke->GetVal(), theta->GetVal(), b->GetVal());

	TFile* output = TFile::Open(outname.c_str(), "RECREATE");

	THashTable* table = new THashTable();

	long blentries = intree->GetEntries();
	long flush = 0.1*blentries;
	long count=0, flush_count=0;

	double dfr_ti, dfl_ti, dbr_ti, dbl_ti;
	double dfr_psd, dfl_psd, dbr_psd, dbl_psd;
	double dfr_tc, dfl_tc, dbr_tc, dbl_tc;
	double x1c, x2c, xac, thetac, thetai;
	double tcheck1, tcheck2, tcheck1c, tcheck2c;
	int this_bin;

	for(int i=0; i<4; i++) {
		for(int j=0; j<nbins; j++) {
			delayBinnedAvgPSD[i][j] = 0.0;
			delayBinnedCounts[i][j] = 0.0;
			calibData[i][j].deltas.clear();
		}
	}

	std::vector<double> delaytimes;
	for(int i=0; i<nbins; i++) {
		delaytimes.push_back(i*BINWIDTH);
	}


	std::cout<<"Generating PSD correction data..."<<std::endl;

	for(long i=0; i<blentries; i++) {
		intree->GetEntry(i);
		count++;
		if(count == flush) {
			count =0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*10<<"%"<<std::flush;
		}

		if(event_add->xavg == -1e6) continue;

		if(!cuts.IsInside(event_add->cathode, event_add->scintLeft, event_add->x1, event_add->x2, event_add->anodeFrontTime, event_add->anodeBackTime)) continue;



		dfr_ti = event_add->delayFrontRightTime - event_add->anodeFrontTime;
		dfl_ti = event_add->delayFrontLeftTime - event_add->anodeFrontTime;
		dbr_ti = event_add->delayBackRightTime - event_add->anodeBackTime;
		dbl_ti = event_add->delayBackLeftTime - event_add->anodeBackTime;

		if(dfr_ti <0 || dfl_ti<0 || dbr_ti<0 || dbl_ti<0) {
			continue;
		}

		tcheck1 = (dfr_ti + dfl_ti)/2.0;
		tcheck2 = (dbr_ti + dbl_ti)/2.0;

		dfr_psd = (1.0 - event_add->delayFrontRightShort/event_add->delayFrontRightE);
		dfl_psd = (1.0 - event_add->delayFrontLeftShort/event_add->delayFrontLeftE);
		dbr_psd = (1.0 - event_add->delayBackRightShort/event_add->delayBackRightE);
		dbl_psd = (1.0 - event_add->delayBackLeftShort/event_add->delayBackLeftE);

		thetai = std::atan((event_add->x2 - event_add->x1)/(Wire_Dist()*10));

		MyFill(table,"x1psdl","x1","PSD Left",600,-300,300, event_add->x1,1000,0,1.0, dfl_psd);
		MyFill(table,"x1psdr","x1","PSD Right",600,-300,300, event_add->x1,1000,0,1.0, dfr_psd);
		MyFill(table,"x2psdl","x2","PSD Left",600,-300,300, event_add->x2,1000,0,1.0, dbl_psd);
		MyFill(table,"x2psdr","x2","PSD Right",600,-300,300, event_add->x2,1000,0,1.0, dbr_psd);
		MyFill(table,"x1","x1",600,-300,300,event_add->x1);
		MyFill(table,"x2","x2",600,-300,300,event_add->x2);
		MyFill(table,"xavg","xavg",600,-300,300,event_add->xavg);
		MyFill(table,"dfr_ti","delayFR Time",2000,-500,1500,dfr_ti);
		MyFill(table,"dfl_ti","delayFL Time",2000,-500,1500,dfl_ti);
		MyFill(table,"dbr_ti","delayBR Time",2000,-500,1500,dbr_ti);
		MyFill(table,"dbl_ti","delayBL Time",2000,-500,1500,dbl_ti);
		MyFill(table,"tcheck1","tcheck1",500,600,650,tcheck1);
		MyFill(table,"tcheck2","tcheck2",500,630,680,tcheck2);
		MyFill(table,"tc1_x1","x1","tcheck1",600,-300,300,event_add->x1,500,600,650, tcheck1);
		MyFill(table,"tc1_dfr","delayFR Time","tcheck1",1500,0,1500,dfr_ti,500,600,650, tcheck1);
		MyFill(table,"tc1_dfl","delayFL Time","tcheck1",1500,0,1500,dfl_ti,500,600,650, tcheck1);
		MyFill(table,"tc2_dbr","delayBR Time","tcheck2",1500,0,1500,dbr_ti,500,630,680, tcheck2);
		MyFill(table,"tc2_dbl","delayBL Time","tcheck2",1500,0,1500,dbl_ti,500,630,680, tcheck2);
		MyFill(table,"tc2_x2","x2","tcheck2",600,-300,300,event_add->x2,500,630,680, tcheck2);
		MyFill(table,"dfr_ti_psd","delayFR Time","PSD",2000,-500,1500,dfr_ti,1000,0,1.0,dfr_psd);
		MyFill(table,"dfl_ti_psd","delayFL Time","PSD",2000,-500,1500,dfl_ti,1000,0,1.0,dfl_psd);
		MyFill(table,"dbr_ti_psd","delayBR Time","PSD",2000,-500,1500,dbr_ti,1000,0,1.0,dbr_psd);
		MyFill(table,"dbl_ti_psd","delayBL Time","PSD",2000,-500,1500,dbl_ti,1000,0,1.0,dbl_psd);
		MyFill(table,"xavg_theta","xavg","theta",600,-300,300,event_add->xavg,180,0,90, thetai*180.0/M_PI);

		this_bin = floor(dfr_ti/BINWIDTH);
		delayBinnedAvgPSD[FrontRight][this_bin] += dfr_psd;
		delayBinnedCounts[FrontRight][this_bin]++;


		this_bin = floor(dfl_ti/BINWIDTH);
		delayBinnedAvgPSD[FrontLeft][this_bin] += dfl_psd;
		delayBinnedCounts[FrontLeft][this_bin]++;
		
		this_bin = floor(dbr_ti/BINWIDTH);
		delayBinnedAvgPSD[BackRight][this_bin] += dbr_psd;
		delayBinnedCounts[BackRight][this_bin]++;

		this_bin = floor(dbl_ti/BINWIDTH);
		delayBinnedAvgPSD[BackLeft][this_bin] += dbl_psd;
		delayBinnedCounts[BackLeft][this_bin]++;
	}
	std::cout<<std::endl;

	for(int i=0; i<4; i++) {
		for(int j=0; j<nbins; j++) {
			if(delayBinnedCounts[i][j] == 0) continue;
			delayBinnedAvgPSD[i][j] /= delayBinnedCounts[i][j];
		}
	}

	PSDsplines[FrontLeft].ReadData(delaytimes, delayBinnedAvgPSD[FrontLeft]);
	PSDsplines[FrontRight].ReadData(delaytimes, delayBinnedAvgPSD[FrontRight]);
	PSDsplines[BackLeft].ReadData(delaytimes, delayBinnedAvgPSD[BackLeft]);
	PSDsplines[BackRight].ReadData(delaytimes, delayBinnedAvgPSD[BackRight]);
	
	std::cout<<"Finished. Applying PSD correction..."<<std::endl;
	std::cout<<"Generating final histograms..."<<std::endl;
	
	flush_count=0;
	count=0;
	for(long i=0; i<blentries; i++) {
		intree->GetEntry(i);
		count++;
		if(count == flush) {
			count =0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*10<<"%"<<std::flush;
		}

		if(event_add->xavg == -1e6) continue;
		if(!cuts.IsInside(event_add->cathode, event_add->scintLeft, event_add->x1, event_add->x2, event_add->anodeFrontTime, event_add->anodeBackTime)) continue;

		
		dfr_ti = event_add->delayFrontRightTime - event_add->anodeFrontTime;
		dfl_ti = event_add->delayFrontLeftTime - event_add->anodeFrontTime;
		dbr_ti = event_add->delayBackRightTime - event_add->anodeBackTime;
		dbl_ti = event_add->delayBackLeftTime - event_add->anodeBackTime;

		dfr_psd = (1.0 - event_add->delayFrontRightShort/event_add->delayFrontRightE);
		dfl_psd = (1.0 - event_add->delayFrontLeftShort/event_add->delayFrontLeftE);
		dbr_psd = (1.0 - event_add->delayBackRightShort/event_add->delayBackRightE);
		dbl_psd = (1.0 - event_add->delayBackLeftShort/event_add->delayBackLeftE);

		dfr_tc = ApplyPSDCorrection(dfr_ti, dfr_psd, FrontRight, lambdaFront);
		dfl_tc = ApplyPSDCorrection(dfl_ti, dfl_psd, FrontLeft, lambdaFront);
		dbr_tc = ApplyPSDCorrection(dbr_ti, dbr_psd, BackRight, lambdaBack);
		dbl_tc = ApplyPSDCorrection(dbl_ti, dbl_psd, BackLeft, lambdaBack);

		this_bin = floor(dfr_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[FrontRight][this_bin].deltas.push_back(dfr_tc - dfr_ti);

		this_bin = floor(dfl_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[FrontLeft][this_bin].deltas.push_back(dfl_tc - dfl_ti);

		this_bin = floor(dbl_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[BackLeft][this_bin].deltas.push_back(dbl_tc - dbl_ti);

		this_bin = floor(dbr_ti/BINWIDTH);
		if(this_bin >= 0 && this_bin < nbins) calibData[BackRight][this_bin].deltas.push_back(dbr_tc - dbr_ti);

		x1c = (dfl_tc - dfr_tc)*0.5*0.505;
		x2c = (dbl_tc - dbr_tc)*0.5*0.476;
		xac = x1c*weights.first + x2c*weights.second;
		thetac = std::atan((x2c - x1c)/(Wire_Dist()*10.0));
		
		tcheck1c = (dfr_tc + dfl_tc)/2.0;
		tcheck2c = (dbr_tc + dbl_tc)/2.0;


		MyFill(table,"x1psdl_cor","x1","PSD Left",600,-300,300,x1c,1000,0.0,1.0, dfl_psd);
		MyFill(table,"x1psdr_cor","x1","PSD Right",600,-300,300,x1c,1000,0.0,1.0, dfr_psd);
		MyFill(table,"x2psdl_cor","x2","PSD Left",600,-300,300,x2c,1000,0.0,1.0, dbl_psd);
		MyFill(table,"x2psdr_cor","x2","PSD Right",600,-300,300,x1c,1000,0.0,1.0, dbr_psd);
		MyFill(table,"x1_cor","x1",600,-300,300,x1c);
		MyFill(table,"x2_cor","x2",600,-300,300,x2c);
		MyFill(table,"xavg_cor","xavg",600,-300,300,xac);
		MyFill(table,"xavg_theta_cor","xavg","theta",600,-300,300,xac,180,0,90,thetac*180.0/M_PI);
		MyFill(table,"tcheck1c","tcheck1",500,600,650,tcheck1c);
		MyFill(table,"tcheck2c","tcheck2",500,630,680,tcheck2c);
		MyFill(table,"tc1_dfrc","delayFR Time","tcheck1",1500,0,1500,dfr_tc,500,600,650, tcheck1c);
		MyFill(table,"tc1_dflc","delayFL Time","tcheck1",1500,0,1500,dfl_tc,500,600,650, tcheck1c);
		MyFill(table,"tc2_dbrc","delayBR Time","tcheck2",1500,0,1500,dbr_tc,500,630,680, tcheck2c);
		MyFill(table,"tc2_dblc","delayBL Time","tcheck2",1500,0,1500,dbl_tc,500,630,680, tcheck2c);
		MyFill(table,"tc1_x1c","x1","tcheck1",600,-300,300,x1c,500,600,650, tcheck1c);
		MyFill(table,"tc2_x2c","x2","tcheck2",600,-300,300,x2c,500,630,680,tcheck2c);
		MyFill(table,"dfr_tc_psd","delayFR Time","PSD",2000,-500,1500,dfr_tc,1000,0,1.0,dfr_psd);
		MyFill(table,"dfl_tc_psd","delayFL Time","PSD",2000,-500,1500,dfl_tc,1000,0,1.0,dfl_psd);
		MyFill(table,"dbr_tc_psd","delayBR Time","PSD",2000,-500,1500,dbr_tc,1000,0,1.0,dbr_psd);
		MyFill(table,"dbl_tc_psd","delayBL Time","PSD",2000,-500,1500,dbl_tc,1000,0,1.0,dbl_psd);
	}
	std::cout<<std::endl;
	std::cout<<"Finished. Writing to disk."<<std::endl;

	TGraph** graphs = new TGraph*[4];
	for(int i=0; i<4; i++) {
		auto& avgs = delayBinnedAvgPSD[i];
		std::string name = "delay"+std::to_string(i)+"_avgPSDgraph";
		graphs[i] = new TGraph(nbins, &(delaytimes[0]), &(avgs[0]));
		graphs[i]->SetName(name.c_str());
		graphs[i]->SetTitle(name.c_str());
		graphs[i]->SetMarkerColor(2);
	}
	TF1* flPSD = new TF1("flPSD",&(PSDsplines[FrontLeft]), &CubicSpline::EvaluateROOT, 200, 1100,0, "CubicSpline", "EvaluateROOT");
	TF1* frPSD = new TF1("frPSD",&(PSDsplines[FrontRight]), &CubicSpline::EvaluateROOT, 200, 1100,0, "CubicSpline", "EvaluateROOT");
	TF1* blPSD = new TF1("blPSD",&(PSDsplines[BackLeft]), &CubicSpline::EvaluateROOT, 200, 1100,0, "CubicSpline", "EvaluateROOT");
	TF1* brPSD = new TF1("brPSD",&(PSDsplines[BackRight]), &CubicSpline::EvaluateROOT, 200, 1100,0, "CubicSpline", "EvaluateROOT");

	TCanvas* c1 = new TCanvas();
	flPSD->Draw();

	TCanvas* c2 = new TCanvas();
	frPSD->Draw();

	TCanvas* c3 = new TCanvas();
	blPSD->Draw();

	TCanvas* c4 = new TCanvas();
	brPSD->Draw();

	input->Close();
	output->cd();
	table->Write();
	c1->Write();
	c2->Write();
	c3->Write();
	c4->Write();
	for(int i=0; i<4; i++) {
		graphs[i]->Write();
	}
	output->Close();
	delete event_add;
	for(int i=0; i<4; i++) {
		delete graphs[i];
	}
	delete[] graphs;
}

void WiggleKiller::GenerateCalibrationFiles(std::string& dflname, std::string& dfrname, std::string& dblname, std::string& dbrname) {

	std::ofstream dflfile(dflname);
	std::ofstream dfrfile(dfrname);
	std::ofstream dblfile(dblname);
	std::ofstream dbrfile(dbrname);

	if(!(dflfile.is_open() && dfrfile.is_open() && dblfile.is_open() && dbrfile.is_open())) {
		std::cerr<<"Error at WiggleKiller::GenerateCalibrationFiles! Unable to open one of the files named: "<<dflname<<" "<<dfrname<<" "<<dblname<<" "<<dbrname<<std::endl;
		return;
	}

	int flcounter=0, frcounter=0, blcounter=0, brcounter=0; 
	for(int i=0; i<nbins; i++) {
		if(calibData[FrontRight][i].GetDelta() != 0.0) {
			dfrfile<<calibData[FrontRight][i].time_i<<" "<<calibData[FrontRight][i].GetDelta()<<std::endl;
			frcounter++;
		} 
		if(calibData[FrontLeft][i].GetDelta() != 0.0) {
			dflfile<<calibData[FrontLeft][i].time_i<<" "<<calibData[FrontLeft][i].GetDelta()<<std::endl;
			flcounter++;
		}
		if(calibData[BackRight][i].GetDelta() != 0.0) {
			dbrfile<<calibData[BackRight][i].time_i<<" "<<calibData[BackRight][i].GetDelta()<<std::endl;
			brcounter++;
		}
		if(calibData[BackLeft][i].GetDelta() != 0.0) {
			dblfile<<calibData[BackLeft][i].time_i<<" "<<calibData[BackLeft][i].GetDelta()<<std::endl;
			blcounter++;
		}
	}

	dfrfile.close();
	dflfile.close();
	dbrfile.close();
	dblfile.close();

	if(frcounter == 0 || flcounter == 0 || brcounter == 0 || blcounter == 0) {
		std::cerr<<"Unexpected behavior at WiggleKiller::GenerateCalibrationFiles! No data written to calibration files... Make sure either OptimizeParameters or ApplyLambdas is run first!"<<std::endl;
	}
}