#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "WiggleKiller.h"
#include "WiggleCorrector.h"


int main(int argc, char** argv) {
	if(argc < 2) {
		std::cerr<<"Incorrect number of commandline arguments! Expects at least one"<<std::endl;
		return 1;
	}
	std::string option;
	std::string configname;
	if(argc == 3) {
		option = argv[1];
		configname = argv[2];
	} else {
		configname = argv[1];
		option = "";
	}

	std::ifstream config(configname);
	if(!config.is_open()) {
		std::cerr<<"Unable to open configuration! Exiting without doing anything."<<std::endl;
		return 1;
	}

	/*
		Valid program options:

		--kill: takes lambda from input file, applies them to data set and saves some histograms
		--optimize: takes lambda from input file, applies the optimization routine to find best lambdas and saves some histograms
		--clean: takes calibration files and applies the calibration to data. Saves histograms and a tree.
		--clean-no-tree: takes calibration files and applies the calibration to data. Saves histograms.
		--spline-test: simply a test to see if the splines actually work.
	*/

	std::string junk, input_optimize, output_optimize, edename, xxname, input_bulk, output_bulk;
	std::string t1lfile, t1rfile, t2lfile, t2rfile;
	double lambdaFront, lambdaBack;

	config>>junk>>input_optimize;
	config>>junk>>output_optimize;
	config>>junk>>edename;
	config>>junk>>xxname;
	config>>junk>>t1lfile;
	config>>junk>>t1rfile;
	config>>junk>>t2lfile;
	config>>junk>>t2rfile;
	config>>junk>>lambdaFront;
	config>>junk>>lambdaBack;
	config>>junk>>input_bulk;
	config>>junk>>output_bulk;


	config.close();

	if(option == "--kill") {
		std::cout<<std::setprecision(4);
		std::cout<<"Running wiggle elimination routine..."<<std::endl;
		std::cout<<"Input file: "<<input_optimize<<" Output file: "<<output_optimize<<std::endl;
		std::cout<<"LambdaFront: "<<lambdaFront<<" LambdaBack: "<<lambdaBack<<std::endl;
		WiggleKiller steveBrown(edename, xxname, lambdaFront, lambdaBack);
		steveBrown.ApplyLambdas(input_optimize, output_optimize);
		steveBrown.GenerateCalibrationFiles(t1lfile, t1rfile, t2lfile, t2rfile);
	} else if(option == "--optimize") {
		std::cout<<std::setprecision(4);
		std::cout<<"Running psd parameter optimization routine..."<<std::endl;
		std::cout<<"Input file: "<<input_optimize<<" Output file: "<<output_optimize<<std::endl;
		std::cout<<"LambdaFront: "<<lambdaFront<<" LambdaBack: "<<lambdaBack<<std::endl;
		WiggleKiller steveBrown(edename, xxname, lambdaFront, lambdaBack);
		steveBrown.OptimizeParameters(input_optimize, output_optimize);
		steveBrown.GenerateCalibrationFiles(t1lfile, t1rfile, t2lfile, t2rfile);
		std::cout<<"Finished."<<std::endl;
	} else if(option == "--clean") {
		std::cout<<"Cleaning up data..."<<std::endl;
		std::cout<<"Input file: "<<input_bulk<<" Output file: "<<output_bulk<<std::endl;
		WiggleCorrector fixer(t1lfile, t1rfile, t2lfile, t2rfile);
		fixer.SetCuts(edename, xxname);
		fixer.ApplyCorrection(input_bulk, output_bulk);
		std::cout<<"Finished."<<std::endl;
	} else if (option == "--clean-no-tree") {
		std::cout<<"Cleaning up data and only writing histograms..."<<std::endl;
		std::cout<<"Input file: "<<input_bulk<<" Output file: "<<output_bulk<<std::endl;
		WiggleCorrector fixer(t1lfile, t1rfile, t2lfile, t2rfile);
		fixer.SetCuts(edename, xxname);
		fixer.ApplyCorrection_NoTree(input_bulk, output_bulk);
		std::cout<<"Finished."<<std::endl;
	} else if(option == "--spline-test") {
		CubicSpline test(t1lfile);
		double value = test.Evaluate(500.0);
		std::cout<<"Testing cubic spline class... For initial time 500.0 calculate a delta of "<<value<<std::endl;
	}
}