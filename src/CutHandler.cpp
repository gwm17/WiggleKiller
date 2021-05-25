#include "CutHandler.h"

CutHandler::CutHandler() :
	edefile(nullptr), xxfile(nullptr), edecut(nullptr), xxcut(nullptr), validFlag(false)
{
}

CutHandler::CutHandler(std::string& edename, std::string& xxname) :
	edefile(nullptr), xxfile(nullptr), edecut(nullptr), xxcut(nullptr), validFlag(false)
{
	SetCuts(edename, xxname);
}

CutHandler::~CutHandler() {
	if(edefile != nullptr && edefile->IsOpen()) edefile->Close();
	if(xxfile != nullptr && xxfile->IsOpen()) xxfile->Close();
}

void CutHandler::SetCuts(std::string& edename, std::string& xxname) {
	edefile = TFile::Open(edename.c_str(), "READ");
	if(edefile->IsOpen()) {
		edecut = (TCutG*) edefile->Get("CUTG");
		edecut->SetName("EdECut");
	}

	xxfile = TFile::Open(xxname.c_str(), "READ");
	if(xxfile->IsOpen()) {
		xxcut = (TCutG*) xxfile->Get("CUTG");
		xxcut->SetName("X1X2Cut");
	}

	if(edecut != nullptr && xxcut != nullptr) {
		validFlag = true;
	} else {
		validFlag = false;
	}
}

bool CutHandler::IsInside(double cath, double scint, double x1, double x2, double anodeFTime, double anodeBTime) {
	if(edecut->IsInside(scint, cath) && xxcut->IsInside(x1, x2) && anodeFTime != -1 && anodeBTime != -1) {
		return true;
	} else {
		return false;
	}
}