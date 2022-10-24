#ifndef CUTHANDLER_H
#define CUTHANDLER_H

#include <string>
#include <TFile.h>
#include <TCutG.h>

class CutHandler
{
public:
    CutHandler();
    CutHandler(const std::string& edename, const std::string& xxname);
    ~CutHandler();
    void SetCuts(const std::string& edename, const std::string& xxname);
    bool IsValid() { return validFlag; };
    bool IsInside(double cath, double scint, double x1, double x2, double anodeFTime, double anodeBTime);

private:
    TFile *edefile, *xxfile;
    TCutG *edecut, *xxcut;

    bool validFlag;
};

#endif