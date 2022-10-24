/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

*/
#include "MassLookup.h"
#include <fstream>
#include <iostream>


//Static instance lives for the lifetime of the program
MassLookup* MassLookup::s_instance = new MassLookup();

/*
  Read in AMDC mass file, preformated to remove excess info. Here assumes that by default
  the file is in a local directory etc/
*/
MassLookup::MassLookup()
{
    std::ifstream massfile("./etc/mass.txt");
    if (massfile.is_open())
    {
        std::string junk, element;
        double atomicMassBig, atomicMassSmall;
        NucData isotope;
        getline(massfile, junk);
        getline(massfile, junk);
        while (massfile >> junk)
        {
            massfile >> isotope.Z >> isotope.A >> element >> atomicMassBig >> atomicMassSmall;
            isotope.mass = (atomicMassBig + atomicMassSmall * 1e-6 - isotope.Z * s_eMass) * s_u2MeV;
            isotope.symbol = std::to_string(isotope.A) + element;
            m_dataMap[SzudzikPairing(isotope.Z, isotope.A)] = isotope;
        }
    }
    else
    {
        std::cerr << "Unable to open mass.txt at MassLookup! Prepare for errors." << std::endl;
    }
}

MassLookup::~MassLookup() {}

// Returns nuclear mass in MeV
double MassLookup::FindMass(uint64_t Z, uint64_t A) const
{
    auto data = m_dataMap.find(SzudzikPairing(Z, A));
    if (data == m_dataMap.end())
    {
        std::cerr << "Invaild nucleus at MassLookup! Returning mass of 0" << std::endl;
        return 0;
    }
    return data->second.mass;
}

// returns element symbol
std::string MassLookup::FindSymbol(uint64_t Z, uint64_t A) const
{
    auto data = m_dataMap.find(SzudzikPairing(Z, A));
    if (data == m_dataMap.end())
    {
        std::cerr << "Invaild nucleus at MassLookup! Returning empty symbol" << std::endl;
        return "";
    }
    return data->second.symbol;
}
