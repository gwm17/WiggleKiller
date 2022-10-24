/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

*/
#ifndef MASS_LOOKUP_H
#define MASS_LOOKUP_H

#include <string>
#include <unordered_map>

static constexpr uint64_t SzudzikPairing(uint64_t i, uint64_t j)
{
    return i >= j ? i*i + i + j : j*j +i;
}

class MassLookup
{
public:

    struct NucData
    {
        std::string symbol = "";
        double mass = 0.0;
        uint64_t Z = 0;
        uint64_t A = 0;
    };

    static const MassLookup& GetInstance() { return *s_instance; }

    ~MassLookup();
    double FindMass(uint64_t Z, uint64_t A) const;
    std::string FindSymbol(uint64_t Z, uint64_t A) const;

private:
    MassLookup();

    static MassLookup* s_instance;

    std::unordered_map<uint64_t, NucData> m_dataMap;
    //constants
    static constexpr double s_u2MeV = 931.4940954;
    static constexpr double s_eMass = 0.000548579909;
    
};



#endif
