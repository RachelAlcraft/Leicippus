#pragma once
/************************************************************************
* RSA 28/11/21
* Static functions manage conversions between various electron density and structure factors
in theoretical or experimental state
************************************************************************/

#include <string>
#include <vector>
#include "StructureFactors.h"
#include "ElectronDensity.h"
#include "Protein.h"

using namespace std;

class ED_SF
{
public:
    ElectronDensity* EdTheo;
    StructureFactorsTheoretical* SfTheo;
    double Resolution;
private:
    VectorThree _startHKL;
    VectorThree _endHKL;
    
public:
    ED_SF();
    void proteinToTheoreticalED(string name, Protein* prot);
    void proteinToTheoreticalSF(string name, Protein* prot);
    
    
    
};