
#pragma once
/*
* RSA 14/11/2021
* The very most basic possioble attempt to do molecular replacement
* All I want to do is create phases as if I knew the structure, and I will start with the sf and the pdb itself
*/

#include <string>
#include <vector>
#include "CifFile.h"
#include "StructureFactors.h"
#include "Protein.h"
#include "ElectronDensity.h"

using namespace std;

class MolRep
{
private:
	CifFile* _cifStructureFactors;
	CifFile* _cifAtomCoordinates;
	Protein* _structure;
	StructureFactorsExperimental* _sfExp;
	StructureFactorsTheoretical* _sfTheo;
	RealLattice _lattice;
public:
	ElectronDensity* EdMolRep;

public:
	MolRep(string coordsFile, string structureFactorsFile);

};