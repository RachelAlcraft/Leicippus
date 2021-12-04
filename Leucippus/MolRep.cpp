/*
* RSA 14/11/2021
* The very most basic possioble attempt to do molecular replacement
* All I want to do is create phases as if I knew the structure, and I will start with the sf and the pdb itself
*/

#include "MolRep.h"

#include <iostream>

MolRep::MolRep(string coordsFile, string structureFactorsFile)
{
	/*_cifStructureFactors = new CifFile(structureFactorsFile);
	_cifAtomCoordinates = new CifFile(coordsFile);
	_sfExp = new StructureFactorsExperimental(_cifStructureFactors);
	_structure = new Protein(_cifAtomCoordinates->LoopElements["_atom_site"]);
	string id = _cifAtomCoordinates->NonLoopElements["_entry"]["id"];
	EdMolRep = new ElectronDensity(id);
	int a = atol(_cifStructureFactors->NonLoopElements["_cell"]["length_a"].c_str());
	int b = atol(_cifStructureFactors->NonLoopElements["_cell"]["length_b"].c_str());
	int c = atol(_cifStructureFactors->NonLoopElements["_cell"]["length_c"].c_str());
	double alpha = atof(_cifStructureFactors->NonLoopElements["_cell"]["angle_alpha"].c_str());
	double beta = atof(_cifStructureFactors->NonLoopElements["_cell"]["angle_beta"].c_str());
	double gamma = atof(_cifStructureFactors->NonLoopElements["_cell"]["angle_gamma"].c_str());
	_lattice = RealLattice(a, b, c, alpha, beta, gamma);		
	cout << "Calculating molecular replacement electron density\n";
	for (unsigned int x = 0; x < _sfExp->SFs.size(); ++x)
	{
		StructureFactor* sf = _sfExp->SFs[x].Second;
		cout << "..." << sf->HKL.A << "," << sf->HKL.B << "," << sf->HKL.C << "\n";
		VectorThree v3;
		VectorThree mnp = v3.getCrystalIndices(sf->HKL.A, sf->HKL.B, sf->HKL.C);
		VectorThree xyz = _lattice.OrthoToCarte(mnp);
		double ed = 0;
		for (unsigned int y = 0; y < _structure->Atoms.size(); ++y)
		{
			Atom* atm = _structure->Atoms[y];
			ed += atm->electronDensityContributionTheoretical(xyz);
		}
		EdMolRep->addVoxel(mnp.A, mnp.B, mnp.C, ed,x);
	}
	_sfTheo = new StructureFactorsTheoretical(id,_sfExp->hStart, _sfExp->hEnd, _sfExp->kStart, _sfExp->kEnd, _sfExp->lStart, _sfExp->lEnd, _structure->Atoms);
	*/

}

