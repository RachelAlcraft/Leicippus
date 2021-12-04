
#include "StructureFactors.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


void StructureFactors::addStructureFactor(StructureFactor* sf)
{
	//we need to accumulate the structure factor onto the key
	stringstream id;
	id << sf->HKL.A << "_" << sf->HKL.B << "_" << sf->HKL.C;
	map<string, StructureFactor*>::iterator iter = SFs.find(id.str());
	if (iter == SFs.end())
	{
		SFs.insert(pair<string, StructureFactor*>(id.str(), sf));
	}
	else
	{
		iter->second->ComplexSF += sf->ComplexSF;
		iter->second->Intensity += sf->Intensity;
	}
	
}
void StructureFactors::print(string filename)
{
	cout << "\nPrinting to cif...\n";
	ofstream sffile;
	sffile.open(filename.c_str());
	sffile << "data_LeucipPlus_SF_" << Name << "\n";
	sffile << "#\n";
	sffile << "loop_\n";
	//edfile << "_ed.id\n";
	sffile << "_sf.index_h\n";
	sffile << "_sf.index_k\n";
	sffile << "_sf.index_l\n";
	sffile << "_sf.intensity\n";
	int i = 0;
	for (map<string, StructureFactor*>::iterator iter = SFs.begin(); iter != SFs.end(); ++iter)
	{
		if (i % 25000 == 0)		
			cout << "...Printing " << i << "/" << SFs.size() << "\n";
		++i;
				
		stringstream out;
		out << setprecision(5);
		//out << iter->first << " ";
		
		if (abs(iter->second->Intensity) > 0.00001)
			out << iter->second->HKL.A << " " << iter->second->HKL.B << " " << iter->second->HKL.C << " " << iter->second->Intensity << "\n";
		sffile << out.str();
	}
	sffile.close();
}

/*******************************************
* EXPERIMENTAL
*******************************************/
StructureFactorsExperimental::StructureFactorsExperimental(CifFile* cif)
{
	_cif = cif;
	Name = _cif->NonLoopElements["_symmetry"]["entry_id"];
	hStart = 0;
	hEnd = 0;
	kStart = 0;
	kEnd = 0;
	lStart = 0;
	lEnd = 0;
	_waveLength = 0;
	size_t num = cif->LoopElements["_refln"]["index_h"].size();
	for (int i = 0; i < int(num); ++i)
	{
		int h = atol(cif->LoopElements["_refln"]["index_h"][i].c_str());
		int k = atol(cif->LoopElements["_refln"]["index_k"][i].c_str());
		int l = atol(cif->LoopElements["_refln"]["index_l"][i].c_str());
		string intsy = cif->LoopElements["_refln"]["F_meas_au"][i];

		StructureFactor* sf = new StructureFactor(VectorThree(h, k, l), atof(intsy.c_str()));
		addStructureFactor(sf);

		if (h < hStart)
			hStart = h;
		if (h > hEnd)
			hEnd = h;
		if (k < kStart)
			kStart = k;
		if (k > kEnd)
			kEnd = k;
		if (l < lStart)
			lStart = l;
		if (l > lEnd)
			lEnd = l;
	}

}
void StructureFactorsExperimental::printRealConversion(string filename)
{
	/*
	map<string, double> real_space_sfs;
	VectorThree v3;
	
	for (unsigned int i = 0; i < SFs.size(); ++i)
	{
		StructureFactor* sf = SFs[i];		
		VectorThree mnp = v3.getCrystalIndices(sf->HKL.A, sf->HKL.B, sf->HKL.C);		
		stringstream id;
		id << mnp.A << " " << mnp.B << " " << mnp.C;
		map<string, double>::iterator iter = real_space_sfs.find(id.str());
		if (iter == real_space_sfs.end())
		{
			real_space_sfs.insert(pair<string, double>(id.str(), sf->Intensity));
		}
		else
		{
			iter->second += sf->Intensity;
		}
	}

	ofstream edfile;
	edfile.open(filename.c_str());
	edfile << "data_LeucipPy_SF_RealSpace" << Name << "\n";
	edfile << "#\n";
	edfile << "loop_\n";
	edfile << "_sfreal.index_m\n";
	edfile << "_sfreal.index_n\n";
	edfile << "_sfreal.index_p\n";
	edfile << "_sfreal.realspace\n";
	for (map<string, double>::iterator iter = real_space_sfs.begin(); iter != real_space_sfs.end(); ++iter)
	{
		stringstream out;
		out << setprecision(5);
		out << iter->first << " " << iter->second << "\n";
		edfile << out.str();
	}
	edfile.close();
	*/

}
RealLattice StructureFactorsExperimental::getLattice()
{
	int a = atol(_cif->NonLoopElements["_cell"]["length_a"].c_str());
	int b = atol(_cif->NonLoopElements["_cell"]["length_b"].c_str());
	int c = atol(_cif->NonLoopElements["_cell"]["length_c"].c_str());
	double alpha = atof(_cif->NonLoopElements["_cell"]["angle_alpha"].c_str());
	double beta = atof(_cif->NonLoopElements["_cell"]["angle_beta"].c_str());
	double gamma = atof(_cif->NonLoopElements["_cell"]["angle_gamma"].c_str());
	return RealLattice(a, b, c, alpha, beta, gamma);
}
/*******************************************
* THEORETICAL
*******************************************/
StructureFactorsTheoretical::StructureFactorsTheoretical(string name, int h_start, int h_end, int k_start, int k_end, int l_start, int l_end, vector<Atom*> atoms)
{//the model for structure factors
	cout << "Calculating Theoretical Structure factors for\n";
	cout << "h :" << h_start << " to " << h_end << "\n";
	cout << "k :" << k_start << " to " << k_end << "\n";
	cout << "l :" << l_start << " to " << l_end << "\n";
	for (int h = h_start; h <= h_end; ++h)
	{
		for (int k = k_start; k <= k_end; ++k)
		{
			for (int l = l_start; l <= l_end; ++l)
			{
				cout << "..." << h << "," << k << "," << l << "\n";
				VectorThree hkl = VectorThree(h, k, l);
				complex<double> im_sf = 0;
				for (unsigned int a = 0; a < atoms.size(); ++a)
				{
					
					complex<double> im_sfa = atoms[a]->structureFactorContributionTheoretical(hkl);
					im_sf = im_sf + im_sfa;
				}
				StructureFactor *sf = new StructureFactor(hkl, im_sf);
				addStructureFactor(sf);
				
			}
		}
	}
}
StructureFactorsTheoretical::StructureFactorsTheoretical(string name)
{
	Name = name;
}

/************************************************************************************************************************************
* Class for a single structure factor
***************************************************************************************************************************************/
StructureFactor::StructureFactor(VectorThree hkl, double intensity)
{
	HKL = hkl;	
	Intensity = intensity;
	Experimental = true;
}
StructureFactor::StructureFactor(VectorThree hkl, complex<double> sf)
{
	HKL = hkl;
	ComplexSF = sf;
	Intensity = sf.real();
	Experimental = false;
}

StructureFactors::StructureFactors()
{
	hStart = 0;
	hEnd = 0;
	kStart = 0;
	kEnd = 0;
	lStart = 0;
	lEnd = 0;
}
