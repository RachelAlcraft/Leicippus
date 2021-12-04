#include<sstream>
#include <iostream>
#include<iomanip>

#include "ED_SF.h"


using namespace std;

//https://www-structmed.cimr.cam.ac.uk/Course/Adv_diff2/Diffraction2.html#reciprocal_space

ED_SF::ED_SF()
{
	Resolution = 0.5;
	_startHKL = VectorThree(-20, -20, -20);
	_endHKL = VectorThree(20, 20, 20);
}

void ED_SF::proteinToTheoreticalED(string name, Protein* prot)
{
	cout << "Calculating theoretical electron density\n";
	// We can choose any kind of lattice and resolution we want
	// Make it orthogonal and A=B=C	
	
	//I need the boundaries of the molecule.
	int xmin = int(prot->MinX - 1);
	int ymin = int(prot->MinY - 1);
	int zmin = int(prot->MinZ - 1);
	int xmax = int(prot->MaxX + 1);
	int ymax = int(prot->MaxY + 1);
	int zmax = int(prot->MaxZ + 1);
	double lengthX = float(xmax) - float(xmin);
	double lengthY = float(ymax) - float(ymin);
	double lengthZ = float(zmax) - float(zmin);

	int numC = (lengthX / Resolution)+1;
	int numR = (lengthY / Resolution)+1;
	int numS = (lengthZ / Resolution)+1;

	EdTheo = new ElectronDensity(name);
	EdTheo->NumsCRS = VectorThree(numC, numR, numS);
	EdTheo->NumsXYZ = VectorThree(numC, numR, numS);
	EdTheo->StartsCRS = VectorThree(xmin, ymin, zmin);
	EdTheo->LengthABC = VectorThree(lengthX,lengthY,lengthZ);
	EdTheo->calculateOrthoMats();
	
	int count = 0;
	for (double x = xmin; x <= xmax; x += Resolution)
	{
		for (double y = ymin; y <= ymax; y += Resolution)
		{
			for (double z = zmin; z <= zmax; z += Resolution)
			{
				double edall = 0;
				complex<double> sfall;
				for (unsigned int a = 0; a < prot->Atoms.size(); ++a)
				{
					Atom* atm = prot->Atoms[a];					
					edall += atm->electronDensityContributionTheoretical(VectorThree(x, y, z));										
				}
				VectorThree crs = EdTheo->getCRSFromXYZ(VectorThree(x, y, z));
				
				
				EdTheo->addVoxel(int(crs.A), int(crs.B), int(crs.C), edall, count);
																 
				count++;
				if (count % 25000)
					cout << "... theoretical (" << x << "," << y << "," << z << ") out of (" << xmax << "," << ymax << "," << zmax << ")\n";
			}
		}
	}	
}

void ED_SF::proteinToTheoreticalSF(string name, Protein* prot)
{
	cout << "Calculating theoretical structure factors\n";
	// We can choose any kind of lattice and resolution we want
		
	SfTheo = new StructureFactorsTheoretical(name);

	int count = 0;
	for (double hh = _startHKL.A; hh <= _endHKL.A; ++hh)
	{
		for (double kk = _startHKL.B; kk <= _endHKL.B; ++kk)
		{
			for (double ll = _startHKL.C; ll <= _endHKL.C; ++ll)
			{
				
				complex<double> sfall;
				VectorThree hkl(hh, kk, ll);
				VectorThree hkl_recip = EdTheo->LattFractional.FracToRecip(VectorThree(hh, kk, ll));
				for (unsigned int v = 0; v < EdTheo->Voxels.size(); ++v)
				{
					Voxel* vox = EdTheo->Voxels[v];
					sfall += EdTheo->structureFactorContributionTheoretical(hkl_recip, vox->CRS, vox->v);
					//sfall += EdTheo->structureFactorContributionTheoretical(hkl, vox->CRS, vox->v);
				}
				StructureFactor* sf = new StructureFactor(hkl, sfall);
				
				SfTheo->addStructureFactor(sf);

				count++;
				if (count % 250)
					cout << "... theoretical (" << hh << "," << kk << "," << ll << "\n";
			}
		}
	}
}



