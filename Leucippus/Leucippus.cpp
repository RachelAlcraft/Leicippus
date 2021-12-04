// Leucippus.cpp : Defines the entry point for the application.
//

#include "Leucippus.h"
#include "Lattice.h"
#include "CifFile.h"
#include "StructureFactors.h"
#include "MolRep.h"
#include "Ccp4File.h"
#include "ED_SF.h"

using namespace std;

int main()
{
	//old();	
	//convertCcp4ToCif();
	calculateTheoreticalElectronDensity();
	return 0;
}

void calculateSlicesList()
{
	string outdir = 
	vector<string> slices;
	//pdb_code, chain1, aa1, rid1, atom1, chain2, aa2, rid2, atom2, givendistance, comment				
	slices.push_back("6zwh, 6zwh, -27.327, 17.986, -15.654, -28.508, 19.956, -14.196, -29.514, 20.298, -13.152");
	slices.push_back("6zwj, 6zwj, -27.321, 17.941, -15.658, -28.456, 19.927, -14.195, -29.465, 20.298, -13.168");
	slices.push_back("6zx4, 6zx4, -27.333, 17.934, -15.65, -28.501, 19.873, -14.129, -29.504, 20.269, -13.113");
	slices.push_back("7b0l, 7b0l, -27.292, 17.909, -15.493, -28.538, 19.787, -14.075, -29.507, 20.246, -13.044");
	slices.push_back("1cwc, 1cwc, 12.259, 42.305, 30.71, 13.491, 44.37, 32.092, 13.132, 44.63, 33.528");
	slices.push_back("1cwf, 1cwf, 12.372, 42.472, 30.961, 13.411, 44.558, 32.008, 13.158, 45.01, 33.396");


}
void calculateTheoreticalElectronDensity()
{
	string pdb_code = "6eex";

	//##############################################################
	cout << " Calculating theoretical electron density\n";
	string infilename = "C:/Dev/Github/ProteinDataFiles/pdb_cif/" + pdb_code + ".cif";
	string outfilename_ed = "C:/Dev/Github/ProteinDataFiles/LeucippusFiles/EDCIF/" + pdb_code + "_ed_theo.cif";
	string outfilename_sf = "C:/Dev/Github/ProteinDataFiles/LeucippusFiles/EDCIF/" + pdb_code + "_sf_theo.cif";

	//Create the protein object
	CifFile* cifCoordinates = new CifFile(infilename);	
	Protein* prot = new Protein(cifCoordinates->LoopElements["_atom_site"]);		
	
	//create the electron density and structure factors together
	ED_SF edsf;
	edsf.proteinToTheoreticalED(pdb_code, prot);
	edsf.EdTheo->print(outfilename_ed);
	edsf.proteinToTheoreticalSF(pdb_code, prot);
	edsf.SfTheo->print(outfilename_sf);

	// A sanity checkt hat it works
	VectorThree mnp = edsf.EdTheo->getCRSFromXYZ(VectorThree(0, 0, 0));
	cout << "\n(0,0,0)=" << "\n" << mnp.A << "," << mnp.B << "," << mnp.C << "\nAnd back=";
	VectorThree xyz = edsf.EdTheo->getXYZFromCRS(mnp);
	cout << xyz.A << "," << xyz.B << "," << xyz.C << "\n";

	mnp = edsf.EdTheo->getCRSFromXYZ(VectorThree(1, 1, 10));
	cout << "\n(1,1,10)=" << "\n" << mnp.A << "," << mnp.B << "," << mnp.C << "\nAnd back=";
	xyz = edsf.EdTheo->getXYZFromCRS(mnp);
	cout << xyz.A << "," << xyz.B << "," << xyz.C << "\n";

}

void convertCcp4ToCif()
{
	string pdb_code = "1ejg";
	
	//##############################################################
	cout << "Converting Ccp4 file to cif file";	
	string infilename = "C:/Dev/Github/ProteinDataFiles/ccp4_data/" + pdb_code + ".ccp4";
	string outfilename = "C:/Dev/Github/ProteinDataFiles/LeucippusFiles/EDCIF/" + pdb_code + "_ed.cif";	
	ElectronDensity ED = ElectronDensity(pdb_code, infilename, "ccp4");
	ED.print(outfilename);	
	
	VectorThree mnp = ED.getCRSFromXYZ(VectorThree(0, 0, 0));
	cout << "\n(0,0,0)=" << "\n" << mnp.A << "," <<  mnp.B << "," << mnp.C << "\nAnd back=";
	VectorThree xyz = ED.getXYZFromCRS(mnp);
	cout << xyz.A << "," << xyz.B << "," << xyz.C << "\n";

	mnp = ED.getCRSFromXYZ(VectorThree(1, 1, 10));
	cout << "\n(1,1,10)=" << "\n" << mnp.A << "," << mnp.B << "," << mnp.C << "\nAnd back="; 
	xyz = ED.getXYZFromCRS(mnp);
	cout << xyz.A << "," << xyz.B << "," << xyz.C << "\n";
	
}

void old()
{
	int choose_function = 2;
	if (choose_function == 0)
	{
		// Diffraction Q 3
		Lattice* lat1 = new RealLattice(100.02, 90.57, 68.33, 90, 104.48, 90);
		Lattice* rec1 = lat1->makeInverseLattice();
		string infoA1 = lat1->printLattice();
		string infoB1 = rec1->printLattice();
		cout << "Answers to Question 3\n";
		cout << infoA1 << endl;
		cout << infoB1 << endl;

		// Diffraction Q 4
		Lattice* lat2 = new RealLattice(45.69, 150.22, 75.87, 90, 90, 90);
		Lattice* rec2 = lat2->makeInverseLattice();
		string infoA2 = lat2->printLattice();
		string infoB2 = rec2->printLattice();
		cout << "Answers to Question 4\n";
		cout << infoA2 << endl;
		cout << infoB2 << endl;

		//Add bragg's law in real and reciprocal space
	}
	else if (choose_function == 2)
	{
		cout << "Attempting to create a molecular replacement object\n";
		string filename_sf = "C:/Dev/Github/ProteinDataFiles/sf_data/1ejg-sf.cif";
		string filename_pdb = "C:/Dev/Github/ProteinDataFiles/sf_data/1ejg.cif";
		string filename_ed_molrep = "C:/Dev/Github/ProteinDataFiles/sf_data/1ejg_ed_molrep_leu.cif";
		//MolRep moll(filename_pdb, filename_sf);		
		//moll.EdMolRep->print(filename_ed_molrep);
		string filename_sf_real = "C:/Dev/Github/ProteinDataFiles/sf_data/1ejg-sf-real.cif";
		CifFile* cf = new CifFile(filename_sf);
		StructureFactorsExperimental* sf = new StructureFactorsExperimental(cf);
		sf->printRealConversion(filename_sf_real);



	}

}
