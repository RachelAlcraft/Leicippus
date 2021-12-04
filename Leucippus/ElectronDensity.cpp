#include "ElectronDensity.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "Ccp4File.h"

#define M_PI 3.14159265358979323846  /* pi */
typedef complex<double> dcomp;

ElectronDensity::ElectronDensity(string name)
{
	EDMin = 10000;
	EDMax = -10000;
	EDTotal = 0;
	EDCount = 0;
	_name = name;
	_mode = 2;
	NumsCRS = VectorThree(0,0,0);
	StartsCRS = VectorThree(0,0,0);
	NumsXYZ = VectorThree(0,0,0);
	LengthABC = VectorThree(0,0,0);
	_alphaBetaGamma = VectorThree(90, 90, 90);
	_mapCRS = VectorThree(0,1,2);
	_mapXYZ = VectorThree(0, 1, 2);
	
}

ElectronDensity::ElectronDensity(string name, string ccp4File, string  type)
{
	EDMin = 10000;
	EDMax = -10000;
	EDTotal = 0;
	EDCount = 0;
	_name = name;
	if (type == "ccp4")
	{		
		Ccp4File ccp4 = Ccp4File(ccp4File);
		_mode = ccp4.W04_Mode;
		NumsCRS = VectorThree(ccp4.W01_NC, ccp4.W02_NR, ccp4.W03_NS);
		StartsCRS = VectorThree(ccp4.W05_NCSTART, ccp4.W06_NRSTART, ccp4.W07_NSSTART);		
		NumsXYZ = VectorThree(ccp4.W08_NX, ccp4.W09_NY, ccp4.W10_NZ);
		LengthABC = VectorThree(ccp4.W11_CELLA_X, ccp4.W12_CELLA_Y, ccp4.W13_CELLA_Z);
		_alphaBetaGamma = VectorThree(ccp4.W14_CELLB_X, ccp4.W15_CELLB_Y, ccp4.W16_CELLB_Z);
		_mapCRS = VectorThree(ccp4.W17_MAPC, ccp4.W18_MAPR, ccp4.W19_MAPS);		
		_mapXYZ.putByIndex(ccp4.W17_MAPC, 0);
		_mapXYZ.putByIndex(ccp4.W18_MAPR, 1);
		_mapXYZ.putByIndex(ccp4.W19_MAPS, 2);

						
		for (unsigned int i = 0; i < ccp4.Matrix.size(); ++i)
		{
			VectorThree mnp = ccp4.getCRS(i); 
			addVoxel(mnp.A, mnp.B, mnp.C, ccp4.Matrix[i],i);
			EDMin = min(EDMin, double(ccp4.Matrix[i]));
			EDMax = max(EDMax, double(ccp4.Matrix[i]));
			EDTotal += double(ccp4.Matrix[i]);
			EDCount += 1;
		}
		
	}
	else if (type == "cif")
	{
		
	}
	calculateOrthoMats();
	
}
void ElectronDensity::calculateOrthoMats()
{
	//Now make tranformation matrices
	double PI = 3.14159265;
	double alpha = PI / 180 * _alphaBetaGamma.A;
	double beta = PI / 180 * _alphaBetaGamma.B;
	double gamma = PI / 180 * _alphaBetaGamma.C;
	double temp = sqrt(1 - pow(cos(alpha), 2) - pow(cos(beta), 2) - pow(cos(gamma), 2) + 2 * cos(alpha) * cos(beta) * cos(gamma));
	_orthoMat.putValue(LengthABC.A, 0, 0);
	_orthoMat.putValue(LengthABC.B * cos(gamma), 0, 1);
	_orthoMat.putValue(LengthABC.C * cos(beta), 0, 2);
	_orthoMat.putValue(0, 1, 0);
	_orthoMat.putValue(LengthABC.B * sin(gamma), 1, 1);
	_orthoMat.putValue(LengthABC.C * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma), 1, 2);
	_orthoMat.putValue(0, 2, 0);
	_orthoMat.putValue(0, 2, 1);
	_orthoMat.putValue(LengthABC.C * temp / sin(gamma), 2, 2);
	_deOrthoMat = _orthoMat.getInverse();
}

void ElectronDensity::addVoxel(int cm, int rn, int sp, double vv, int id)
{
	Voxel* vox = new Voxel();
	vox->Id = id;
	vox->CRS = VectorThree(cm,rn,sp);	
	vox->v = vv;
	EDMin = min(EDMin, vv);
	EDMax = max(EDMax, vv);
	EDTotal += vv;
	EDCount += 1;
	
	map<int, Voxel*>::iterator iter = Voxels.find(id);
	if (iter == Voxels.end())
	{
		Voxels.insert(pair<int,Voxel*>(id,vox));
	}
	else
	{
		iter->second->v += vv;
	}
}

void ElectronDensity::print(string filename)
{
	cout << "\nPrinting to cif...\n";
	ofstream edfile;
	edfile.open(filename.c_str());
	edfile << "data_LeucipPlus_ED_" << _name << "\n";
	edfile << "#\n";
	edfile << "_transform.mode " << _mode << "\n";
	edfile << "_transform.originX " << "0.00" << "\n";
	edfile << "_transform.originY " << "0.00" << "\n";
	edfile << "_transform.originZ " << "0.00" << "\n";	
	edfile << "_transform.numC " << NumsCRS.A << "\n";
	edfile << "_transform.numR " << NumsCRS.B << "\n";
	edfile << "_transform.numS " << NumsCRS.C << "\n";	
	edfile << "_transform.mapC " << _mapCRS.A << "\n";
	edfile << "_transform.mapR " << _mapCRS.B << "\n";
	edfile << "_transform.mapS " << _mapCRS.C << "\n";
	edfile << "_transform.startC " << StartsCRS.A << "\n";
	edfile << "_transform.startR " << StartsCRS.B << "\n";
	edfile << "_transform.startS " << StartsCRS.C << "\n";
	edfile << "_transform.numX " << NumsXYZ.A << "\n";
	edfile << "_transform.numY " << NumsXYZ.B << "\n";
	edfile << "_transform.numZ " << NumsXYZ.C << "\n";
	edfile << "_transform.lengthX " << LengthABC.A << "\n";
	edfile << "_transform.lengthY " << LengthABC.B << "\n";
	edfile << "_transform.lengthZ " << LengthABC.C << "\n";
	edfile << "_transform.alpha " << _alphaBetaGamma.A << "\n";
	edfile << "_transform.beta " << _alphaBetaGamma.B << "\n";
	edfile << "_transform.gamma " << _alphaBetaGamma.C << "\n";			
	edfile << "#\n";
	edfile << "_stats.min " << EDMin << "\n";
	edfile << "_stats.max " << EDMax << "\n";
	edfile << "_stats.mean " << EDTotal/EDCount << "\n";
	edfile << "#\n";
	edfile << "loop_\n";
	//edfile << "_ed.id\n";
	edfile << "_ed.index_c\n";
	edfile << "_ed.index_r\n";
	edfile << "_ed.index_s\n";
	edfile << "_ed.voxel\n";
	int count = 0;
	for (map<int, Voxel*>::iterator iter = Voxels.begin(); iter!=Voxels.end(); ++iter)
	{
		if (count % 25000 == 0)
		{
			cout << "...Printing " << count << "/" << Voxels.size() << "\n";
		}
		++count;
		stringstream out;
		out << setprecision(5);
		//out << iter->first << " ";
		if (abs(iter->second->v) > 0.00001)
			out << iter->second->CRS.A << " " << iter->second->CRS.B << " " << iter->second->CRS.C << " " << iter->second->v << "\n";
		edfile << out.str();
	}	
	edfile.close();
}

VectorThree ElectronDensity::getCRSFromXYZ(VectorThree vXYZ)
{	
	if (_alphaBetaGamma.A == 90 && _alphaBetaGamma.B == 90 && _alphaBetaGamma.C == 90)
	{
		double xCell = (float(NumsXYZ.A) - 1.0) / float(LengthABC.A);
		double yCell = (float(NumsXYZ.B) - 1.0) / float(LengthABC.B);
		double zCell = (float(NumsXYZ.C) - 1.0) / float(LengthABC.C);
		
		VectorThree shiftedxyz;
		shiftedxyz.A = vXYZ.A - StartsCRS.getByIndex(int(_mapXYZ.A));
		shiftedxyz.B = vXYZ.B - StartsCRS.getByIndex(int(_mapXYZ.B));
		shiftedxyz.C = vXYZ.C - StartsCRS.getByIndex(int(_mapXYZ.C));

		VectorThree stretchedxyz;
		stretchedxyz.A = shiftedxyz.A * xCell;
		stretchedxyz.B = shiftedxyz.B * yCell;
		stretchedxyz.C = shiftedxyz.C * zCell;
		
		VectorThree crs(0, 0, 0);
		crs.A = stretchedxyz.getByIndex(int(_mapCRS.A));
		crs.B = stretchedxyz.getByIndex(int(_mapCRS.B));
		crs.C = stretchedxyz.getByIndex(int(_mapCRS.C));
		return crs;
	}
	else
	{
		VectorThree vFraction = _deOrthoMat.multiply(vXYZ, true);
		vFraction.A = (vFraction.A * NumsXYZ.A) - StartsCRS.getByIndex(int(_mapXYZ.A));
		vFraction.B = (vFraction.B * NumsXYZ.B) - StartsCRS.getByIndex(int(_mapXYZ.B));
		vFraction.C = (vFraction.C * NumsXYZ.C) - StartsCRS.getByIndex(int(_mapXYZ.C));
		VectorThree crs(0, 0, 0);
		crs.A = vFraction.getByIndex(int(_mapCRS.A));
		crs.B = vFraction.getByIndex(int(_mapCRS.B));
		crs.C = vFraction.getByIndex(int(_mapCRS.C));
		return crs;
	}
	
}
VectorThree ElectronDensity::getXYZFromCRS(VectorThree vCRS)
{
	if (_alphaBetaGamma.A == 90 && _alphaBetaGamma.B == 90 && _alphaBetaGamma.C == 90)
	{
		double xCell = (float(NumsXYZ.A) - 1.0) / float(LengthABC.A);
		double yCell = (float(NumsXYZ.B) - 1.0) / float(LengthABC.B);
		double zCell = (float(NumsXYZ.C) - 1.0) / float(LengthABC.C);

		VectorThree movedCRS;
		movedCRS.A = vCRS.getByIndex(int(_mapXYZ.A));
		movedCRS.B = vCRS.getByIndex(int(_mapXYZ.B));
		movedCRS.C = vCRS.getByIndex(int(_mapXYZ.C));

		VectorThree stretched;
		stretched.A = movedCRS.A / xCell;
		stretched.B = movedCRS.B / yCell;
		stretched.C = movedCRS.C / zCell;

		VectorThree shifted;
		shifted.A = stretched.A + StartsCRS.getByIndex(int(_mapCRS.A));
		shifted.B = stretched.B + StartsCRS.getByIndex(int(_mapCRS.B));
		shifted.C = stretched.C + StartsCRS.getByIndex(int(_mapCRS.C));
				
		return shifted;
	}
	else
	{
		VectorThree v1(0, 0, 0);
		v1.A = vCRS.A + StartsCRS.A;
		v1.B = vCRS.B + StartsCRS.B;
		v1.C = vCRS.C + StartsCRS.C;

		VectorThree v2(0, 0, 0);
		v2.putByIndex(int(_mapCRS.A), v1.A);
		v2.putByIndex(int(_mapCRS.B), v1.B);
		v2.putByIndex(int(_mapCRS.C), v1.C);

		VectorThree v3(0, 0, 0);
		v3.A = v2.A / NumsXYZ.A;
		v3.B = v2.B / NumsXYZ.B;
		v3.C = v2.C / NumsXYZ.C;

		VectorThree ret = _orthoMat.multiply(v3, false);
		return ret;
	}
}

complex<double> ElectronDensity::structureFactorContributionTheoretical(VectorThree hkl, VectorThree crs, double ed)
{
	/*
	* F(S) = O.f.T.G where O is occupancy, f is scatter, T is temp fact and G is
	* G = e^(2.pi.i.S.r) is the geometry factor
	* Or, given the bragg planes we end up with
	* F(hjkl) = f . e^(2.pi.i(hx + ky + lz)
	*/	
	
	
	/*dcomp i;
	dcomp a;
	double pi;
	pi = 2 * asin(1);
	i = -1;
	i = sqrt(i);
	a = exp(2 * pi * i);
	cout << "i is " << i << "and Euler was right: e(i pi) = " << a << endl;*/

	dcomp i;	
	i = -1;
	i = sqrt(i);	
	complex<double> fhkl;
	
	if (ed > 0)
	{
		double dot_p = hkl.getDotProduct(crs);				
		fhkl = ed * exp(2 * M_PI * i * dot_p);
	}
	return fhkl;
}
