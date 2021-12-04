#pragma once

#include <string>
#include <vector>
#include <map>
#include <complex>
#include "Lattice.h"
#include "CifFile.h"
#include "MatrixThreeThree.h" 

using namespace std;

class Voxel
{

public:
    int Id;
    VectorThree CRS;    
    double v;    
};

class ElectronDensity
{
public:
    VectorThree NumsCRS;
    VectorThree StartsCRS;
    VectorThree NumsXYZ;
    VectorThree LengthABC;
    RealLattice LattFractional;    
    double EDMin;
    double EDMax;
    double EDTotal;
    int EDCount;
private:
    string _name;       
    int _mode;              
    VectorThree _alphaBetaGamma;
    VectorThree _mapCRS;
    VectorThree _mapXYZ;
    MatrixThreeThree _orthoMat;
    MatrixThreeThree _deOrthoMat;
public:
    map<int,Voxel*> Voxels;
    
    ElectronDensity(string name);
    ElectronDensity(string name, string ccp4File, string type);
    void calculateOrthoMats();
    
    VectorThree getCRSFromXYZ(VectorThree vXYZ);
    VectorThree getXYZFromCRS(VectorThree vMNP);

    complex<double> structureFactorContributionTheoretical(VectorThree hkl, VectorThree crs, double ed);

    void addVoxel(int cm, int rn, int sp, double v, int id);
    void print(string filename);
};
