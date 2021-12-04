#pragma once
/*
* Class to support crystallography learning
*/

#include <string>
#include <vector>
#include <complex>
#include <map>
#include "Xray.h"

using namespace std;

class Atom
{
private:
    double _x;
    double _y;
    double _z;
    double _bfactor;
    double _occupancy;
    string _atomType;
    string _atomClass;
        
public:
    //public functions
    Atom(double x, double y, double z, double bfactor, double occupancy, string atomType, string atomClass);    
    double calculateScatterFactor(Xray xry, double theta);
    double calculateTempFactor(Xray xry, double theta);
    void addAnisotropicTempFactors(int u11, int u22, int u33, int i12, int u13, int u23);//ANISOU  220  N   ASN A  12      256    227    172     22     -7    -24       N  http://skuld.bmsc.washington.edu/parvati/pdb_anisou.html
    //***** structure factors ****
    complex<double> structureFactorContributionTheoretical(VectorThree hkl);
    //we measure them structure factors directly so no experimental version needed
    // **** electron density ****
    double electronDensityContributionTheoretical(VectorThree xyz, string model="IAM");
    complex<double> electronDensityContributionFromSF(VectorThree hkl, double intensity);
    
private:
    double getIAMDensityInternal(VectorThree ABC, VectorThree XYZ, double occupancy);
    double getDensityComponent(double d, double x, double y);
};