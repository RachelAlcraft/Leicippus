#pragma once
/*
* Rachel Alcraft, 31/10/2021
* Class written as training excercise in crystallography course
*/

#include <string>
#include <vector>
#include <complex>
#include <map>

#include "VectorThree.h"
#include "Atom.h"
#include "Lattice.h"
#include "CifFile.h"

using namespace std;
class StructureFactor
{
public://Lazy public interface
    bool Experimental;
    VectorThree HKL;    
    double Intensity; //from experiment we only get the real intensities    
    complex<double> ComplexSF;//from model we get complex numbers
protected:
    
        
public:
    //Functions

    StructureFactor(VectorThree hkl, double intensity);
    StructureFactor(VectorThree hkl, complex<double> sf);
};


class StructureFactors
{    
public:
    string Name;
    map<string,StructureFactor*> SFs;//from model we get complex numbers
    int hStart;
    int hEnd;
    int kStart;
    int kEnd;
    int lStart;
    int lEnd;
protected:
    RealLattice _lattice;    
public:
    //Functions
    StructureFactors();
    void addStructureFactor(StructureFactor* sf);
    void print(string filename);
};
/********************************************************
* Create EXPERIMENTAL Structure Factors
*********************************************************/
class StructureFactorsExperimental : public StructureFactors
{
private:
    double _waveLength;
    CifFile* _cif;

public:
    StructureFactorsExperimental(CifFile* cif);    
    void printRealConversion(string filename);
    RealLattice getLattice();
};
/********************************************************
* Create THEORETICAL Structure Factors
*********************************************************/
class StructureFactorsTheoretical : public StructureFactors
{
public:
    StructureFactorsTheoretical(string name,int h_start, int h_end, int k_start, int k_end, int l_start, int l_end, vector<Atom*> atoms);
    StructureFactorsTheoretical(string name);
    















};




