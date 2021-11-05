#pragma once
/*
* Rachel Alcraft, 31/10/2021
* Class written as training excercise in crystallography course
*/

#include <string>
#include <vector>
#include <complex>

#include "VectorThree.h"
#include "Atom.h"
#include "Lattice.h"

using namespace std;
class StructureFactor
{
public://Lazy public interface
    bool Experimental;
protected:
    complex<double> _sf;//from model we get complex numbers
    double _intensity; //from experiment we only get the real intensities
    int _h;
    int _k;
    int _l;
public:
    //Functions

    StructureFactor(int h, int k, int l, double intensity);
    StructureFactor(int h, int k, int l, complex<double> sf);
};


class StructureFactors
{    
protected:
    vector<StructureFactor*> _sfs;//from model we get complex numbers
    RealLattice _lattice;    
public:
    //Functions
    StructureFactors() {};
};
/********************************************************
* Create EXPERIMENTAL Structure Factors
*********************************************************/
class StructureFactorsExperimental : public StructureFactors
{
private:
    double _waveLength;

public:
    StructureFactorsExperimental(string fileName);    
};
/********************************************************
* Create THEORETICAL Structure Factors
*********************************************************/
class StructureFactorsTheoretical : public StructureFactors
{
public:
    StructureFactorsTheoretical(int h_start, int h_end, int k_start, int k_end, int l_start, int l_end, vector<Atom*> atoms);
};




