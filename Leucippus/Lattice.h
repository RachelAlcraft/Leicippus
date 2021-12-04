#pragma once
/*
* Rachel Alcraft, 20/10/2021
* Class written as training excercise in crystallography course
*/

#include <string>
#include <vector>
#include "VectorThree.h"
//#include "ReciprocalLattice.h"

using namespace std;

class Lattice
{
public:
    //cell dimensions
    double a;
    double b;
    double c;
    //cell angles
    double alpha;
    double beta;
    double gamma;

protected:
    //cell type
    string _latticeKind;

    // cell volume
    double _v;
        
public:    
    //Functions
    Lattice();
    Lattice(double a, double b, double c, double alpha, double beta, double gamma);
    string printLattice();
    virtual Lattice* makeInverseLattice() = 0;    
};

//***************************** Main LATTICE ********************************************************
class RealLattice :public Lattice
{
public:
    //Functions
    RealLattice() {}
    RealLattice(double a, double b, double c, double alpha, double beta, double gamma);
    Lattice* makeInverseLattice() override;
    VectorThree CarteToOrtho(VectorThree xyz);
    VectorThree OrthoToCarte(VectorThree mnp);
    VectorThree FracToRecip(VectorThree mnp);
};
//***************************** RECIPROCAL LATTICE ********************************************************
class ReciprocalLattice :public Lattice
{
public:
    //Functions
    ReciprocalLattice() {}
    ReciprocalLattice(double a, double b, double c, double alpha, double beta, double gamma);    
    Lattice* makeInverseLattice() override;
};


