#pragma once
/*
* Rachel Alcraft, 20/10/2021
* Class written as training excercise in crystallography course
*/

#include <string>
#include <vector>
//#include "ReciprocalLattice.h"

using namespace std;

class Lattice
{
protected:
    //cell type
    string _latticeKind;
    //cell dimensions
    double _a; 
    double _b;
    double _c;
    //cell angles
    double _alpha;
    double _beta;
    double _gamma;
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
    RealLattice(double a, double b, double c, double alpha, double beta, double gamma);
    Lattice* makeInverseLattice() override;
};
//***************************** RECIPROCAL LATTICE ********************************************************
class ReciprocalLattice :public Lattice
{
public:
    //Functions
    ReciprocalLattice(double a, double b, double c, double alpha, double beta, double gamma);    
    Lattice* makeInverseLattice() override;
};


