#pragma once

#include <string>
#include <vector>
#include <map>

#include "Atom.h"

using namespace std;

class Protein
{
public:
    vector<Atom*> Atoms;
    double MinX;
    double MaxX;
    double MinY;
    double MaxY;
    double MinZ;
    double MaxZ;


public:
    //public functions    
    Protein(map<string, vector<string> > atomList);
};