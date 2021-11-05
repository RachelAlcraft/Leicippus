#pragma once
#include <string>
#include <vector>

using namespace std;

class PeriodicTable
{
public:
    static vector<double> getCromerMannCoefficients(string atomType);
    static string getElement(string atomType);
};