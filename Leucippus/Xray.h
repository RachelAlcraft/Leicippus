#pragma once
/*
* Rachel Alcraft, 31/10/2021
* Class written as training excercise in crystallography course
*/

#include <string>
#include <vector>
#include "VectorThree.h"

using namespace std;

class Xray
{
public://Lazy public interface
    double Lamda;

protected:    
    double _amplitude;
    double _phase;
    VectorThree _direction;


public:
    //Functions
    Xray(double amplitude, double phase, double lamda, VectorThree direction);
    double displacement(double x);   
    Xray scatter(VectorThree position, double theta, double electron_density);//but how do you know the angle??
};



