
#include <cmath>
#include <sstream>
#include <iomanip>

#include "VectorThree.h"

using namespace std;

VectorThree::VectorThree()
{
    A = 0;
    B = 0;
    C = 0;
}

VectorThree::VectorThree(double a, double b, double c)
{
    A = a;
    B = b;
    C = c;    
}

double VectorThree::getByIndex(int idx)
{
    if (idx == 0)
        return A;
    else if (idx == 1)
        return B;
    else // (idx == 2)
        return C;
}

void VectorThree::putByIndex(int idx, double val)
{
    if (idx == 0)
        A = val;
    else if (idx == 1)
        B = val;
    else // (idx == 2)
        C = val;
}

double VectorThree::distance(VectorThree ABC)
{
    double dis = (A - ABC.A) * (A - ABC.A) + (B - ABC.B) * (B - ABC.B) + (C - ABC.C) * (C - ABC.C);
    return sqrt(dis);
}

double VectorThree::getMagnitude()
{
    double mag = (A * A) + (B * B) + (C * C);
    return sqrt(mag);
}
VectorThree VectorThree::operator+(VectorThree const& obj)
{
    A += obj.A;
    B += obj.B;
    C += obj.C;
    return VectorThree(A, B, C);
}
VectorThree VectorThree::operator-(VectorThree const& obj)
{
    A -= obj.A;
    B -= obj.B;
    C -= obj.C;
    return VectorThree(A, B, C);
}

VectorThree VectorThree::operator/(double val)
{
    A /= val;
    B /= val;
    C /= val;
    return VectorThree(A, B, C);
}
double VectorThree::getAngle(VectorThree vec)
{
    VectorThree BA(0 - A, 0 - B, 0 - C);
    VectorThree BC(0 - vec.A, 0 - vec.B, 0 - vec.C);
    double dot = BA.getDotProduct(BC);
    double magBA = BA.getMagnitude();
    double magBC = BC.getMagnitude();
    double cosTheta = dot / (magBA * magBC);
    double theta = acos(cosTheta);
    return theta; //in radians
}

VectorThree VectorThree::getMillerIndices(int m, int n, int p)
{//find the miller indizes from the a,b,c, lattice points in real space, which are m,n,p intgers of a,b,c    
    int lcm = getLCM(m,n,p);
    
    int h = 0;
    int k = 0;
    int l = 0;

    if (m > 0)
        h = lcm / m;
    if (n > 0)
        k = lcm / n;
    if (p > 0)
        l = lcm / p;

    //reflection indices may be multiplesof these 
    return VectorThree(h,k,l);
}
VectorThree VectorThree::getCrystalIndices(int h, int k, int l)
{    
    int lcm = getLCM(h,k,l);
    
    int m = 0;
    int n = 0;
    int p = 0;

    if (h > 0)
        m = lcm / h;
    if (k > 0)
        n= lcm / k;
    if (l > 0)
        p = lcm / l;

    //reflection indices may be multiplesof these 
    return VectorThree(h, k, l);
}

int VectorThree::getLCM(int a, int b, int c)
{//if any of the inputs are -1 then that means infinity
    double lcm = 1;
    if (a > 0)
        lcm *= a;
    if (b > 0)
        lcm *= b;
    if (c > 0)
        lcm *= c;

    return lcm;
}
double VectorThree::getDotProduct(VectorThree vec)
{
    double px = A * vec.A;
    double py = B * vec.B;
    double pz = C * vec.C;
    return px + py + pz;
}


