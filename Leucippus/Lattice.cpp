
/*
* Rachel Alcraft 20/10/2021
*/

#include <sstream>
#include "Lattice.h"

#define M_PI 3.14159265358979323846  /* pi */

Lattice::Lattice()
{
	a = 1;
	b = 1;
	c = 1;	
	_v = 0;
	alpha = 90;
	gamma = 90;	
	beta = 90;
}

Lattice::Lattice(double aa, double bb, double cc, double aalpha, double bbeta, double ggamma)
{
	a = aa;
	b = bb;
	c = cc;
	//alpha and gamma must = 90
	alpha = aalpha;	
	gamma = ggamma;
	//if beta also = 90 the cell is orthorhombic
	// if beta > 90 (it can't be smaller) then it is monoclinic
	beta = bbeta;
	_v = 0;
}

string Lattice::printLattice()
{
	stringstream info;
	info << "a=" << a << " b=" << b << " c=" << c << " alpha=" << alpha << " beta=" << beta << " gamma=" << gamma << " volume=" << _v << " kind=" << _latticeKind;
	return info.str();	
}
//**** OVERRIDE FOR REAL LATTICE************************************************************************************************************************************
RealLattice::RealLattice(double a, double b, double c, double alpha, double beta, double gamma) :Lattice(a, b, c, alpha, beta, gamma)
{
	if (alpha == 90 && gamma == 90 && beta == 90)
	{
		_v = a * b * c;
		_latticeKind = "orthorhombic";
	}
	else if (beta > 90)
	{
		_v = a * b * c * sin(beta * M_PI / 180);//convert to radians for sine functions
		_latticeKind = "monoclinic";
	}	
}
Lattice* RealLattice::makeInverseLattice()
{
	if (alpha == 90 && gamma == 90 && beta == 90)
	{
		double ap = 1 / a;
		double bp = 1 / b;
		double cp = 1 / c;
		double alphap = 90;
		double betap = 90;
		double gammap = 90;
		return new ReciprocalLattice(ap, bp, cp, alphap, betap, gammap);
	}
	else if (beta > 90)
	{
		double ap = 1 / (a * sin(beta * M_PI / 180));
		double bp = 1 / b;
		double cp = 1 / (c * sin(beta * M_PI / 180));
		double alphap = 90;
		double betap = 180 - beta;
		double gammap = 90;
		return new ReciprocalLattice(ap, bp, cp, alphap, betap, gammap);
	}
	else
		return new ReciprocalLattice();
}

/* Conversion between orthogonal and fractional coordinates
https://ipfs.io/ipfs/QmXoypizjW3WknFiJnKLwHCnL72vedxjQkDDP1mXWo6uco/wiki/Fractional_coordinates.html
*/

VectorThree RealLattice::CarteToOrtho(VectorThree xyz)
{
	return VectorThree();
}

VectorThree RealLattice::OrthoToCarte(VectorThree mnp)
{
	double xFrac = mnp.A * (1 / a) - mnp.B * cos(gamma) / (a * sin(gamma) + mnp.C * b * c * (cos(alpha) * cos(gamma) - cos(beta) / _v * sin(gamma)));
	double yFrac = mnp.B * (1 / (b * sin(gamma))) + mnp.C* a * c * ((cos(beta) * cos(gamma) - cos(alpha))) / (_v * sin(gamma));
	double zFrac = mnp.C * a * b * sin(gamma) / _v;
	return VectorThree(xFrac,yFrac,zFrac);
}
VectorThree RealLattice::FracToRecip(VectorThree mnp)
{
	return VectorThree(1 / mnp.A, 1 / mnp.B, 1 / mnp.C);
}

//************** Override for RECIPROCAL LATTICE **********************************************************************************************************
ReciprocalLattice::ReciprocalLattice(double a, double b, double c, double alpha, double beta, double gamma):Lattice(a,b,c,alpha,beta,gamma)
{	
	if (alpha == 90 && gamma == 90 && beta == 90)
	{
		_v = a * b * c;
		_latticeKind = "reciprocal orthorhombic";
	}
	else if (beta < 90)
	{
		_v = a * b * c * sin(beta * M_PI / 180);//convert to radians for sine functions
		_latticeKind = "reciprocal monoclinic";
	}
}
Lattice* ReciprocalLattice::makeInverseLattice()
{

	if (alpha == 90 && gamma == 90 && beta == 90)
	{
		double ap = 1 / a;
		double bp = 1 / b;
		double cp = 1 / c;
		double alphap = 90;
		double betap = 90;
		double gammap = 90;
		return new RealLattice(ap, bp, cp, alphap, betap, gammap);
	}
	else if (beta < 90)
	{
		double ap = 1 / (a * sin(beta * M_PI / 180));
		double bp = 1 / b;
		double cp = 1 / (c * sin(beta * M_PI / 180));
		double alphap = 90;
		double betap = 180 - beta;
		double gammap = 90;
		return new RealLattice(ap, bp, cp, alphap, betap, gammap);
	}
	else
		return new RealLattice();
}
//*******************************************************************************************************