
/*
* Rachel Alcraft 20/10/2021
*/

#include <sstream>
#include "Lattice.h"

#define M_PI 3.14159265358979323846  /* pi */

Lattice::Lattice()
{
	_a = 0;
	_b = 0;
	_c = 0;	
	_alpha = 0;
	_gamma = 0;	
	_beta = 0;
}

Lattice::Lattice(double a, double b, double c, double alpha, double beta, double gamma)
{
	_a = a;
	_b = b;
	_c = c;
	//alpha and gamma must = 90
	_alpha = alpha;	
	_gamma = gamma;
	//if beta also = 90 the cell is orthorhombic
	// if beta > 90 (it can't be smaller) then it is monoclinic
	_beta = beta;
	_v = 0;
}

string Lattice::printLattice()
{
	stringstream info;
	info << "a=" << _a << " b=" << _b << " c=" << _c << " alpha=" << _alpha << " beta=" << _beta << " gamma=" << _gamma << " volume=" << _v << " kind=" << _latticeKind;
	return info.str();	
}

//**** OVERRIDE FOR REAL LATTICE************************************************************************************************************************************
RealLattice::RealLattice(double a, double b, double c, double alpha, double beta, double gamma) :Lattice(a, b, c, alpha, beta, gamma)
{
	if (_alpha == 90 && _gamma == 90 && _beta == 90)
	{
		_v = _a * _b * _c;
		_latticeKind = "orthorhombic";
	}
	else if (beta > 90)
	{
		_v = _a * _b * _c * sin(_beta * M_PI / 180);//convert to radians for sine functions
		_latticeKind = "monoclinic";
	}	
}
Lattice* RealLattice::makeInverseLattice()
{
	if (_alpha == 90 && _gamma == 90 && _beta == 90)
	{
		double ap = 1 / _a;
		double bp = 1 / _b;
		double cp = 1 / _c;
		double alphap = 90;
		double betap = 90;
		double gammap = 90;
		return new ReciprocalLattice(ap, bp, cp, alphap, betap, gammap);
	}
	else if (_beta > 90)
	{
		double ap = 1 / (_a * sin(_beta * M_PI / 180));
		double bp = 1 / _b;
		double cp = 1 / (_c * sin(_beta * M_PI / 180));
		double alphap = 90;
		double betap = 180 - _beta;
		double gammap = 90;
		return new ReciprocalLattice(ap, bp, cp, alphap, betap, gammap);
	}
}
//************** Override for RECIPROCAL LATTICE **********************************************************************************************************
ReciprocalLattice::ReciprocalLattice(double a, double b, double c, double alpha, double beta, double gamma):Lattice(a,b,c,alpha,beta,gamma)
{	
	if (_alpha == 90 && _gamma == 90 && _beta == 90)
	{
		_v = _a * _b * _c;
		_latticeKind = "reciprocal orthorhombic";
	}
	else if (beta < 90)
	{
		_v = _a * _b * _c * sin(_beta * M_PI / 180);//convert to radians for sine functions
		_latticeKind = "reciprocal monoclinic";
	}
}
Lattice* ReciprocalLattice::makeInverseLattice()
{

	if (_alpha == 90 && _gamma == 90 && _beta == 90)
	{
		double ap = 1 / _a;
		double bp = 1 / _b;
		double cp = 1 / _c;
		double alphap = 90;
		double betap = 90;
		double gammap = 90;
		return new RealLattice(ap, bp, cp, alphap, betap, gammap);
	}
	else if (_beta < 90)
	{
		double ap = 1 / (_a * sin(_beta * M_PI / 180));
		double bp = 1 / _b;
		double cp = 1 / (_c * sin(_beta * M_PI / 180));
		double alphap = 90;
		double betap = 180 - _beta;
		double gammap = 90;
		return new RealLattice(ap, bp, cp, alphap, betap, gammap);
	}
}
//*******************************************************************************************************