

#include <cmath>


#include "Atom.h"
#include "PeriodicTable.h"



#define M_PI 3.14159265358979323846  /* pi */

Atom::Atom(double x, double y, double z, double bfactor, double occupancy, string AtomType)
{
	_x = x;
	_y = y;
	_z = z;
	_bfactor = bfactor;	
	_occupancy = occupancy;
	_atomType = AtomType;
}



double Atom::calculateScatterFactor(Xray xry, double theta)
{
	//    TODO the atomic scattering factor is dependednt on the theta and lamba where when theta=0 = number of electrons and shrinks to 9
	//  from International Tables for Crystallography (Volume C, 1995) https://it.iucr.org/Cb/contents/
	// how can it be both calculable and in a table???

	return 0; //TODO 
}

void Atom::addAnisotropicTempFactors(int u, int u11, int u22, int, int u33, int i12)
{
	/*
	* http://skuld.bmsc.washington.edu/parvati/pdb_anisou.html scaled by 10^4
	* http://legacy.ccp4.ac.uk/html/pdbformat.html
	* The isotropic temperature factor defined in the ATOM card is defined as:
      Biso = 8pi² × (U(1,1) + U(2,2) + U(3,3))/3

	  https://www.cgl.ucsf.edu/pipermail/chimera-users/2012-February/007246.html
	  
	  The six terms actually describe a 3x3 matrix that is symmetric about  
      the diagonal.  The eigenvectors and eigenvalues of that matrix are the  
      atomic displacement axes and the mean squares of the displacements  
      respectively.  The first thing to know is that the values shown in the  
      ANISOU records are scaled by 10**4.  The other thing to know is that  
      the six numbers correspond to these positions in the matrix:  1,1 2,2  
      3,3 1,2 1,3 2,3.  Since the eigenvalues are mean squares, one  
      typically works with the square roots of the eigenvalues.
	  
	  from numpy.linalg import svd
	  ignore, lengths, axes = svd(atom.anisoU)
      from numpy import sqrt
      lengths2 = sqrt(lengths)



	*/

}

double Atom::calculateTempFactor(Xray xry, double theta)
{	
	double T = exp(-1 * _bfactor * sin(theta) * sin(theta) / (xry.Lamda * xry.Lamda));
	return T;
}

complex<double> Atom::structureFactorContributionTheoretical(VectorThree hkl)
{
	/*
	* F(S) = O.f.T.G where O is occupancy, f is scatter, T is temp fact and G is
	* G = e^(2.pi.i.S.r) is the geometry factor
	* Or, given the bragg planes we end up with
	* F(hjkl) = f . e^(2.pi.i(hx + ky + lz)
	*/
	VectorThree xyz = VectorThree(_x, _y, _z);
	complex<double> im_f = electronDensityContributionTheoretical(hkl);
	complex<double> i;
	i = sqrt(-1);
	complex<double> fhkl = im_f * exp(2 * M_PI * i * hkl.getDotProduct(xyz));
	return fhkl;
}

complex<double> Atom::electronDensityContributionExperimental(VectorThree hkl, double intensity) 
{
	
	VectorThree xyz = VectorThree(_x, _y, _z);
	double f = intensity; // a function of bfactor and temp factor and scatter factor;TODO
	complex<double> i;
	i = sqrt(-1);
	complex<double> im_ed = f * exp(-2.0 * M_PI * i * hkl.getDotProduct(xyz));
	//and divide by V
	return im_ed;
}

double Atom::electronDensityContributionTheoretical(VectorThree xyz, string model)
{//we only have the IAM model at the moment
	return getIAMDensityInternal(VectorThree(_x,_y,_z), xyz, _occupancy);
}

double Atom::getIAMDensityInternal(VectorThree ABC, VectorThree XYZ, double occupancy)
{
	double DISTANCE_CAP = 5;
	/*
		https://www.phenix-online.org/presentations/latest/pavel_maps_2.pdf
		https://github.com/project-gemmi/gemmi/blob/master/include/gemmi/dencalc.hpp
		https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Modules_and_Websites_(Inorganic_Chemistry)/Crystallography/X-rays/CromerMann_coefficients

		rho(r) = sum(i = 1 to 4)
		a(i) * [4 * pi / (bi + B)] ^ 1.5
		* exp[-4 * pi ^ 2 * r ^ 2 / (bi + B)]

		+ c * [4 * pi / B] ^ 1.5
		* exp[-4 * pi ^ 2 * r ^ 2 / B]

	*/
	double distance = XYZ.distance(ABC);

	//let's decide that at a certain distance there is no need to do all thi calculation
	if (distance < DISTANCE_CAP || DISTANCE_CAP == 0)
	{
		double density = 0;
		vector<double> cromerMann = PeriodicTable::getCromerMannCoefficients(_atomType);
		if (cromerMann.size() > 1)
		{

			density += getDensityComponent(distance, cromerMann[0], (cromerMann[4] + _bfactor));
			density += getDensityComponent(distance, cromerMann[1], (cromerMann[5] + _bfactor));
			density += getDensityComponent(distance, cromerMann[2], (cromerMann[6] + _bfactor));
			density += getDensityComponent(distance, cromerMann[3], (cromerMann[7] + _bfactor));
			double c = cromerMann[8];
			c = getDensityComponent(distance, c, _bfactor);
			if (!isnan((double)c))
				density += c;

			return occupancy * density;
		}
		else
			return 0;
	}
	else
	{
		return 0;
	}
}
double Atom::getDensityComponent(double d, double x, double y)
{
	/*
	
	*/
	if (abs(y) < 0.00001)
		return 0;

	double bottom = 4 * M_PI / y;
	double raised = pow(bottom, 3);
	raised = sqrt(raised);

	double index = 4 * pow(M_PI, 2) * pow(d, 2) / y;
	double exponent = exp(index);
	if (abs(exponent) < 0.00001)
		return 0;

	exponent = 1 / exponent;
	return x * raised * exponent;
}