
#include "StructureFactors.h"

/*******************************************
* EXPERIMENTAL
*******************************************/
StructureFactorsExperimental::StructureFactorsExperimental(string fileName)
{

}
/*******************************************
* THEORETICAL
*******************************************/
StructureFactorsTheoretical::StructureFactorsTheoretical(int h_start, int h_end, int k_start, int k_end, int l_start, int l_end, vector<Atom*> atoms)
{//the model for structure factors
	for (int h = h_start; h <= h_end; ++h)
	{
		for (int k = k_start; k <= k_end; ++k)
		{
			for (int l = l_start; l <= l_end; ++l)
			{
				complex<double> im_sf = 0;
				for (unsigned int a = 0; a < atoms.size(); ++a)
				{
					VectorThree hkl = VectorThree(h, k, l);					
					complex<double> im_sfa = atoms[a]->structureFactorContributionTheoretical(hkl);
					im_sf = im_sf + im_sfa;
				}
				StructureFactor *sf = new StructureFactor(h, k, l, im_sf);
				_sfs.push_back(sf);
				
			}
		}
	}
}

/************************************************************************************************************************************
* Class for a single structure factor
***************************************************************************************************************************************/
StructureFactor::StructureFactor(int h, int k, int l, double intensity)
{
	_h = h;
	_k = k;
	_l = l;
	_intensity = intensity;
	Experimental = true;
}
StructureFactor::StructureFactor(int h, int k, int l, complex<double> sf)
{
	_h = h;
	_k = k;
	_l = l;
	_sf = sf;
	_intensity = 0;
	Experimental = false;
}