#pragma once

/************************************************************************
* RSA 28/11/21
* https://ftp.ebi.ac.uk/pub/databases/emdb/doc/Map-format/current/EMDB_map_format.pdf
************************************************************************/

#include <string>
#include <vector>
#include <map>
#include "MatrixThreeThree.h"

using namespace std;

class Ccp4File
{
private:

public:	
	int W04_Mode;
	int W01_NC;
	int W02_NR;
	int W03_NS;
	int W05_NCSTART;
	int W06_NRSTART;
	int W07_NSSTART;	
	int W08_NX;
	int W09_NY;
	int W10_NZ;	
	float W11_CELLA_X;
	float W12_CELLA_Y;
	float W13_CELLA_Z;
	float W14_CELLB_X;
	float W15_CELLB_Y;
	float W16_CELLB_Z;
	int W17_MAPC;
	int W18_MAPR;
	int W19_MAPS;
	
private:
	vector<string> _wordsList;
	vector<string> _wordsDataStr;
	vector<int> _wordsDataInt;
	vector<float> _wordsDataFloat;
	double PI;
	

	float _w22_DMEAN;

	//Calculation data
	MatrixThreeThree _orthoMat;
	MatrixThreeThree _deOrthoMat;
	vector<int> _map2xyz;
	vector<int> _map2crs;
	vector<float>_cellDims;
	vector<int> _axisSampling;
	vector<int> _crsStart;
	vector<int> _dimOrder;
	VectorThree _origin;


	//Helper functioms	
	void calculateOrthoMat(float w11_CELLA_X, float w12_CELLA_Y, float w13_CELLA_Z, float w14_CELLB_X, float w15_CELLB_Y, float w16_CELLB_Z);
	void calculateOrigin(int w05_NXSTART, int w06_NYSTART, int w07_NZSTART, int w17_MAPC, int w18_MAPR, int w19_MAPS);
	
public:
	//The matrix data lazily as a public accessor
	vector<float> Matrix;	

public:
	Ccp4File(string ccp4File);
	VectorThree Ccp4File::getCRS(int position);
	
private:
	void loadFile(string ccp4File);
	void createWordsData(string ccp4File);
	void createWordsList(int symmetry, int length, int nCnRnS);
};



