
/************************************************************************
* RSA 28.11.21
************************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <cstring> // memcpy
#include "VectorThree.h"
#include <iomanip>

#include "Ccp4File.h"
#include "Helper.h"


using namespace std;

typedef unsigned char uchar;

Ccp4File::Ccp4File(string ccp4File)
{
          
    loadFile(ccp4File);
    
}

void Ccp4File::loadFile(string ccp4File)
{    
    PI = 3.14159265; 
    
    createWordsData(ccp4File);

    W01_NC = _wordsDataInt[0];
    W02_NR = _wordsDataInt[1];
    W03_NS = _wordsDataInt[2];
    W04_Mode = _wordsDataInt[3];
    W05_NCSTART = _wordsDataInt[4];
    W06_NRSTART = _wordsDataInt[5];
    W07_NSSTART = _wordsDataInt[6];
    W08_NX = _wordsDataInt[7];
    W09_NY = _wordsDataInt[8];
    W10_NZ = _wordsDataInt[9];
    W11_CELLA_X = _wordsDataFloat[10];
    W12_CELLA_Y = _wordsDataFloat[11];
    W13_CELLA_Z = _wordsDataFloat[12];
    W14_CELLB_X = _wordsDataFloat[13];
    W15_CELLB_Y = _wordsDataFloat[14];
    W16_CELLB_Z = _wordsDataFloat[15];
    W17_MAPC = _wordsDataInt[16];
    W18_MAPR = _wordsDataInt[17];
    W19_MAPS = _wordsDataInt[18];
    W17_MAPC -= 1;
    W18_MAPR -= 1;
    W19_MAPS -= 1;
    float w20_DMIN = _wordsDataFloat[19];
    float w21_DMAX = _wordsDataFloat[20];
    _w22_DMEAN = _wordsDataFloat[21];
    
    int len = W01_NC * W02_NR * W03_NS;
    int startBulk = (int)_wordsDataFloat.size() - len;
    int count = 0;
    for (unsigned int i = startBulk; i < _wordsDataFloat.size(); ++i)
    {
        float mtx = _wordsDataFloat[i];
        Matrix.push_back(mtx);
        count++;
        if (i % 500 == 0)
        {
            stringstream strc;
            strc << i << "/" << _wordsDataFloat.size();            
        }
    }
    calculateOrthoMat(W11_CELLA_X, W12_CELLA_Y, W13_CELLA_Z, W14_CELLB_X, W15_CELLB_Y, W16_CELLB_Z);
    calculateOrigin(W05_NCSTART, W06_NRSTART, W07_NSSTART, W17_MAPC, W18_MAPR, W19_MAPS);

    _map2xyz.push_back(0);
    _map2xyz.push_back(0);
    _map2xyz.push_back(0);
    _map2xyz[W17_MAPC] = 0;
    _map2xyz[W18_MAPR] = 1;
    _map2xyz[W19_MAPS] = 2;

    _map2crs.push_back(0);
    _map2crs.push_back(0);
    _map2crs.push_back(0);
    _map2crs[0] = W17_MAPC;
    _map2crs[1] = W18_MAPR;
    _map2crs[2] = W19_MAPS;

    _cellDims.push_back(0.0);
    _cellDims.push_back(0.0);
    _cellDims.push_back(0.0);
    _cellDims[0] = W11_CELLA_X;
    _cellDims[1] = W12_CELLA_Y;
    _cellDims[2] = W13_CELLA_Z;

    _axisSampling.push_back(0);
    _axisSampling.push_back(0);
    _axisSampling.push_back(0);
    _axisSampling[0] = W08_NX;
    _axisSampling[1] = W09_NY;
    _axisSampling[2] = W10_NZ;

    _crsStart.push_back(0);
    _crsStart.push_back(0);
    _crsStart.push_back(0);
    _crsStart[0] = W05_NCSTART;
    _crsStart[1] = W06_NRSTART;
    _crsStart[2] = W07_NSSTART;

    _dimOrder.push_back(0);
    _dimOrder.push_back(0);
    _dimOrder.push_back(0);
    _dimOrder[0] = W01_NC;
    _dimOrder[1] = W02_NR;
    _dimOrder[2] = W03_NS;


}

VectorThree Ccp4File::getCRS(int position)
{
    int sliceSize = W01_NC * W02_NR;
    int i = position / sliceSize;
    int remainder = position % sliceSize;
    int j = remainder / W01_NC;
    int k = remainder % W01_NC;
    VectorThree CRS(i,j,k);
    return CRS;
}

void Ccp4File::createWordsData(string ccp4File)
{
    ifstream infile;        
    infile.open(ccp4File.c_str(), ios::binary | ios::in);
    unsigned char temp[sizeof(float)];
    int count = 0;
    while (infile.read(reinterpret_cast<char*>(temp), sizeof(float)))
    {
        ++count;
        if (count % 25000 == 0)
        {
            stringstream strc;            
            cout << "...word line " << count << "\n";
        }

        string ss(reinterpret_cast<char const*>(temp));
        float sf = reinterpret_cast<float&>(temp);
        int si = reinterpret_cast<int&>(temp);

        if (ss.size() > 8)
            ss = ss.substr(8);


        int pos = ss.find(",");
        if (pos > 0)
        {
            vector<string> sv = Helper::stringToVector(ss, ",");
            ss = "";
            for (unsigned int c = 0; c < sv.size(); ++c)
            {
                ss += sv[c] + ".";
            }
        }

        _wordsDataStr.push_back(ss);
        _wordsDataFloat.push_back(sf);
        _wordsDataInt.push_back(si);
    }
    infile.close();
    int symmetry = _wordsDataInt[23]; //the length of the symmetry data is held here
    int length = (int)_wordsDataInt.size();
    int nCnRnS = _wordsDataInt[0] * _wordsDataInt[1] * _wordsDataInt[2];
    createWordsList(symmetry, length, nCnRnS);
}

void Ccp4File::createWordsList(int symmetry, int length, int nCnRnS)
{//https://ftp.ebi.ac.uk/pub/databases/emdb/doc/Map-format/current/EMDB_map_format.pdf
    //int unless otherwise stated
    _wordsList.push_back("1_NC");
    _wordsList.push_back("2_NR");
    _wordsList.push_back("3_NS");

    _wordsList.push_back("4_MODE");

    _wordsList.push_back("5_NCSTART");
    _wordsList.push_back("6_NRSTART");
    _wordsList.push_back("7_NSSTART");

    _wordsList.push_back("8_NX");
    _wordsList.push_back("9_NY");
    _wordsList.push_back("10_NZ");

    _wordsList.push_back("11_X_LENGTH");//float
    _wordsList.push_back("12_Y_LENGTH");//float
    _wordsList.push_back("13_Z_LENGTH");//float

    _wordsList.push_back("14_ALPHA");//float
    _wordsList.push_back("15_BETA");//float
    _wordsList.push_back("16_GAMMA");//float

    _wordsList.push_back("17_MAPC");
    _wordsList.push_back("18_MAPR");
    _wordsList.push_back("19_MAPS");

    _wordsList.push_back("20_AMIN");//float
    _wordsList.push_back("21_AMAX");//float
    _wordsList.push_back("22_AMEAN");//float

    _wordsList.push_back("23_ISPG");

    _wordsList.push_back("24_NYYMBT");//num of bytes in symmetry table

    _wordsList.push_back("25_LSKFLG");//skew flag

    _wordsList.push_back("26_SKWMAT_S11");//float, skew matrix
    _wordsList.push_back("27_SKWMAT_S12");//float, skew matrix
    _wordsList.push_back("28_SKWMAT_S13");//float, skew matrix
    _wordsList.push_back("29_SKWMAT_S21");//float, skew matrix
    _wordsList.push_back("30_SKWMAT_S22");//float, skew matrix
    _wordsList.push_back("31_SKWMAT_S23");//float, skew matrix
    _wordsList.push_back("32_SKWMAT_S31");//float, skew matrix
    _wordsList.push_back("33_SKWMAT_S32");//float, skew matrix
    _wordsList.push_back("34_SKWMAT_S33");//float, skew matrix

    _wordsList.push_back("35_SKWTRN_T1");//float, skew turn
    _wordsList.push_back("36_SKWTRN_T2");//float, skew turn
    _wordsList.push_back("37_SKWTRN_T3");//float, skew turn

    for (unsigned int i = 38; i < 53; ++i)
    {
        stringstream word;
        word << i << "_EXTRA";//binary
        _wordsList.push_back(word.str());
    }

    _wordsList.push_back("53_MAP");//char MRC or CCP4 I think
    _wordsList.push_back("54_MACHST");//binary machine stamp
    _wordsList.push_back("55_RMS");//float root mean square deviation
    _wordsList.push_back("56_NLABL");// num of labels

    for (unsigned int i = 57; i < 257; ++i)
    {
        stringstream word;
        word << i << "_LABEL";//binary
        _wordsList.push_back(word.str());
    }

    //symmetry info not EDMS XRAY only
    cout << "Create Voxels...\n";
    int startVoxels = length - nCnRnS + 1;
    int symCount = 0;
    for (int i = 257; i < startVoxels; ++i)
    {                
        ++symCount;
        stringstream word;
        word << symCount << "_SYM";//binary
        _wordsList.push_back(word.str());        
    }
    //voxels we have to work out backwards from the data
    int voxCount = 0;
    for (int i = startVoxels; i < length + 1; ++i)
    {
        ++voxCount;
        stringstream word;
        word << voxCount << "_VOXEL";//float
        _wordsList.push_back(word.str());
        if (_wordsList.size() % 25000 == 0)
            cout << "...voxels " << i << "/" << length << "\n";
    }
}


void Ccp4File::calculateOrthoMat(float w11_CELLA_X, float w12_CELLA_Y, float w13_CELLA_Z, float w14_CELLB_X, float w15_CELLB_Y, float w16_CELLB_Z)
{
    // Cell angles is w14_CELLB_X, w15_CELLB_Y, w16_CELLB_Z
    // Cell lengths is w11_CELLA_X , w12_CELLA_Y , w13_CELLA_Z 
    double alpha = PI / 180 * w14_CELLB_X;
    double beta = PI / 180 * w15_CELLB_Y;
    double gamma = PI / 180 * w16_CELLB_Z;
    double temp = sqrt(1 - pow(cos(alpha), 2) - pow(cos(beta), 2) - pow(cos(gamma), 2) + 2 * cos(alpha) * cos(beta) * cos(gamma));

    double v00 = w11_CELLA_X;
    double v01 = w12_CELLA_Y * cos(gamma);
    double v02 = w13_CELLA_Z * cos(beta);
    double v10 = 0;
    double v11 = w12_CELLA_Y * sin(gamma);
    double v12 = w13_CELLA_Z * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
    double v20 = 0;
    double v21 = 0;
    double v22 = w13_CELLA_Z * temp / sin(gamma);

    _orthoMat.putValue(w11_CELLA_X, 0, 0);
    _orthoMat.putValue(w12_CELLA_Y * cos(gamma), 0, 1);
    _orthoMat.putValue(w13_CELLA_Z * cos(beta), 0, 2);
    _orthoMat.putValue(0, 1, 0);
    _orthoMat.putValue(w12_CELLA_Y * sin(gamma), 1, 1);
    _orthoMat.putValue(w13_CELLA_Z * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma), 1, 2);
    _orthoMat.putValue(0, 2, 0);
    _orthoMat.putValue(0, 2, 1);
    _orthoMat.putValue(w13_CELLA_Z * temp / sin(gamma), 2, 2);
    _deOrthoMat = _orthoMat.getInverse();    
}

void Ccp4File::calculateOrigin(int w05_NXSTART, int w06_NYSTART, int w07_NZSTART, int w17_MAPC, int w18_MAPR, int w19_MAPS)
{
    /****************************
    * These comments are from my C# version and I have no idea currently what they mean (RSA 6/9/21)
    * ******************************
     *TODO I am ignoring the possibility of passing in the origin for nowand using the dot product calc for non orthoganality.
     *The origin is perhaps used for cryoEM only and requires orthoganility
     *CRSSTART is w05_NXSTART, w06_NYSTART, w07_NZSTART
     *Cell dims w08_MX, w09_MY, w10_MZ;
     *Map of indices from crs to xyz is w17_MAPC, w18_MAPR, w19_MAPS
     */

    VectorThree oro;

    for (int i = 0; i < 3; ++i)
    {
        int startVal = 0;
        if (w17_MAPC == i)
            startVal = w05_NXSTART;
        else if (w18_MAPR == i)
            startVal = w06_NYSTART;
        else
            startVal = w07_NZSTART;

        oro.putByIndex(i, startVal);
    }
    oro.putByIndex(0, oro.getByIndex(0) / W08_NX);
    oro.putByIndex(1, oro.getByIndex(1) / W09_NY);
    oro.putByIndex(2, oro.getByIndex(2) / W10_NZ);
    _origin = _orthoMat.multiply(oro, true);
}