#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

#include "vapor/MatWaveWavedec.h"
using namespace VAPoR;


const int NX = 1;
const int NY = 1;
bool KeepAppCoeff = false;
float cratio = 16.0;


void test1d(string wavename, const float *srcarr, float *dstarr, int NZ )
{


	//
	// Now do 1D forward DWT
	//
	MatWaveWavedec mw(wavename);

	size_t nlevels1d = mw.wmaxlev(NZ);

	size_t clen1d = mw.coefflength(NZ, nlevels1d);
	assert(clen1d == NZ);

	// Space for coefficients resulting from NX*NY 1D transforms
	//
    float* C = new float[ clen1d ];

	size_t L1d[(nlevels1d)+2];
	assert(L1d);
	mw.computeL(NZ, nlevels1d, L1d);

    mw.ResetFXformTimer();
    int rc = mw.wavedec( srcarr, NZ, nlevels1d, C, L1d );
    struct timeval timev = mw.GetFXformTimer();
    double seconds = timev.tv_sec + (timev.tv_usec/1000000.0);
    cerr << "1D forward transform took " << seconds << " seconds " << endl;

    for( int i = 0; i < 5; i++ ) {
        cerr << "\t" << C[i] << ",  " << C[ clen1d - i ] << endl;
    }

/*
    cerr << "printing original signal: " << endl;
    for( int i = 0; i < 9; i++ )
        cerr << srcarr[i] << endl;
    cerr << "printing 1D coeffs: " << endl;
    for( int i = 0; i < 9; i++ )
        cerr << C[i] << endl;
*/

    delete[] C;
}
	

int main(int argc, char **argv) {


	assert(argc == 3);
	string file = argv[1];
    int tmp = stoi( argv[2] );
    int NZ = tmp * 1024 * 1024;

	float *srcarr = new float[NX*NY*NZ];
	float *dstarr = NULL;



	FILE *fp = fopen(file.c_str(), "rb");
    int rc = fread(srcarr, sizeof(float), NX*NY*NZ, fp);
    fclose(fp);
    cerr << "read elements: " << rc << endl;

	string wname = "bior4.4";


    cout << "Test 1d\n";
    test1d( wname, srcarr, dstarr, NZ );


    delete[] srcarr;
    if( dstarr )        delete[] dstarr;
}

