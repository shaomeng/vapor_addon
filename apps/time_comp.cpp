#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

#include "vapor/MatWaveWavedec.h"
using namespace VAPoR;


const int NX = 504;
const int NY = 504;
const int NZ = 504;
bool KeepAppCoeff = false;
float cratio = 4.0;

const size_t NumCoeff = (float) (NX*NY*NZ) / cratio;

inline bool my_compare(const float &x1, const float &x2) {
    return(fabsf(x1) > fabsf(x2));
}

//
// Compress and decompress a 3D volume using a 3D wavelet transform
//
void test3d(string wavename, const float *srcarr, float *dstarr) {

	for (size_t i=0; i<NX*NY*NZ; i++) dstarr[i] = 0.0;

	MatWaveWavedec mw(wavename);

	// Compute # of wavelet decompositions. Same # along each 
	// dimension
	//
	size_t nlevels = min(
		min(mw.wmaxlev(NX), mw.wmaxlev(NY)), mw.wmaxlev(NZ)
	);

	size_t clen = mw.coefflength3(NX, NY, NZ, nlevels);
	float *C = new float[clen];
	assert(C);

	size_t L[ (21*nlevels) + 6];

	// Compute the bookkeeping vector, L.
	//
	// N.B. L will be recomputed by wavedec(), but not
	// waverec();
	//
	mw.computeL3(NX, NY, NZ, nlevels, L);

	size_t startCoeffIdx = 0;
	if (KeepAppCoeff) startCoeffIdx = L[0]*L[1]*L[2];
	int rc = mw.wavedec3(srcarr, NX, NY, NZ, nlevels, C, L);
	assert (rc>=0);

	assert(NumCoeff >= startCoeffIdx);

	//
    // sort the coefficients  and find the threshold for culling
	// coefficients.
    //
   	vector <float> sortedC; 
    for (size_t i=startCoeffIdx; i<clen; i++)  sortedC.push_back(C[i]);
    sort(sortedC.begin(), sortedC.end(), my_compare);

	cout << "sortedC.size() = " << sortedC.size() << endl;

	size_t ti = NumCoeff - startCoeffIdx - 1;

	double threshold  = sortedC[ti];

	// Zero out wavelet coefficients that are smaller than threshold
	//
	cout << "ti " << ti << endl;
	cout << "threshold  " << threshold  << endl;
    size_t should_squash = clen - clen / cratio;    // should squash this many
	size_t squashed = 0;
	for (size_t i=startCoeffIdx; i<clen; i++) {
		if (fabs(C[i]) <= fabs(threshold) && squashed <= should_squash) {
			squashed += 1;
			C[i] = 0;
		}
	}
	cout << "squashed " << squashed << endl;
	
	// Inverse DWT
	//
	mw.waverec3(C,L,nlevels,dstarr);	

    delete[] C;
}

// 
// Compress and decompress a 3D volume using a 2D wavelet transform
// transform along X and Y axes, followed by a 1D transform along
// the Z axis.
//
void test2dp1d(string wavename, const float *srcarr, float *dstarr) 
{

	for (size_t i=0; i<NX*NY*NZ; i++) dstarr[i] = 0.0;

	// 
	// 2D transform first
	//
	MatWaveWavedec mw(wavename);

	size_t nlevels2d = min(mw.wmaxlev(NX), mw.wmaxlev(NY));

	size_t clen2d = mw.coefflength2(NX, NY, nlevels2d);
	assert (clen2d = NX*NY);

	// space for 2D wavelet coefficients
	//
	vector <float *> C2d;
	for (int i=0; i<NZ; i++) C2d.push_back(new float[clen2d]);

	size_t L2d[(6*nlevels2d)+4];

	// Compute the bookkeeping vector, L2d.
	//
	mw.computeL2(NX, NY, nlevels2d, L2d);

	// perform 2D transforms, one for each of the NZ slices
	//
	for (int i=0; i<NZ; i++) {
		int rc = mw.wavedec2(srcarr+i*NX*NY, NX, NY, nlevels2d, C2d[i], L2d);
		assert (rc>=0);
	}

	//
	// Now do 1D forward DWT
	//

	size_t nlevels1d = mw.wmaxlev(NZ);

	size_t clen1d = mw.coefflength(NZ, nlevels1d);
	assert(clen1d == NZ);

	// Space for coefficients resulting from NX*NY 1D transforms
	//
	vector <float *> buf1d, C1d;
	for (int i=0; i<NX*NY; i++) {
		C1d.push_back(new float[clen1d]);
		buf1d.push_back(new float[clen1d]);
	}

	size_t L1d[(nlevels1d)+2];
	assert(L1d);

	mw.computeL(NZ, nlevels1d, L1d);

	// Copy 2D coefficients to buf1d, tranpose array in process
	//
	for (int i=0; i<NX*NY; i++) {
	for (int j=0; j<NZ; j++) {
		buf1d[i][j] = C2d[j][i];
	}
	}

// Sam
//cerr << "Sam: printing the first 18 elements before DWT: " << endl;
//for( int i = 0; i < NZ; i++ )
//    cerr << "\t" << buf1d[0][i] << endl;

	// Forward 1D DWTs
	//
	for (int i=0; i<NX*NY; i++) {
		int rc = mw.wavedec(buf1d[i], NZ, nlevels1d, C1d[i], L1d);
		assert (rc>=0);
	}

// Sam
//cerr << "Sam: printing the first 18 elements after DWT: " << endl;
//for( int i = 0; i < NZ; i++ )
//    cerr << "\t" << C1d[0][i] << endl;

	// Copy coefficients to single vector for easy sorting
	//
	vector <float> sortedC;
	for (int i=0; i<NX*NY; i++) {
		for (size_t j=0; j<clen1d; j++)  {
			sortedC.push_back(C1d[i][j]);
		}
	}

	cout << "sortedC.size() = " << sortedC.size() << endl;

	//
    // sort the coefficients 
    //
    sort(sortedC.begin(), sortedC.end(), my_compare);

	size_t ti = NumCoeff - 1;

	double threshold  = sortedC[ti];

	// Kill coefficients that are smaller than threshold
	//
	cout << "ti " << ti << endl;
	cout << "threshold  " << threshold  << endl;
	size_t squashed = 0;
	for (int i=0; i<NX*NY; i++) {
	for (size_t j=0; j<clen1d; j++) {
		if (fabs(C1d[i][j]) <= fabs(threshold )) {
			squashed += 1;
			C1d[i][j] = 0;
		}
	}
	}
	cout << "squashed " << squashed << endl;

	//
	// Now perform IDWTs. First 1D, than 2D transforms
	//

	for (int i=0; i<NX*NY; i++) {
		int rc = mw.waverec(C1d[i], L1d, nlevels1d, buf1d[i]);
		assert (rc>=0);
	}

	// copy and tranpose coefficients
	//
	for (int i=0; i<NX*NY; i++) {
	for (int j=0; j<NZ; j++) {
		C2d[j][i] = buf1d[i][j];
	}
	}

	for (int i=0; i<NZ; i++) {
		int rc = mw.waverec2(C2d[i], L2d, nlevels2d, dstarr+i*NX*NY);
		assert (rc>=0);
	}
	
    for( int i = 0; i < C2d.size(); i++ )
        if( C2d[i] )    delete[] C2d[i];
    for( int i = 0; i < C1d.size(); i++ )
        if( C1d[i] )    delete[] C1d[i];
    for( int i = 0; i < buf1d.size(); i++ )
        if( buf1d[i] )    delete[] buf1d[i];
}
	

void compute_error(
    const float *data, const float *cdata, int nx, int ny, int nz,
    double &l1, double &l2, double &lmax, double &rms
) {
    l1 = 0.0;
    l2 = 0.0;
    lmax = 0.0;
    rms = 0.0;
    float sum = 0.0;
    float c = 0.0;
    for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                double delta = fabs (data[k*nx*ny + j*nx + i] - cdata[k*nx*ny + j*nx + i]);
                l1 += delta;
                l2 += delta * delta;
                lmax = delta > lmax ? delta : lmax;
                float y = delta * delta - c;
                float t = sum + y;
                c = (t - sum ) - y;
                sum = t;
            }
        }
    }
//cerr << "sum, l2 are " << sum << "  " << l2 << endl;
    l2 = sqrt(l2);
//    rms = l2 / (nx*ny*nz);
    rms = sqrt( sum / (1.0 * nx * ny * nz ));
}

int main(int argc, char **argv) {

	assert(argc == 2);
	string file = argv[1];
//    string file = "/Users/samuel/Backyard/256cubes/e0.float";

	float *srcarr = new float[NX*NY*NZ];
	float *dstarr = new float[NX*NY*NZ];



	FILE *fp = fopen(file.c_str(), "r");
    int rc = fread(srcarr, sizeof(*srcarr), NX*NY*NZ, fp);
    fclose(fp);
    assert (rc == NX*NY*NZ);

	string wname = "bior4.4";


	cout << "Test 3d\n";
	test3d(wname, srcarr, dstarr);

    double l1, l2, lmax, rms;
	compute_error(
		srcarr, dstarr, NX, NY, NZ, l1, l2, lmax, rms
	);

    cout << "L1 = " << l1 << endl;
    cout << "L2 = " << l2 << endl;
    cout << "LMax = " << lmax << endl;
    cout << "RMS = " << rms << endl;
	cout << endl;
	cout << endl;
	
/*
	cout << "Test 2dp1d \n";
	test2dp1d(wname, srcarr, dstarr);

	compute_error(
		srcarr, dstarr, NX, NY, NZ, l1, l2, lmax, rms
	);

    cout << "L1 = " << l1 << endl;
    cout << "L2 = " << l2 << endl;
    cout << "LMax = " << lmax << endl;
    cout << "RMS = " << rms << endl;
*/

    delete[] srcarr;
    delete[] dstarr;
}

