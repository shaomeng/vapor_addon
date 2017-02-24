#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <sys/time.h>

#include "vapor/MatWaveWavedec.h"


const long int NX = 10;
const long int NY = 10;
const long int NZ = 10;
bool KeepAppCoeff = false;


template < typename T>
inline bool my_compare(const T &x1, const T &x2) 
{
    return(fabs(x1) > fabs(x2));
}

template< typename T>
void printPlaneX( const T* buf, long int x ) 		// prints a plane perpendicular to X axis.
{
	for( long int y = NY - 1; y >= 0; y-- )
	{
		printf("y = %ld:    ", y );
		for( long int z = 0; z < NZ; z++ )
		{
			printf("%f ", buf[z * NX * NY + y * NX + x] );	
		}	
		printf("\n");
	}
}


//
// Compress and decompress a 3D volume using a 3D wavelet transform
//
template< typename T>
void test3d(std::string wavename, const T *srcarr, T *dstarr, float cratio ) 
{
	for (size_t i=0; i<NX*NY*NZ; i++) dstarr[i] = 0.0;

	VAPoR::MatWaveWavedec mw(wavename);

	// Compute # of wavelet decompositions. Same # along each 
	// dimension
	//
	size_t nlevels = min(
		min(mw.wmaxlev(NX), mw.wmaxlev(NY)), mw.wmaxlev(NZ)
	);
  nlevels = 1;

	size_t clen = mw.coefflength3(NX, NY, NZ, nlevels);
	T *C = new T[clen];
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

	printPlaneX( C, 0 );

  size_t NumCoeff = (float) (NX*NY*NZ) / cratio;
	assert(NumCoeff >= startCoeffIdx);

	//
	// sort the coefficients  and find the threshold for culling
	// coefficients.
	//
	vector <T> sortedC; 
	for (size_t i=startCoeffIdx; i<clen; i++)  sortedC.push_back(C[i]);
	sort(sortedC.begin(), sortedC.end(), my_compare<T>);

	cout << "sortedC.size() = " << sortedC.size() << endl;

	size_t ti = NumCoeff - startCoeffIdx - 1;

	T threshold  = sortedC[ti];

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



template< typename T >
void compute_error(
    const T *data, const T *cdata, long int nx, long int ny, long int nz )
    //double &l1, double &l2, double &lmax, double &rms,
  	//double &min, double &max ) 
{
    double l1 = 0.0;
    double l2 = 0.0;
    double lmax = 0.0;
    double rms = 0.0;
    double sum = 0.0;
    double c = 0.0;
	  double min = data[0];
	  double max = data[0];
    for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                double delta = fabs (data[k*nx*ny + j*nx + i] - cdata[k*nx*ny + j*nx + i]);
                l1 += delta;
                l2 += delta * delta;
                lmax = delta > lmax ? delta : lmax;
                double y = delta * delta - c;
                double t = sum + y;
                c = (t - sum ) - y;
                sum = t;

                if( data[k*nx*ny + j*nx + i] > max )
                  max = data[ k*nx*ny + j*nx + i ]; 
                if( data[k*nx*ny + j*nx + i] < min )
                  min = data[ k*nx*ny + j*nx + i ]; 
            }
        }
    }
    l2 = sqrt(l2);
    rms = sqrt( sum / (1.0 * nx * ny * nz ));
    double range = max-min; 
    cout << "L1              = " << l1 << endl;
    cout << "L2              = " << l2 << endl;
    cout << "Data range      = " << range << endl;
    cout << "L-infy norm     = " << lmax << ", after normalization  = " << lmax/range << endl;
    cout << "RMSE            = " << rms  << ", after normalization  = " << rms/range << endl;
    cout << endl;
}


template< typename T >
void RectangleCopy( const T* src, size_t srcDimX, size_t srcDimY,
                          T* dst, size_t dstDimX, size_t dstDimY,
                          size_t srcStartX, size_t srcStartY,
                          size_t dstStartX, size_t dstStartY,
                          size_t copyDimX,  size_t copyDimY )
{
  for(   size_t y = 0; y < copyDimY; y++ )
    for( size_t x = 0; x < copyDimX; x++ )
    {
      size_t srcIdx = (y + srcStartY) * srcDimX + x + srcStartX;
      size_t dstIdx = (y + dstStartY) * dstDimX + x + dstStartX;
      dst[dstIdx]   = src[srcIdx];
    }
}




int main(int argc, char* argv[] ) {

	assert(argc == 3);
  float cratio = atof( argv[1] );
	std::string file = argv[2];

	typedef double T;
	T *srcarr = new T[NX * NY * NZ];
	T *dstarr = new T[NX * NY * NZ];

	FILE *fp = fopen(file.c_str(), "r");
  size_t rc = fread(srcarr, sizeof(T), NX * NY * NZ, fp);
  fclose(fp);
  assert (rc == NX * NY * NZ);

	std::string wname = "bior4.4";
  
  //test2d( wname, srcarr, dstarr, cratio, dimX, dimY );
  struct timeval start, end;
  gettimeofday( &start, NULL );
  test3d( wname, srcarr, dstarr, cratio );
  gettimeofday( &end, NULL );
  double t = (double)( (end.tv_sec * 1000000 + end.tv_usec) -
             (start.tv_sec * 1000000 + start.tv_usec) )/1000000.0;
  cout << "total time      = " << t << endl;

	compute_error( srcarr, dstarr, NX, NY, NZ ); 


  delete[] srcarr;
  delete[] dstarr;
}

