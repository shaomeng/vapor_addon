// This class (SamSliceGroup3) differs from SamSliceGroup2 in that
// this one puts all coefficients in a single array, and 
// performs all the transforms in place.
//
// A class that performs DWT on one dimension
// It   1) stores POINTERS to the coefficient arrays of each file
//      2) perform 1D DWT on all coefficients 
//      3) sorts all coefficients and picks out big ones
//      4) restores the coefficients by performing 1D IDWT


#ifndef _SamSliceGroup3_h_
#define _SamSliceGroup3_h_

#include <iostream>
#include <algorithm>
#include <cassert>
#include "vapor/MatWaveWavedec.h"

namespace VAPoR{

class SamSliceGroup3{

public:
    SamSliceGroup3( string wavename, size_t rawlen, size_t nslices );
    ~SamSliceGroup3();

    void UpdateRawPtr( vector< float* > &rawarr  );

    void Decompose();

    // Decompress the coefficients from Compress1()
    void Reconstruct( int ratio );

    // Get the ith coeff array. This is for the ith file.
    void FillReconstructedPtr( int idx, float* ptr );

/*
    void FreeReconstructed( int i );
    void FreeCoeffs( );
    void FreeRaw( int i );

    void Print1DRaw();
    void Print1DCoeffs();
*/

protected:
    string _wavename;
    size_t _rawlen;                    // length of each raw data array
    size_t _nslices;
    size_t _nlevels1d;
    size_t _clen1d;

    MatWaveWavedec* _mw;

    size_t* _L1d;           // bookkeeping array, shared among all 1D arrays
    float* _buf;            // huge array that has everything happening in place


    float FindCoeffThreshold( int ratio );
};

}

#endif
