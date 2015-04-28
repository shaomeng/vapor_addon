// This class (SamSliceGroup2) differs from SamSliceGroup in that
// this one calls MatWaveWavedec directly.
//
// A class that performs DWT on one dimension
// It   1) stores POINTERS to the coefficient arrays of each file
//      2) perform 1D DWT on all coefficients 
//      3) sorts all coefficients and picks out big ones
//      4) restores the coefficients by performing 1D IDWT


#ifndef _SamSliceGroup2_h_
#define _SamSliceGroup2_h_

#include <iostream>
#include <algorithm>
#include <cassert>
#include "vapor/MatWaveWavedec.h"

namespace VAPoR{

class SamSliceGroup2{

public:
    SamSliceGroup2( string wavename, size_t rawlen, size_t nslices );
    ~SamSliceGroup2();

    // Expose the ith raw pointer.
    // The purpose was to let SamSlice2 hand it over.
    float* ExposeRawPtr( size_t i );

    void Decompose();

    // Decompress the coefficients from Compress1()
    void Reconstruct( int ratio );

    // Get the ith coeff array. This is for the ith file.
    float* GetReconstructedPtr( int i );

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

    vector< float* > _raw;              // raw data to process by this class
    vector< float* > _C1d;              // coefficients from DWT on raw arrays
    size_t* _L1d;             // bookkeeping array, shared among all 1D arrays

    vector< float* > _reconstructed;   

    float FindCoeffThreshold( int ratio );
};

}

#endif
