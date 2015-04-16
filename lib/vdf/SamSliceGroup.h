
// A class that performs DWT on one dimension
// It   1) stores POINTERS to the coefficient arrays of each file
//      2) perform 1D DWT on all coefficients 
//      3) sorts all coefficients and picks out big ones
//      4) restores the coefficients by performing 1D IDWT


#ifndef _SamSliceGrouop_h_
#define _SamSliceGrouop_h_

#include <iostream>
#include "vapor/Compressor.h"

namespace VAPoR{

class SamSliceGroup{

public:
    SamSliceGroup();
    void Setup( int nslices, size_t rawlen, const vector< float* > &raw );
    ~SamSliceGroup();

    int SetupCompressor1 ( string wavename);
    
    int Compress1();

    // Decompress the coefficients from Compress1()
    int Decompress1( int ratio );


    // Get the ith coeff array. This is for the ith file.
    float* GetReconstructedPtr( int i );

    void FreeReconstructed( int i );
    void FreeCoeffs( );

    void Print1DRaw();
    void Print1DCoeffs();

protected:

    int _nslices;
    size_t _rawlen;                                     // length of each raw data array

    Compressor* _c1;

    vector< float* > _raw;                              // raw data ro process by this class

    vector< float* > _coeffs;                           // coefficients from DWT on raw arrays
    vector< SignificanceMap >*  _sigmapGroup;           // significance maps for coeffs

    vector< float* > _reconstructed;                    // reconstructed arrays from Decompress()

    

    float FindCoeffThreshold( int ratio );
};

}

#endif
