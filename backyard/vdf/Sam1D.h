
// A class that performs DWT on one dimension
// It   1) stores POINTERS to the coefficient arrays of each file
//      2) perform 1D DWT on all coefficients 
//      3) sorts all coefficients and picks out big ones
//      4) restores the coefficients by performing 1D IDWT


#ifndef _Sam1D_h_
#define _Sam1D_h_

#include <iostream>
//#include <algorithm>
#include "vapor/Compressor.h"

namespace VAPoR{

class Sam1D{

public:
    Sam1D( int nfiles, size_t rawlen, const vector< float* > &raw );
    ~Sam1D();

    int SetupCompressor1 ( string wavename);
    
    int Compress1();

    // Decompress the coefficients from Compress1()
    int Decompress1( int ratio );


    // Get the ith coeff array. This is for the ith file.
    float* GetDecompressedCoeff( size_t i );

    // print the thresholds in raw arrays that are 
    // cutting-off magnitude values
    void PrintRawThresholds( int ratio );


protected:
//    static bool mycompare (float i, float j) { return (i>j); }
    float FindCoeffThreshold( int ratio );

    int _nfiles;

    Compressor* _c1;

    vector< float* > _raw;     // raw data ro process by this class
    size_t _rawlen;            // length of each raw data array

    vector< float* > _coeffs;       // coefficients from DWT on raw arrays
    vector< SignificanceMap* > _sigmaps;   // significance maps for coeffs

    vector< float* > _reconstructed;   // reconstructed arrays from Decompress()

    

};

}

#endif
