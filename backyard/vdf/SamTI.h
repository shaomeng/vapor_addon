
// A class contains a number of time intervals
// It   1) stores POINTERS to the coefficient arrays of each file
//      2) picks out a number of time intervals to perform DWT
//      3) stores coefficients from DWT on the time intervals
//      4) restores the coefficients of time intervals by performing IDWT


#ifndef _SamTI_h_
#define _SamTI_h_

#include <iostream>
#include "vapor/Compressor.h"

namespace VAPoR{

class SamTI{

public:
    SamTI( int nfiles, size_t coeffArraySize, const vector< float* > &arrays );
    ~SamTI();

    int SetupTCompressor ( string wavename);
    
    // 1st way to pick up coefficients (and do DWT)
    // nTCoeffArray: performs DWT on first n coeffs 
    // TCoeffSize:   keep this number of coefficients for each DWT 
    int Compress1( size_t nTCoeffArray, size_t TCoeffSize  );

    // Decompress the coefficients from Compress1()
    int Decompress1();


    // Get the coefficients from decompression
    vector< float* >* GetDecompressedCoeffs()
    {
        return &_decompressedCoeffs;
    }    
    // Get the length of the decompressed arrays. 
    // This should be the same as the raw data size (totallen)
    size_t GetDecompressedArraySize() { return _decompressedCoeffSize; } 

    // Get the ith coeff array. This is for the ith file.
    float* GetDecompressedCoeff( size_t i );


protected:
    int _nfiles;

    vector< float* > _coeffArrays;  // coefficients in their positions 
                                    // _nfiles arrays in total
    size_t _coeffArraySize;         // size of each coefficient array


    vector< float* > _TCoeffArrays; // coefficients from DWT on Time Intervals
                                    // _nTCoeffArray in total
    size_t  _nTCoeffArray;          // number of positions to do DWT
                                    // this is <= _coeffArraySize
    int     _TCoeffSize;            // size of coefficients on each position
                                    // this is usually 2, 4, 8.


    vector< SignificanceMap* > _TSigmaps;   // _nTCoeffArray in total


    vector< float* > _decompressedCoeffs;   // _nfiles arrays, stores coefficients
    size_t _decompressedCoeffSize;          // should == _coeffArraySize

    Compressor* _c1;
    

};

}

#endif
