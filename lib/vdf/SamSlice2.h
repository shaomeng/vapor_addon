// SamSlice2 differs from SamSlice in that it doesn't call class Compressor
// but call class MatWaveWavedec directly.
// This effort can significantly reduce the memory footprint.

// A Time Slice
// A class that contains all the data from a file
// It   1) reads raw data from a file
//      2) stores coefficients from DWT
//      3) stores reconstructed grid

#ifndef _SamSlice2_h_
#define _SamSlice2_h_

#include <iostream>
#include <cstdio>
#include <algorithm>
#include "vapor/MatWaveWavedec.h"
#include <cassert>


namespace VAPoR{

class SamSlice2{
public:
    SamSlice2( string filename, string wavename, 
              size_t NX, size_t NY, size_t NZ );
    ~SamSlice2();
    int Decompose();
    int Reconstruct(int ratio);

    void UpdateCoeffs( const float* update );
/*
    void FreeRaw();
    void FreeCoeffs();
    void FreeReconstructed();


    float* GetRawPtr();
    float* GetCoeffsPtr();
    float* GetReconstructedPtr();

    
    void Print10Coeffs();
    void Print10Raws();
*/


protected:
    void ReadFile();
    string _filename;
    string _wavename;
    size_t _NX;
    size_t _NY;
    size_t _NZ;

    MatWaveWavedec* _mw;
    size_t _nlevels;
    size_t _clen;
    float* _C;
    size_t* _L;

    float* _raw;
    float* _reconstructed;

    float FindCoeffThreshold( int ratio );
};

}

#endif
