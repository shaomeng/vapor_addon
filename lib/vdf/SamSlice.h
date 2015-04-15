// A Time Slice
// A class that contains all the data from a file
// It   1) reads raw data from a file
//      2) stores coefficients from DWT
//      3) stores reconstructed grid

#ifndef _SamSlice_h_
#define _SamSlice_h_

#include <iostream>
#include <cstdio>
#include "vapor/Compressor.h"
#include <cassert>


namespace VAPoR{

using std::vector;
using std::string;

class SamSlice{
public:
    SamSlice( string filename );
    SamSlice( );
    ~SamSlice();
    void ReadFile( string filename );
    int SetupCompressor1( vector< size_t > &dims, string wavename );
    int Compress1();
    int Decompress1(int ratio);
    void PutCoeffsInPosition();
    void ReplaceCoeffs( float* replacement );


    float* GetRawPtr()      { return _raw; };
    float* GetCoeffsPtr()   { return _coeffs; }
    float* GetReconstructedPtr()                { return _reconstructed; }

    size_t GetRawSize()                         { return _rawlen; }
    
    void Print10Coeffs();
    void Print10Raws();

protected:
    string _filename;
    vector< size_t > _dims;
    Compressor* _c1;

    float* _raw;
    size_t _rawlen;

    float* _coeffs;

    vector< SignificanceMap > _sigmaps; // This vector is supposed to have only 1 element.

    float* _reconstructed;

    float FindCoeffThreshold( int ratio );
};

}

#endif
