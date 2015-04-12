// A Time Slice
// A class that contains all the data from a file
// It   1) reads raw data from a file
//      2) stores coefficients from DWT
//      3) stores reconstructed grid

#ifndef _SamTS_h_
#define _SamTS_h_

#include <iostream>
// #include <cstdio>
#include "vapor/Compressor.h"
//#include "vapor/SignificanceMap.h"


namespace VAPoR{

using std::vector;
using std::string;

class SamTS{
public:
    SamTS( string filename );
    ~SamTS();
    void Print1stRaw();
    void SetupSigmap( size_t cvectorsize, int nsigmaps);
    void SetupReconstructedArray();
    void PutCoeffsInPosition( );
//    void SetLOD( int lod )  { _reconstructedLOD = lod; }

    float* GetRawPtr()      { return _raw; };
    float* GetCoeffsPtr()   { return _coeffs; }
    float* GetReconstructedPtr()                { return _reconstructed; }
    float* GetCoeffsInPosition()                { return _coeffsInPosition; }

    size_t GetRawSize()                         { return _rawsize; }
    size_t GetCoeffsSize()                      { return _cvectorsize; }
    size_t GetReconstructedSize()               { return _reconstructedSize;}

    vector< SignificanceMap >*  GetOriginSigMapsPtr()    { return &_originSigmaps; }
    
    float* GetReconstructed2Ptr();

    // Recover the original array based on passed in compressor and coefficients.
//    int Decompress( Compressor* c, float* recoveredCoeffs ); 


// Data fields of a SamTS object:
protected:
    string _filename;

    float* _raw;
    size_t _rawsize;

    float* _coeffs;
    size_t _cvectorsize;        // should be the same as _rawsize

    // Sigmaps for the orignial coefficient array
    vector< SignificanceMap > _originSigmaps;

    float* _reconstructed;
    size_t _reconstructedSize;  // should be the same as _rawsize.
//    int    _reconstructedLOD;   // the lod used to reconstruct.


    float* _coeffsInPosition;   // the same values as in _coeffs, 
                                // in non-prioritized positions/ordering.

    float* _reconstructed2;       // generated from the decompressed coefficients
                                // should be the same lengh as _rawsize;
//    SignificanceMap* _decompressSigmap;     // manually created sigmap
};

}

#endif
