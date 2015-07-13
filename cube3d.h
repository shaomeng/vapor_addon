/*
 * This class, Cube3D, is forked from the previous SamSlice2.
 * This class aims to be simpler: it doesn't store the raw data,
 * but only the latest coefficients from DWT.
 * As a result, it mainly handles file I/O, 
 * and calls MatWaveWavedec to perform 3D DWT/IDWT.
 *
 * Programmer: Samuel Li
 * Date: 7/12/2015
 */


#ifndef _Cube3D_
#define _Cube3D_

#include <iostream>
#include <cstdio>
#include <algorithm>
#include "vapor/MatWaveWavedec.h"
#include <cassert>


namespace VAPoR{

class Cube3D{
public:
    Cube3D( string filename, string wavename, 
              size_t NX, size_t NY, size_t NZ );
    ~Cube3D();
    int Decompose();
    int Reconstruct(int ratio);

    void Print10Elements();


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

/*
    float* _raw;
    float* _reconstructed;
*/

    float FindCoeffThreshold( int ratio );
};

}

#endif
