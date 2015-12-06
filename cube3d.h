/*
 * This class, Cube3D, is forked from the previous SamSlice2.
 * This class aims to be simpler: it doesn't store the raw data,
 * but only the latest coefficients from DWT.
 * As a result, it mainly handles file I/O, 
 * and calls MatWaveWavedec to perform 3D DWT/IDWT.
 *
 * Programmer: Samuel Li
 * Date: 7/12/2015
 *
 * Revision: 11/27/2015
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
    Cube3D( string filename, string wavename, 
              size_t NX, size_t NY, size_t NZ,
			  size_t NX_total, size_t NY_total, size_t NZ_total,
			  size_t startX, size_t endX,
			  size_t startY, size_t endY,
			  size_t startZ, size_t endZ );
    ~Cube3D();
    int Decompose();
    int Reconstruct(int ratio);

    float GetCoeff( size_t idx );
    void PutCoeff( size_t idx, float c );
    void CullCoeffs( float t ); // the threshold, t, must be positive.
    void Evaluate( double& rms, double& lmax );

    void Print10Elements();

    void ReloadInputFile()  { ReadFile( _C ); }
    size_t GetCoeffLen()    { return _clen; }

protected:
    void ReadFile( float* buf );
    void ReadFileChunck( float* buf );
    string _filename;
    string _wavename;

    size_t _NX;			// dimension of this cube working on
    size_t _NY;
    size_t _NZ;

	size_t _NX_total;	// dimension of the total data brick
	size_t _NY_total;
	size_t _NZ_total;
	size_t _startX;		// Index range of current cube
	size_t _endX;
	size_t _startY;
	size_t _endY;
	size_t _startZ;
	size_t _endZ;

    MatWaveWavedec* _mw;
    size_t _nlevels;
    size_t _clen;
    float* _C;
    size_t* _L;

    float FindCoeffThreshold( int ratio );
};

}

#endif
