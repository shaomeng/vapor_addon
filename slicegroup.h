/* 
 * This class, slicegroup, is forked from previous SamSliceGroup3.
 *
 *
 * It   1) stores POINTERS to the coefficient arrays of each file
 *      2) perform 1D DWT on all coefficients 
 *      3) sorts all coefficients and picks out big ones
 *      4) restores the coefficients by performing 1D IDWT
 *
 * Programmer: Samuel Li
 * Date: 7/14/2015
 */

#ifndef _SLICEGROUP_
#define _SLICEGROUP_

#include <iostream>
#include <algorithm>
#include <cassert>
#include "vapor/MatWaveWavedec.h"
#include "cube3d.h"

namespace VAPoR{

class SliceGroup{

public:
    SliceGroup( string wavename );
    ~SliceGroup();

    void AddSlice( Cube3D* cube );
    void Initialize();

    void Decompose();

    void Reconstruct( int ratio );

    // Get the ith coeff array. This is for the ith file.
    void UpdateSlices( );


protected:
    string _wavename;
    size_t _sliceLen;                    // length of each slice of data 
    size_t _nslices;
    size_t _nlevels1d;

    MatWaveWavedec* _mw;

    vector< Cube3D* > _sliceVec;

    float* _buf;            // huge array that stores all the intermediate numbers 

    float FindCoeffThreshold( int ratio );
};

}

#endif
