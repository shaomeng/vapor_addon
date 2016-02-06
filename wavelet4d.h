/*
 * This is the main library that Sam attemps to apply 4D wavelet widely.
 * It also employs parallel mechanisms to speed up the calculation.
 * Currently it divides big data size_to blocks of 64^3, 
 * and uses OpenMP to assign each block to a CPU core.
 *
 * Programmer: Samuel Li
 * Date: 12/8/2015
 */

#ifndef _Wavelet4D_
#define _Wavelet4D_

#include <iostream>
#include <cstdio>
#include <cassert>
#include <string>
#include <vector>
#include <omp.h>

#include "cube3d.h"
#include "slicegroup.h"


class Wavelet4D{

public:
	Wavelet4D( size_t NX, size_t NY, size_t NZ, size_t NT );
	~Wavelet4D();

	/*
	 * assume all file names start with the same prefix,
	 * and differ by an index.
	 */
	void GenerateFilenames( const std::string &path, int startIdx, const std::string &var );
	int ParallelExec();

	void PrintBlockIndices();
	void PrintFilenames();
	double FindMax( const double* arr, size_t len );
	double FindRMS( const double* arr, size_t len);


protected:
	size_t _NX, _NY, _NZ;	// spatial dimensions
	size_t _NT;				// temporal dimension
	size_t _BLOCKDIM;		// dimension of small blocks
	size_t _BLOCKNUM;		// total number of blocks in one time step
	int    _cratio;			// compression ratio
	
	std::vector<std::string> _filenames;
	std::string  _wavename;

	size_t* _block_indices;	// stores indices for each block at one time step
							// 6 indices to specify a block:
							// startX, endX, startY, endY, startZ, endZ

	void CalcBlockIndices();
};

#endif
