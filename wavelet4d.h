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

#define EVALUATE
#define INOTIFY

#ifdef INOTIFY
#include <sys/inotify.h>
#include <unistd.h>
#define BUF_LEN (100 * (sizeof(struct inotify_event) + 255 + 1))
#endif

class Wavelet4D{

public:
	Wavelet4D( long NX, long NY, long NZ, long NT );
	~Wavelet4D();

	void SetPath( const std::string &path );
	/*
	 * Filenames differ at the index number at the last.
	 */
	void GenerateFilenames( const std::string &name, long startIdx );
	int ParallelExec();

	void PrintBlockIndices();
	void PrintFilenames();
	double FindMax( const double* arr, long len );
	double FindRMS( const double* arr, long len);

	void SetCRatio( int i )	{ _cratio = i; }

#ifdef INOTIFY
	void StartMonitor();
#endif

protected:
	long _NX, _NY, _NZ;	// spatial dimensions
	long _NT;				// temporal dimension
	long _BLOCKDIM;		// dimension of small blocks
	long _BLOCKNUM;		// total number of blocks
	int    _cratio;			// compression ratio
	
	std::vector<std::string> _filenames;
	std::string  _wavename;

    std::string _path;
	

	long* _block_indices;	// stores indices for each block at one time step
							// 6 indices to specify a block:
							// startX, endX, startY, endY, startZ, endZ

	void CalcBlockIndices();

#ifdef INOTIFY
	/*
	 * Test if the last file in _filenames has finished writting.
	 */
	bool FinishWriteLast(struct inotify_event *i);
#endif

};

#endif
