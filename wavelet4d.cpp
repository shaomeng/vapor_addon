#include "wavelet4d.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

Wavelet4D::Wavelet4D( size_t NX, size_t NY, size_t NZ, size_t NT )
{

	_BLOCKDIM = 64;
	_wavename = "bior4.4";
	_cratio = 8;

	_NX = NX;
	_NY = NY;
	_NZ = NZ;
	_NT = NT;

	_filenames = NULL;
	_block_indices = NULL;
	CalcBlockIndices();
}

Wavelet4D::~Wavelet4D()
{
	if( _filenames ) {
		delete[] _filenames;
		_filenames = NULL;
	}
	if( _block_indices ) {
		delete[] _block_indices;
		_block_indices = NULL;
	}
}

void
Wavelet4D::SetFilePath( string path )
{
	_filepath= path;
}

void
Wavelet4D::SetFileStartIndex( size_t idx )
{
	_filenames = new string[ _NT ];
	char buf[256];

	for( size_t i = 0; i < _NT; i++ ) {
		sprintf( buf, "/vx.%04lu.out", i + idx );
		_filenames[i] = _filepath + buf;
	}	
}

void
Wavelet4D::CalcBlockIndices()
{
	assert( _NX % _BLOCKDIM == 0 );
	assert( _NY % _BLOCKDIM == 0 );
	assert( _NZ % _BLOCKDIM == 0 );
	
	size_t x_block = _NX / _BLOCKDIM;				// number of blocks in X
    size_t y_block = _NY / _BLOCKDIM;				// number of blocks in X
    size_t z_block = _NZ / _BLOCKDIM;				// number of blocks in X
    _BLOCKNUM = x_block * y_block * z_block;   

	_block_indices = new size_t[ _BLOCKNUM * 6 ];
	size_t counter = 0;

	/*
	 * 6 indices for a cell are stored together.
	 */
	for( size_t k = 0; k < z_block; k++ )
		for( size_t j = 0; j < y_block; j++ )
			for( size_t i = 0; i < x_block; i++ )
			{
                _block_indices[ counter++ ] = _BLOCKDIM *  i;		// startX
                _block_indices[ counter++ ] = _BLOCKDIM * (i+1);	// endX
                _block_indices[ counter++ ] = _BLOCKDIM *  j;		// startY
                _block_indices[ counter++ ] = _BLOCKDIM * (j+1);	// endY
                _block_indices[ counter++ ] = _BLOCKDIM *  k;		// startZ
                _block_indices[ counter++ ] = _BLOCKDIM * (k+1);	// endZ
			}
}

void
Wavelet4D::PrintBlockIndices()
{
	if( _block_indices == NULL )
		cerr << "_block_indices == NULL " << endl;
	else
		for( size_t i = 0; i < _BLOCKNUM; i++ )
			printf("%lu->%lu, %lu->%lu, %lu->%lu\n", 
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ], 
					_block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ], 
					_block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );
		
}

void
Wavelet4D::PrintFilenames()
{
	if( _filenames == NULL )
		cerr << "_filenames == NULL " << endl;
	else
		for( size_t i = 0; i < _NT; i++ )
			cerr << _filenames[i] << endl;
}



int
Wavelet4D::ParallelExec()
{
	/*
  	 * RMS and LMAX with same _BLOCKNUM but different _NT are stored together
	 */
	double* rms = new double[ _BLOCKNUM * _NT ];
	double* lmax = new double[ _BLOCKNUM * _NT ];

//printf("_BLOCKNUM = %ld, _NT = %ld\n", _BLOCKNUM, _NT );

	#pragma omp parallel for schedule( dynamic )
	for( size_t i = 0; i < _BLOCKNUM; i++ )
	{
//cerr << "i = " << i << endl;
		Cube3D** slices = new Cube3D*[ _NT ];

		for( size_t t = 0; t < _NT; t++ )
		{
//cerr << "t = " << t << endl;
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );
//cerr << "finish construction" << endl;
			slices[t] -> Decompose();
//cerr << "finish decompose" << endl;
			slices[t] -> Reconstruct (_cratio);
//cerr << "finish reconstruct" << endl;
			slices[t] -> Evaluate( rms[ i*_NT + t ], lmax[ i*_NT + t ] );
//cerr << "finish evaluate" << endl;
		}							

		for( size_t t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}

	for( size_t i = 0; i < _NT; i++ )
		cout << "RMS, LMAX: " << rms[i + 2*_NT] << "\t" << lmax[i + 2*_NT] << endl;

	delete[] rms;
	delete[] lmax;

	return 0;

/*
	size_t hw_concurrency = std::thread::hardware_concurrency();
	omp_set_dynamic(0);		// disables dynamic thread adjustment
	#pragma omp parallel num_threads( nthreads )
	{
		int my_tid = omp_get_thread_num();		// 0 to (nthreads - 1)
		printf("Hello World from thread = %d\n", tid);
		if (tid == 0) {
    		printf("Number of threads = %d\n", omp_get_num_threads());
    	} 
	}
*/
}


int main()
{
	Wavelet4D wav( 128, 128, 128, 20);
	string filepath = "/flash_buffer/Sam/HD_128";
	wav.SetFilePath( filepath );
	wav.SetFileStartIndex( 380 );
cerr << "finish set file start index" << endl;
	
	wav.ParallelExec();

}
