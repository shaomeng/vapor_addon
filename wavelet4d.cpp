#include "wavelet4d.h"


using std::string;
using std::cout;
using std::cerr;
using std::endl;


static float
FindMin( const float* arr, long len ) {
    float min = arr[0];
    for( long i = 1; i < len; i++ )
        if( arr[i] < min )
            min = arr[i];
    return min;
}

static float
FindMax( const float* arr, long len ) {
    float max = arr[0];
    for( long i = 1; i < len; i++ )
        if( arr[i] > max )
            max = arr[i];
    return max;
}

static double 
FindMax( const double* arr, long len ) {
    double max = arr[0];
    for( long i = 1; i < len; i++ )
        if( arr[i] > max )
            max = arr[i];
    return max;
}

static double 
FindRMS( const double* arr, long len)
{
    double sum = 0.0;
    double c = 0.0;
    for( long i = 0; i < len; i++ ) {
        double y = arr[i] * arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= double(len);
    sum = sqrt( sum );
    return sum;
}
static double 
FindRMS2( const double* arr, long len, long denominator)
{
    double sum = 0.0;
    double c = 0.0;
    for( long i = 0; i < len; i++ ) {
        double y = arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= double(denominator);
    sum = sqrt( sum );
    return sum;
}


/* constructor */
Wavelet4D::Wavelet4D( long NX, long NY, long NZ, long NT )
{

	_BLOCKDIM = 70;
	_wavename = "bior4.4";
	_cratio = 1;

	_NX = NX;
	_NY = NY;
	_NZ = NZ;
	_NT = NT;

	_block_indices = NULL;
	_path.clear();
	CalcBlockIndices();
	_reconstruct_buf = NULL;
}

Wavelet4D::~Wavelet4D()
{
	if( _block_indices ) {
		delete[] _block_indices;
		_block_indices = NULL;
	}
	if( _reconstruct_buf ) {
		for( long i = 0; i < _NT; i++ )
			if( _reconstruct_buf[i] ) {
				delete[] _reconstruct_buf[i];
				_reconstruct_buf[i] = NULL;
			}
		delete[] _reconstruct_buf;
		_reconstruct_buf = NULL;
	}
}

void 
Wavelet4D::SetPath( const std::string &path )
{
	_path.clear();
	_path = path;
}

void
Wavelet4D::GenerateFilenames( const string &name, long idx)
{
	assert( !_path.empty());
	char buf[256];
	_filenames.clear();

	int step = 2;

	for( long i = 0; i < _NT; i++ ) {
		sprintf( buf, "%04ld.float", i*step + idx);
		string f = _path + "/" + name + buf;
		_filenames.push_back( f );
	}	
}

void
Wavelet4D::GenerateBkpFilenames( const string &path, const string &name, long idx)
{
	char buf[256];
	_bkpFilenames.clear();

	for( long i = 0; i < _NT; i++ ) {
		sprintf( buf, "%02ld.float", i+idx);
		string f = path + "/" + name + buf;
		_bkpFilenames.push_back( f );
	}	
}

void
Wavelet4D::CalcBlockIndices()
{
	assert( _NX % _BLOCKDIM == 0 );
	assert( _NY % _BLOCKDIM == 0 );
	assert( _NZ % _BLOCKDIM == 0 );
	
	long x_block = _NX / _BLOCKDIM;				// number of blocks in X
    long y_block = _NY / _BLOCKDIM;				// number of blocks in X
    long z_block = _NZ / _BLOCKDIM;				// number of blocks in X
    _BLOCKNUM = x_block * y_block * z_block;   

	_block_indices = new long[ _BLOCKNUM * 6 ];
	long counter = 0;

	/*
	 * 6 indices for a cell are stored together.
	 */
	for( long k = 0; k < z_block; k++ )
		for( long j = 0; j < y_block; j++ )
			for( long i = 0; i < x_block; i++ )
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
		for( long i = 0; i < _BLOCKNUM; i++ )
			printf("%lu->%lu, %lu->%lu, %lu->%lu\n", 
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ], 
					_block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ], 
					_block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );
		
}

void
Wavelet4D::PrintFilenames()
{
	for( unsigned int i = 0; i < _filenames.size(); i++ )
		cout << _filenames[i] << endl;
}



int
Wavelet4D::ParallelExec()
{
#ifdef EVALUATE
	/*
  	 * RMS and LMAX with same block index but different _NT are stored together
	 */
	double* rms3d = new double[ _BLOCKNUM * _NT ];
	double* lmax3d = new double[ _BLOCKNUM * _NT ];
	double* rms4d = new double[ _BLOCKNUM * _NT ];
	double* lmax4d = new double[ _BLOCKNUM * _NT ];
	float*  min4d  = new float[ _BLOCKNUM * _NT ];
	float*  max4d  = new float[ _BLOCKNUM * _NT ];
	double*  squaresum3d  = new double[ _BLOCKNUM * _NT ];
	double*  squaresum4d  = new double[ _BLOCKNUM * _NT ];
	long*   numNonSpecial = new long[ _BLOCKNUM * _NT ];
#endif

	/*
	 * Each thread takes care of one index from _BLOCKNUM.
	 */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );
	//printf("system thread = %d, _BLOCKNUM = %d\n", nthreadSys, _BLOCKNUM );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];
		SliceGroup* group = new SliceGroup( _wavename );

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

#ifdef EVALUATE
			slices[t] -> GetMinMax( min4d[ i*_NT + t ], max4d[ i*_NT + t ] );

			/* individual slice reconstruction and evaluation */
			slices[t] -> Decompose();
			slices[t] -> Reconstruct (_cratio);
			slices[t] -> Evaluate( rms3d[ i*_NT + t ], lmax3d[ i*_NT + t ] );
			//slices[t] -> EvaluateWithAnotherFile( _bkpFilenames[t], numNonSpecial[ i*_NT + t],
			//									 squaresum3d[ i*_NT + t ], lmax3d[ i*_NT + t ] );
			slices[t] -> ReloadInputFile();
#endif
			slices[t] -> Decompose();
			group -> AddSlice( slices[t] );
		}							
		/* Temporal compression */
		group -> Initialize();
		group -> Decompose();		// All coefficients stored now!

		/* Output coefficients 
		string coeff_name = _filenames[0] + ".block" + to_string(i) + ".coeff";
		group -> OutputFile( coeff_name, _cratio );
		*/

#ifdef EVALUATE
		group -> Reconstruct( _cratio );
		group -> UpdateSlices();
		for( long t = 0; t < _NT; t++ ) {
			slices[t] -> Reconstruct(1);
			slices[t] -> Evaluate( rms4d[ i*_NT + t ], lmax4d[ i*_NT + t ] );
			//slices[t] -> EvaluateWithAnotherFile( _bkpFilenames[t], numNonSpecial[ i*_NT + t],
			//									 squaresum4d[ i*_NT + t ], lmax4d[ i*_NT + t ] );
		}
#endif

		if( group )				delete group;
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */

#ifdef EVALUATE
	float min = FindMin( min4d, _BLOCKNUM * _NT );
	float max = FindMax( max4d, _BLOCKNUM * _NT );
	printf("In %ld steps, min=%f, max=%f\n", _NT, min, max ); 
	/* Compile results from Evaluate() */
	/*
	printf("3D DWT RMS = %e, LMAX = %e\n", FindRMS( rms3d,  _BLOCKNUM * _NT),
										   FindMax( lmax3d, _BLOCKNUM * _NT));
	printf("4D DWT RMS = %e, LMAX = %e\n", FindRMS( rms4d,  _BLOCKNUM * _NT),
										   FindMax( lmax4d, _BLOCKNUM * _NT));
	*/

	/* Compile results from EvaluateWithAnotherFile() */
 	/*
	long total_non_special = 0;
	for( long i = 0; i < _BLOCKNUM * _NT; i++ )
		total_non_special += numNonSpecial[i];
	printf("3D DWT RMS = %e, LMAX = %e\n", FindRMS2( squaresum3d, _BLOCKNUM * _NT, total_non_special),
										   FindMax( lmax3d, _BLOCKNUM * _NT));
	printf("4D DWT RMS = %e, LMAX = %e\n", FindRMS2( squaresum4d, _BLOCKNUM * _NT, total_non_special),
										   FindMax( lmax4d, _BLOCKNUM * _NT));
	*/

	double range = max - min;
	printf("3D DWT NRMS = %e, NLMAX = %e\n", FindRMS( rms3d, _BLOCKNUM * _NT) / range,
										    FindMax( lmax3d, _BLOCKNUM * _NT) / range);
	printf("4D DWT NRMS = %e, NLMAX = %e\n", FindRMS( rms4d, _BLOCKNUM * _NT) / range,
										    FindMax( lmax4d, _BLOCKNUM * _NT) / range);

	delete[]  squaresum3d;
	delete[]  squaresum4d;
	delete[]  numNonSpecial;
	delete[] rms3d;
	delete[] lmax3d;
	delete[] rms4d;
	delete[] lmax4d;
	delete[] min4d;
	delete[] max4d;
#endif

	return 0;
}

void
Wavelet4D::Output4DReconstruct()
{
	if( _reconstruct_buf == NULL ) 
	{
		_reconstruct_buf = new float*[ _NT ];
		for( long i = 0; i < _NT; i++ )
			_reconstruct_buf[i] = new float[ _NX * _NY * _NZ ];
	}

	/* * Each thread takes care of one index from _BLOCKNUM.  */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];
		SliceGroup* group = new SliceGroup( _wavename );

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

			slices[t] -> Decompose();
			group -> AddSlice( slices[t] );
		}							
		/* Temporal compression */
		group -> Initialize();
		group -> Decompose();		// All coefficients stored now!
		group -> Reconstruct( _cratio );
		group -> UpdateSlices();
		for( long t = 0; t < _NT; t++ ) {
			slices[t] -> Reconstruct(1);

			for( long z = 0; z < _BLOCKDIM; z++ )
				for( long y = 0; y < _BLOCKDIM; y++ )
					for( long x = 0; x < _BLOCKDIM; x++ )
					{
						float v = slices[t] -> GetCoeff( x, y, z );
						PutValToWrite( x + _block_indices[6*i],
									   y + _block_indices[6*i+2],
									   z + _block_indices[6*i+4],
									   t, v );
					}
		}

		if( group )				delete group;
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */


	/* outputs files from 4D reconstruction onto disk */
	for( long t = 0; t < _NT; t++ )
	{
		string ratio_str = to_string( _cratio );
		string name3d = _filenames[t] + "." + ratio_str + "c4d";
		FILE* f = fopen( name3d.c_str(), "wb" );
		fwrite( _reconstruct_buf[t], sizeof(float), _NX*_NY*_NZ, f );
		fclose(f);
	}
}



void
Wavelet4D::Output3DReconstruct()
{
	if( _reconstruct_buf == NULL ) 
	{
		_reconstruct_buf = new float*[ _NT ];
		for( long i = 0; i < _NT; i++ )
			_reconstruct_buf[i] = new float[ _NX * _NY * _NZ ];
	}

	/* * Each thread takes care of one index from _BLOCKNUM.  */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];
		SliceGroup* group = new SliceGroup( _wavename );

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

			slices[t] -> Decompose();
			slices[t] -> Reconstruct (_cratio);
		
			for( long z = 0; z < _BLOCKDIM; z++ )
				for( long y = 0; y < _BLOCKDIM; y++ )
					for( long x = 0; x < _BLOCKDIM; x++ )
					{
						float v = slices[t] -> GetCoeff( x, y, z );
						PutValToWrite( x + _block_indices[6*i],
									   y + _block_indices[6*i+2],
									   z + _block_indices[6*i+4],
									   t, v );
					}
		}							

		if( group )				delete group;
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */

	/* outputs files from 3D reconstruction onto disk */
	for( long t = 0; t < _NT; t++ )
	{
		string ratio_str = to_string( _cratio );
		string name3d = _filenames[t] + "." + ratio_str + "c3d";
		FILE* f = fopen( name3d.c_str(), "wb" );
		fwrite( _reconstruct_buf[t], sizeof(float), _NX*_NY*_NZ, f );
		fclose(f);
	}
}

int
Wavelet4D::TimeDimParallelExec()
{
	/*
  	 * RMS and LMAX with same block index but different _NT are stored together
	 */
	double* rms4d = new double[ _BLOCKNUM * _NT ];
	double* lmax4d = new double[ _BLOCKNUM * _NT ];
	float*  min4d  = new float[ _BLOCKNUM * _NT ];
	float*  max4d  = new float[ _BLOCKNUM * _NT ];

	/*
	 * Each thread takes care of one index from _BLOCKNUM.
	 */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];
		SliceGroup* group = new SliceGroup( _wavename );

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

			slices[t] -> GetMinMax( min4d[ i*_NT + t ], max4d[ i*_NT + t ] );

			group -> AddSlice( slices[t] );
		}							
		/* Temporal compression */
		group -> Initialize();
		group -> Decompose();		// All coefficients stored now!


		group -> Reconstruct( _cratio );

		group -> UpdateSlices();


		for( long t = 0; t < _NT; t++ ) {
			slices[t] -> Evaluate( rms4d[ i*_NT + t ], lmax4d[ i*_NT + t ] );
		}

		if( group )				delete group;
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */

	float min = FindMin( min4d, _BLOCKNUM * _NT );
	float max = FindMax( max4d, _BLOCKNUM * _NT );
	printf("In %ld steps, min=%f, max=%f\n", _NT, min, max ); 

	double range = max - min;
	printf("1D DWT Time NRMS = %e, NLMAX = %e\n", FindRMS( rms4d, _BLOCKNUM * _NT) / range,
										    FindMax( lmax4d, _BLOCKNUM * _NT) / range);

	delete[] rms4d;
	delete[] lmax4d;
	delete[] min4d;
	delete[] max4d;

	return 0;
}

void
Wavelet4D::XDimParallelExec()
{
	/*
  	 * RMS and LMAX with same block index but different _NT are stored together
	 */
	double* rmse1d = new double[ _BLOCKNUM * _NT ];
	double* lmax1d = new double[ _BLOCKNUM * _NT ];
	float*  min1d  = new float[ _BLOCKNUM * _NT ];
	float*  max1d  = new float[ _BLOCKNUM * _NT ];

	/*
	 * Each thread takes care of one index from _BLOCKNUM.
	 */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

			slices[t] -> GetMinMax( min1d[ i*_NT + t ], max1d[ i*_NT + t ] );

			/* individual slice reconstruction and evaluation */
			slices[t] -> DecomposeX();
			slices[t] -> ReconstructX (_cratio);
			slices[t] -> Evaluate( rmse1d[ i*_NT + t ], lmax1d[ i*_NT + t ] );
		}							
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */

	float min = FindMin( min1d, _BLOCKNUM * _NT );
	float max = FindMax( max1d, _BLOCKNUM * _NT );
	printf("In %ld steps, min=%f, max=%f\n", _NT, min, max ); 

	double range = max - min;
	printf("1D DWT in X: NRMS = %e, NLMAX = %e\n", FindRMS( rmse1d, _BLOCKNUM * _NT) / range,
										    FindMax( lmax1d, _BLOCKNUM * _NT) / range);

	delete[] rmse1d;
	delete[] lmax1d;
	delete[] min1d;
	delete[] max1d;
}

void
Wavelet4D::YDimParallelExec()
{
	/*
  	 * RMS and LMAX with same block index but different _NT are stored together
	 */
	double* rmse1d = new double[ _BLOCKNUM * _NT ];
	double* lmax1d = new double[ _BLOCKNUM * _NT ];
	float*  min1d  = new float[ _BLOCKNUM * _NT ];
	float*  max1d  = new float[ _BLOCKNUM * _NT ];

	/*
	 * Each thread takes care of one index from _BLOCKNUM.
	 */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

			slices[t] -> GetMinMax( min1d[ i*_NT + t ], max1d[ i*_NT + t ] );

			/* individual slice reconstruction and evaluation */
			slices[t] -> DecomposeY();
			slices[t] -> ReconstructY (_cratio);
			slices[t] -> Evaluate( rmse1d[ i*_NT + t ], lmax1d[ i*_NT + t ] );
		}							
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */

	float min = FindMin( min1d, _BLOCKNUM * _NT );
	float max = FindMax( max1d, _BLOCKNUM * _NT );
	printf("In %ld steps, min=%f, max=%f\n", _NT, min, max ); 

	double range = max - min;
	printf("1D DWT in Y: NRMS = %e, NLMAX = %e\n", FindRMS( rmse1d, _BLOCKNUM * _NT) / range,
										    FindMax( lmax1d, _BLOCKNUM * _NT) / range);

	delete[] rmse1d;
	delete[] lmax1d;
	delete[] min1d;
	delete[] max1d;
}

void
Wavelet4D::ZDimParallelExec()
{
	/*
  	 * RMS and LMAX with same block index but different _NT are stored together
	 */
	double* rmse1d = new double[ _BLOCKNUM * _NT ];
	double* lmax1d = new double[ _BLOCKNUM * _NT ];
	float*  min1d  = new float[ _BLOCKNUM * _NT ];
	float*  max1d  = new float[ _BLOCKNUM * _NT ];

	/*
	 * Each thread takes care of one index from _BLOCKNUM.
	 */
	int nthreadSys = omp_get_max_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );

	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );

			slices[t] -> GetMinMax( min1d[ i*_NT + t ], max1d[ i*_NT + t ] );

			/* individual slice reconstruction and evaluation */
			slices[t] -> DecomposeZ();
			slices[t] -> ReconstructZ (_cratio);
			slices[t] -> Evaluate( rmse1d[ i*_NT + t ], lmax1d[ i*_NT + t ] );
		}							
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}
	/* OMP parallel section finish */

	float min = FindMin( min1d, _BLOCKNUM * _NT );
	float max = FindMax( max1d, _BLOCKNUM * _NT );
	printf("In %ld steps, min=%f, max=%f\n", _NT, min, max ); 

	double range = max - min;
	printf("1D DWT in Z: NRMS = %e, NLMAX = %e\n", FindRMS( rmse1d, _BLOCKNUM * _NT) / range,
										    FindMax( lmax1d, _BLOCKNUM * _NT) / range);

	delete[] rmse1d;
	delete[] lmax1d;
	delete[] min1d;
	delete[] max1d;
}



void
Wavelet4D::StartMonitor()
{
	long size_should = _NX * _NY * _NZ * 4;
	struct stat buffer;
	for( long i = 0; i < _NT; i++ )
	{
		const char* name = _filenames[i].c_str();
		int status;
		while(true)
		{
			status = stat(name, &buffer );
			if( status==0 && size_should==(long)buffer.st_size )
				break;
		}	

	}
}



void 
Wavelet4D::PutValToWrite( long x, long y, long z, long t, float v )
{
	assert( _reconstruct_buf    != NULL );
	assert( _reconstruct_buf[t] != NULL );
	long idx = z*_NX*_NY + y*_NX + x;
	_reconstruct_buf[t][idx] = v;
}
