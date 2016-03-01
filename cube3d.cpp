#include "cube3d.h"

using std::vector;
using std::string;
using std::cerr;
using std::endl;

// Constructor for reading a whole file as current cube to work on
Cube3D::Cube3D( string filename, string wavename, 
                    size_t NX, size_t NY, size_t NZ )
{
    _filename = filename;
    _wavename = wavename;
    _NX = _NX_total = _endX = NX;
    _NY = _NY_total = _endY = NY;
    _NZ = _NZ_total = _endZ = NZ;
    _L = NULL;

	_startX = _startY = _startZ = 0;
    
    _mw = new VAPoR::MatWaveWavedec( wavename );
    _nlevels = min( min(_mw->wmaxlev(NX), _mw->wmaxlev(NY)), _mw->wmaxlev(NZ));
    _clen = _mw->coefflength3( NX, NY, NZ, _nlevels );
    assert( _clen == NX * NY * NZ );
    _C = new float[ _clen ];
    ReadFileChunck( _C );
}

// Constructor for reading a block out of a big chunck as current cube.
Cube3D::Cube3D( string filename, string wavename, 
				size_t NX, size_t NY, size_t NZ,
				size_t NX_total, size_t NY_total, size_t NZ_total,
				size_t startX, size_t endX,
				size_t startY, size_t endY,
				size_t startZ, size_t endZ )
{
    _filename = filename;
    _wavename = wavename;
    _NX = NX;
    _NY = NY;
    _NZ = NZ;
    _L = NULL;

	_NX_total = NX_total;
	_NY_total = NY_total;
	_NZ_total = NZ_total;
	_startX = startX;
	_endX   = endX;
	_startY = startY;
	_endY   = endY;
	_startZ = startZ;
	_endZ   = endZ;
	
    _mw = new VAPoR::MatWaveWavedec( wavename );
    _nlevels = min( min(_mw->wmaxlev(NX), _mw->wmaxlev(NY)), _mw->wmaxlev(NZ));
    _clen = _mw->coefflength3( NX, NY, NZ, _nlevels );
    assert( _clen == NX * NY * NZ );
    _C = new float[ _clen ];
    ReadFileChunck( _C );
}

// Reads the whole file into buffer
/*
void
Cube3D::ReadFile( float* buf )
{
    FILE* f = fopen( _filename.c_str(), "rb" );
    if( f != NULL ) {
        fseek( f, 0, SEEK_END );
        size_t size = ftell( f );

        size_t totallen = _NX_total * _NY_total * _NZ_total;
        assert (size == sizeof(float) * totallen );
        assert ( buf != NULL );

        fseek( f, 0, SEEK_SET );
        size_t rsize = fread( buf, sizeof(float), totallen, f);
        fclose( f );
        if( rsize != totallen ) {
            cerr << "read size error: " << rsize << endl;
            exit(1);
        }
    }
    else{
        cerr << "file open error: " << _filename << endl;
        exit(1);
    }
}
*/

// Reads a chunck of file specified by ranges
void
Cube3D::ReadFileChunck( float* buf )
{
    FILE* f = fopen( _filename.c_str(), "rb" );
    if( f != NULL ) {
        fseek( f, 0, SEEK_END );
        size_t size = ftell( f );

        size_t totallen = _NX_total * _NY_total * _NZ_total;
        assert (size == sizeof(float) * totallen );
        assert ( buf != NULL );
		float* tmp = new float[ _NX ];

		size_t counter = 0;
		for( size_t k = _startZ; k < _endZ; k++ )
			for( size_t j = _startY; j < _endY; j++ )
			{
				size_t offset = k*_NX_total*_NY_total + j*_NX_total + _startX;
				int rt = fseek( f, sizeof(float)*offset, SEEK_SET );
				assert ( rt == 0 );
				size_t rsize = fread( tmp, sizeof(float), _NX, f );
				assert( rsize == _NX );
				memcpy( (void*)(buf + _NX * counter), (void*)tmp, 
						sizeof(float) * _NX );
				counter++;
			}
		delete[] tmp;
        fclose( f );
    }
    else{
        cerr << "file open error: " << _filename << endl;
        exit(1);
    }
}

Cube3D::~Cube3D()
{
    if( _mw )           delete _mw;
    if( _C )            delete[] _C;
    if( _L )            delete[] _L;
}

int
Cube3D::Decompose()
{
    float* buf = new float[ _clen ];
    assert( buf != NULL );
    if( _L == NULL )    _L = new size_t[ (21 * _nlevels) + 6 ];
    _mw -> computeL3( _NX, _NY, _NZ, _nlevels, _L );
    int rc = _mw->wavedec3( _C, _NX, _NY, _NZ, _nlevels, buf, _L );
    assert (rc >= 0);

    delete[] _C;
    _C = buf;

    return rc;
}

size_t
Cube3D::xyz2idx( size_t x, size_t y, size_t z )
{
	return (z * _NX * _NY + y * _NX + x);
}

void
Cube3D::DecomposeX()
{
	VAPoR::MatWaveWavedec mw1d( _wavename );
	size_t nlevels1d = mw1d.wmaxlev( _NX );
	size_t l1d[ nlevels1d+2 ];

	float  src[ _NX ];
	float* buf = new float[ _clen ];

	for( size_t z = 0; z < _NZ; z++ )
		for( size_t y = 0; y < _NY; y++ )
		{
			for( size_t x = 0; x < _NX; x++ )
				src[x] = _C[ xyz2idx(x, y, z) ];

			float* dst = buf + xyz2idx(0, y, z);

			int rc = mw1d.wavedec( src, _NX, nlevels1d, dst, l1d );
			assert( rc >= 0 );
		}

    delete[] _C;
    _C = buf;
}

void
Cube3D::ReconstructX( int ratio )
{
    if( ratio > 1 ) {   						// cull coefficients
        float nth = FindCoeffThreshold(ratio);  // nth largest, indexing from 1.
        CullCoeffs( nth );
    }
	VAPoR::MatWaveWavedec mw1d( _wavename );
	size_t nlevels1d = mw1d.wmaxlev( _NX );
	size_t l1d[ nlevels1d+2 ];
	mw1d.computeL( _NX, nlevels1d, l1d );

    float* buf = new float[ _NX * _NY * _NZ ];
	float  dst[ _NX ];

	for( size_t z = 0; z < _NZ; z++ )
		for( size_t y = 0; y < _NY; y++ )
		{
			float* src = _C + xyz2idx(0, y, z);
			int rc = mw1d.waverec( src, l1d, nlevels1d, dst );
			assert( rc >= 0 );

			for( size_t x = 0; x < _NX; x++ )
				buf[ xyz2idx(x, y, z) ] = dst[x];
		}

    delete[] _C;
    _C = buf;
}


void
Cube3D::DecomposeY()
{
	VAPoR::MatWaveWavedec mw1d( _wavename );
	size_t nlevels1d = mw1d.wmaxlev( _NY );
	size_t l1d[ nlevels1d+2 ];

	float  src[ _NY ];
	float  dst[ _NY ];
	float* buf = new float[ _clen ];

	for( size_t z = 0; z < _NZ; z++ )
		for( size_t x = 0; x < _NX; x++ )
		{
			for( size_t y = 0; y < _NY; y++ )
				src[y] = _C[ xyz2idx(x, y, z) ];

			int rc = mw1d.wavedec( src, _NY, nlevels1d, dst, l1d );
			assert( rc >= 0 );

			for( size_t y = 0; y < _NY; y++ )
				buf[ xyz2idx(x, y, z) ] = dst[y];
		}

    delete[] _C;
    _C = buf;
}

void
Cube3D::ReconstructY( int ratio )
{
    if( ratio > 1 ) {   						// cull coefficients
        float nth = FindCoeffThreshold(ratio);  // nth largest, indexing from 1.
        CullCoeffs( nth );
    }
	VAPoR::MatWaveWavedec mw1d( _wavename );
	size_t nlevels1d = mw1d.wmaxlev( _NY );
	size_t l1d[ nlevels1d+2 ];
	mw1d.computeL( _NY, nlevels1d, l1d );

    float* buf = new float[ _NX * _NY * _NZ ];
	float  src[ _NY ];
	float  dst[ _NY ];

	for( size_t z = 0; z < _NZ; z++ )
		for( size_t x = 0; x < _NX; x++ )
		{
			for( size_t y = 0; y < _NY; y++ )
				src[y] = _C[ xyz2idx(x, y, z) ];
	
			int rc = mw1d.waverec( src, l1d, nlevels1d, dst );
			assert( rc >= 0 );

			for( size_t y = 0; y < _NY; y++ )
				buf[ xyz2idx(x, y, z) ] = dst[y];
		}

    delete[] _C;
    _C = buf;
}

void
Cube3D::DecomposeZ()
{
	VAPoR::MatWaveWavedec mw1d( _wavename );
	size_t nlevels1d = mw1d.wmaxlev( _NZ );
	size_t l1d[ nlevels1d+2 ];

	float  src[ _NZ ];
	float  dst[ _NZ ];
	float* buf = new float[ _clen ];

	for( size_t y = 0; y < _NY; y++ )
		for( size_t x = 0; x < _NX; x++ )
		{
			for( size_t z = 0; z < _NZ; z++ )
				src[z] = _C[ xyz2idx(x, y, z) ];

			int rc = mw1d.wavedec( src, _NZ, nlevels1d, dst, l1d );
			assert( rc >= 0 );

			for( size_t z = 0; z < _NZ; z++ )
				buf[ xyz2idx(x, y, z) ] = dst[z];
		}

    delete[] _C;
    _C = buf;
}

void
Cube3D::ReconstructZ( int ratio )
{
    if( ratio > 1 ) {   						// cull coefficients
        float nth = FindCoeffThreshold(ratio);  // nth largest, indexing from 1.
        CullCoeffs( nth );
    }
	VAPoR::MatWaveWavedec mw1d( _wavename );
	size_t nlevels1d = mw1d.wmaxlev( _NZ );
	size_t l1d[ nlevels1d+2 ];
	mw1d.computeL( _NZ, nlevels1d, l1d );

    float* buf = new float[ _NX * _NY * _NZ ];
	float  src[ _NZ ];
	float  dst[ _NZ ];

	for( size_t y = 0; y < _NY; y++ )
		for( size_t x = 0; x < _NX; x++ )
		{
			for( size_t z = 0; z < _NZ; z++ )
				src[z] = _C[ xyz2idx(x, y, z) ];
	
			int rc = mw1d.waverec( src, l1d, nlevels1d, dst );
			assert( rc >= 0 );

			for( size_t z = 0; z < _NZ; z++ )
				buf[ xyz2idx(x, y, z) ] = dst[z];
		}

    delete[] _C;
    _C = buf;
}


int
Cube3D::Reconstruct( int ratio )
{
    if( ratio > 1 ) {   // cull coefficients
        float nth = FindCoeffThreshold(ratio);  // nth largest, indexing from 1.
        CullCoeffs( nth );
    }
    float* buf = new float[ _NX * _NY * _NZ ];
    int rc = _mw -> waverec3( _C, _L, _nlevels, buf );

    delete[] _C;
    _C = buf;

    return rc;
}

float
Cube3D::GetCoeff( size_t idx )
{
    if( _C == NULL ) {
        cerr << "\tCube3D::GetCoeff( int idx ): _C == NULL!" << endl;
        exit (1);
    }
    else if( idx >= _clen ) {
        cerr << "\tCube3D::GetCoeff( int idx ): idx out of range:" << endl;
        exit (1);
    }
    
    return _C[idx];
}

void 
Cube3D::PutCoeff( size_t idx, float c )
{
    if( _C == NULL ) {
        cerr << "\tCube3D::PutCoeff( int idx, float c ): _C == NULL!" << endl;
        exit (1);
    }
    else if( idx >= _clen ) {
        cerr << "\tCube3D::PutCoeff( int idx, float c ): idx out of range:" << endl;
        exit (1);
    }

    _C[idx] = c;
}

void
Cube3D::CullCoeffs( float t )
{
    if( t < 0 )
        cerr << "\tCube3D::CullCoeffs( float t ): t < 0!" << endl;
    else if( _C == NULL )
        cerr << "\tCube3D::CullCoeffs( int idx ): _C == NULL!" << endl;
    else{
        float nt = -1.0 * t;
        for( size_t i = 0; i < _clen; i++ )
            if( _C[i] < t && _C[i] > nt )
                _C[i] = 0.0;
    }
}

void 
Cube3D::Evaluate( double &rms, double &lmax )
{
    float* raw = new float[ _clen ];
    ReadFileChunck( raw );

    double sum = 0.0;
    double c = 0.0;
    double max = 0.0;
    double tmp;
    for( size_t i = 0; i < _clen; i++) {
        tmp = (double)raw[i] - (double)_C[i];
        if (tmp < 0)        tmp *= -1.0;
        if (tmp > max)      max = tmp;
        double y = tmp * tmp - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= double(_clen);
    sum = sqrt( sum );
    
    rms = sum;
    lmax = max;

	delete[] raw;
}

void 
Cube3D::EvaluateWithAnotherFile( string anotherFileName, long &numNonSpecial,
								 double &squaresum, double &lmax )
{
	/* temporarily redirect filename of current cube */
	string origFilename = _filename;
	_filename = anotherFileName;

	/* perform evaluation excluding special values */
    float* raw = new float[ _clen ];
    ReadFileChunck( raw );

    double sum = 0.0;
    double c = 0.0;
    double max = 0.0;
    double tmp;
	long no_special_count = 0;
    for( size_t i = 0; i < _clen; i++) 
		if( raw[i] < 1e+34 )
		{
			tmp = (double)raw[i] - (double)_C[i];
			double y = tmp * tmp - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
			if (fabs(tmp) > max)      max = fabs(tmp);
			no_special_count++;
		}
    
	numNonSpecial = no_special_count;
    squaresum = sum;
    lmax = max;

	delete[] raw;

	/* restore the original filename */
	_filename = origFilename;
}

void
Cube3D::GetMinMax( float &min, float &max )
{
	assert( _C != NULL );
	assert( _clen > 0 );
	float lmin = _C[0];
	float lmax = _C[0];
	for( size_t i = 1; i < _clen; i++ )
	{
		if( _C[i] > lmax )		lmax = _C[i];
		if( _C[i] < lmin )		lmin = _C[i];
	}
	min = lmin;
	max = lmax;
}

float
Cube3D::FindCoeffThreshold( int ratio )
{
    size_t n = _clen  / ratio - 1;    // Find the nth largest, indexing from 1.

    vector<float> allCoeffs( _clen, 0.0); 
    for( size_t i = 0; i < _clen; i++ )
        if ( _C[i] > 0 )
            allCoeffs[i] = ( -1.0 * _C[i] );
        else
            allCoeffs[i] = _C[i];

    std::nth_element( allCoeffs.begin(), allCoeffs.begin()+n, allCoeffs.end() );
    float nth = -1.0 * (allCoeffs[n]);

    return nth;
}

void
Cube3D::Print10Elements()
{
    cerr << "here are the first 10 elements: " << endl;
    for( int i = 0; i < 10; i++ )
        cerr << "\t" << _C[i] << endl;
/*
	for( int k = _startZ; k < _startZ + 5; k++ )
		for( int j = _startY+5; j < _startY+10; j++ )
			for( int i = _startX + 10; i < _startX + 15; i++ )
			{
				printf("V\(%d, %d, %d\) should be %d, ", i, j, k, i+j+k);
				int idx = _NX*_NY*(k-_startZ) + _NY*(j-_startY) + (i-_startX);
				cout << "It actually is: " << _C[ idx ] << endl;
			}
*/
}
	

int main( int argc, char* argv[] )
{
	string filename = argv[1];
	int cratio = atoi( argv[2] );

	int NX = 40;
	int NY = 40;
	int NZ = 40;

	int NX_total = 480;
	int NY_total = 480;
	int NZ_total = 280;

	int startX = 40;
	int endX   = 80;
	int startY = 40;
	int endY   = 80;
	int startZ = 40;
	int endZ   = 80;

	Cube3D* slice = new Cube3D( filename, "bior4.4", NX, NY, NZ, 
								NX_total, NY_total, NZ_total,
								startX, endX, startY, endY, startZ, endZ );

	slice->Print10Elements();
	
	slice->DecomposeZ();
	slice->ReconstructZ( cratio );

	slice->Print10Elements();
	
	double rms, lmax;
	slice->Evaluate( rms, lmax );
	printf("rms = %e, lmax=%e\n", rms, lmax );


	delete slice;
}
