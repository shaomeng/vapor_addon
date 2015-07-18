#include "cube3d.h"

using namespace VAPoR;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

// Constructor
Cube3D::Cube3D( string filename, string wavename, 
                    size_t NX, size_t NY, size_t NZ )
{
    _filename = filename;
    _wavename = wavename;
    _NX = NX;
    _NY = NY;
    _NZ = NZ;
    _L = NULL;
    
    _mw = new MatWaveWavedec( wavename );
    _nlevels = min( min(_mw->wmaxlev(NX), _mw->wmaxlev(NY)), _mw->wmaxlev(NZ));
    _clen = _mw->coefflength3( NX, NY, NZ, _nlevels );
    assert( _clen == NX * NY * NZ );
    _C = new float[ _clen ];
    ReadFile( _C );
}

void
Cube3D::ReadFile( float* buf )
{
    FILE* f = fopen( _filename.c_str(), "rb" );
    if( f != NULL ) {
        fseek( f, 0, SEEK_END );
        size_t size = ftell( f );

        size_t totallen = _NX*_NY*_NZ;
        assert (size == 4*_NX*_NY*_NZ );
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
Cube3D::Evaluate( double& rms, double& lmax )
{
    float* raw = new float[ _clen ];
    ReadFile( raw );

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
}

