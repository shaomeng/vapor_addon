#include "vapor/SamSlice2.h"

using namespace VAPoR;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

// Constructor
SamSlice2::SamSlice2( string filename, string wavename, 
                    size_t NX, size_t NY, size_t NZ )
{
    _filename = filename;
    _wavename = wavename;
    _NX = NX;
    _NY = NY;
    _NZ = NZ;
    _raw = NULL;
    _C = NULL;
    _L = NULL;
    _reconstructed = NULL;
    
    ReadFile();
    _mw = new MatWaveWavedec( wavename );
    _nlevels = min( min(_mw->wmaxlev(NX), _mw->wmaxlev(NY)), _mw->wmaxlev(NZ));
    _clen = _mw->coefflength3( NX, NY, NZ, _nlevels );
    assert( _clen == NX * NY * NZ );
}

void
SamSlice2::ReadFile(  )
{
    FILE* f = fopen( _filename.c_str(), "rb" );
    if( f != NULL ) {
        fseek( f, 0, SEEK_END );
        size_t size = ftell( f );

        size_t totallen = _NX*_NY*_NZ;
        assert (size == 4*_NX*_NY*_NZ );
        if( _raw == NULL )  _raw = new float[ totallen ];

        fseek( f, 0, SEEK_SET );
        size_t rsize = fread( _raw, sizeof(float), totallen, f);
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

SamSlice2::~SamSlice2()
{
    if( _mw )           delete _mw;
    if( _raw )          delete[] _raw;
    if( _C )            delete[] _C;
    if( _L )            delete[] _L;
    if( _reconstructed) delete[] _reconstructed;
}

int
SamSlice2::Decompose()
{
    if( _C == NULL )    _C = new float[ _clen ];
    if( _L == NULL )    _L = new size_t[ (21 * _nlevels) + 6 ];
    _mw -> computeL3( _NX, _NY, _NZ, _nlevels, _L );
    assert( _raw != NULL );
    int rc = _mw->wavedec3( _raw, _NX, _NY, _NZ, _nlevels, _C, _L );
    assert (rc >= 0);
    return rc;
}

int
SamSlice2::Reconstruct( int ratio )
{
    if( _reconstructed == NULL )
        _reconstructed = new float[ _NX * _NY * _NZ ];
    int rc = 0;

    if( ratio > 1 ) {   // cull coefficients
        float nth = FindCoeffThreshold(ratio);  // nth largest, indexing from 1.
        float nnth = -1.0 * nth;
        
        float* culledCoeffs = new float[ _clen ];

        for( size_t i = 0; i < _clen; i++ )
            if( (_C[i] >= nth || _C[i] <= nnth) ) {
                culledCoeffs[i] = _C[i];
            }
            else  culledCoeffs[i] = 0.0;

        rc = _mw -> waverec3( culledCoeffs, _L, _nlevels, _reconstructed );
        delete[] culledCoeffs;
    }
    else            // Decompress using the same coeffs, no compression
        rc = _mw -> waverec3( _C, _L, _nlevels, _reconstructed );

    return rc;
}


float
SamSlice2::FindCoeffThreshold( int ratio )
{
    size_t n = _clen  / ratio - 1;    // Find the nth largest, indexing from 1.

    vector<float> allCoeffs; 
    for( size_t i = 0; i < _clen; i++ )
        if ( _C[i] > 0 )
            allCoeffs.push_back( -1.0 * _C[i] );
        else
            allCoeffs.push_back( _C[i] );

    std::nth_element( allCoeffs.begin(), allCoeffs.begin()+n, allCoeffs.end() );
    float nth = -1.0 * (allCoeffs[n]);

    return nth;
}

void
SamSlice2::UpdateCoeffs( const float* update )
{
    if( _C == NULL )    _C = new float[ _clen ];
    memcpy( (void*)_C, (void*)update, sizeof(float) * _clen );
}
float*
SamSlice2::HandoverCoeffs(  )
{
    float* c = _C;
    _C = NULL;
    return c;
}

float*
SamSlice2::GetRawPtr()
{
    if( _raw != NULL )      return _raw;
    else{
        cerr << "Raw pointer is NULL when asked" << endl;
        return NULL;
    }
}
float*
SamSlice2::GetCoeffsPtr()
{
    if( _C != NULL )      return _C;
    else{
        cerr << "Coeffs pointer is NULL when asked" << endl;
        return NULL;
    }
}
float*
SamSlice2::GetReconstructedPtr()
{
    if( _reconstructed != NULL )      return _reconstructed;
    else{
        cerr << "Reconstructed pointer is NULL when asked" << endl;
        return NULL;
    }
}
void
SamSlice2::FreeRaw()
{
    if( _raw )              { delete[] _raw; _raw = NULL; }
}
void 
SamSlice2::FreeCoeffs()
{
    if( _C )           { delete[] _C; _C = NULL; }
}
void
SamSlice2::FreeReconstructed()
{
    if( _reconstructed )    { delete[] _reconstructed; _reconstructed = NULL; }
}

void
SamSlice2::Print10Coeffs()
{
    cerr << "here are the first 10 coeffs: " << endl;
    for( int i = 0; i < 10; i++ )
        cerr << "\t" << _C[i] << endl;
}

void
SamSlice2::Print10Raws()
{
    cerr << "here are the first 10 raw inputs: " << endl;
    for( int i = 0; i < 10; i++ )
        cerr << "\t" << _raw[i] << endl;
}
