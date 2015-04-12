#include "vapor/SamSlice.h"

using namespace VAPoR;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

SamSlice::SamSlice()
{
    _raw = NULL;
    _rawlen = 0;
    _coeffs = NULL;
    _reconstructed = NULL;
    _c1 = NULL;
}

SamSlice::SamSlice( string filename )
{
    _filename = filename;
    _raw = NULL;
    _rawlen = 0;
    _coeffs = NULL;
    _reconstructed = NULL;
    _c1 = NULL;
    ReadFile( filename );
}

void
SamSlice::ReadFile( string filename )
{
    _filename = filename;
    FILE* f = fopen( filename.c_str(), "rb" );
    if( f != NULL ) {

        fseek( f, 0, SEEK_END );
        size_t size = ftell( f );
        _rawlen = size / 4;
        if( _raw == NULL )
            _raw = new float[ _rawlen ];

        fseek( f, 0, SEEK_SET );
        size_t rsize = fread( _raw, sizeof(float), _rawlen, f);
        fclose( f );
        if( rsize != _rawlen ) {
            cerr << "read size error: " << rsize << endl;
            exit(1);
        }
    }
    else{
        cerr << "file open error: " << filename << endl;
        exit(1);
    }
}

SamSlice::~SamSlice(){
    if( _c1 )           delete _c1;
    if( _raw )          delete[] _raw;
    if( _coeffs )       delete[] _coeffs;
    if( _reconstructed) delete[] _reconstructed;
}

int 
SamSlice::SetupCompressor1( vector< size_t > &dims, string wavename )
{
    size_t tmp = 1;
    for( int i = 0; i < dims.size(); i++ )
        tmp *= dims[i];
    assert( tmp == _rawlen );
    _dims = dims;
    if( wavename.compare("bior1.1") == 0 ) {
        _c1 = new Compressor( dims, wavename, "symh" );
        return 0;
    }
    else if( wavename.compare( "bior3.3" ) == 0) {
        _c1 = new Compressor( dims, wavename, "symh" );
        return 0;
    }
    else if( wavename.compare("bior4.4") == 0 ) {
        _c1 = new Compressor(dims, wavename, "symw" );
        return 0;
    }
    else{
        cerr << "wavename not recognized: " << wavename << endl;
        return 1;
    }
}

int 
SamSlice::Compress1()
{
    if( _coeffs == NULL )
        _coeffs = new float[ _rawlen ];
    vector< size_t > dst_arr_len;
    dst_arr_len.push_back( _rawlen );
    _sigmaps.resize( 1 );

    return _c1 -> Decompose( _raw, _coeffs, dst_arr_len, _sigmaps );
}

int 
SamSlice::Decompress1( int ratio )
{
    if( _rawlen  % ratio != 0 )
        cerr << "WARNING: requesting compression ratio is not divisible: "
             << ratio << endl;

    float nth = FindCoeffThreshold( ratio ); // nth largest, indexing from 1.

    if( _reconstructed == NULL )
        _reconstructed = new float[ _rawlen ];
    float* culledCoeffs = new float[ _rawlen ];
    size_t inCount = 0;

    for( size_t i = 0; i < _rawlen; i++ ) {
        if( _coeffs[i] >= nth || _coeffs[i] <= -1.0 * nth ) {
            culledCoeffs[i] = _coeffs[i];
            inCount++;
        }
        else
            culledCoeffs[i] = 0.0;
    }

    int rc = _c1 -> Reconstruct( culledCoeffs, _reconstructed, _sigmaps, -1 );
    
    delete[] culledCoeffs;

    size_t n = _rawlen / ratio; // num of coeffs to use
    if( n != inCount )
        cerr << "WARNING: should use " << n << " coeffs, but actually used : "
             << inCount << endl;

    return rc;
}

void
SamSlice::PutCoeffsInPosition()
{
    assert( _rawlen == _sigmaps[0].GetNumSignificant() );
    float* coeffsInPosition = new float[ _rawlen ];

    size_t idx;
    _sigmaps[0].GetNextEntryRestart();
    for( size_t i = 0; i < _rawlen; i++ ) {
        _sigmaps[0].GetNextEntry( &idx );
        coeffsInPosition[ idx ] = _coeffs[i];
    }

/*
    for( size_t i = 0; i < _rawlen; i++ ) {
        size_t idx;
        _sigmaps[0].GetCoordinates(i, &idx);
        coeffsInPosition[ idx ] = _coeffs[i];
    }
*/

    _sigmaps[0].Clear();
    for( size_t i = 0; i < _rawlen; i++ )
        _sigmaps[0].Set( i );

    delete[] _coeffs;
    _coeffs = coeffsInPosition;
}

float
SamSlice::FindCoeffThreshold( int ratio )
{
    size_t n = _rawlen  / ratio - 1;    // Find the nth largest, indexing from 1.

    vector<float>* allCoeffs = new vector<float>( );
    for( size_t i = 0; i < _rawlen; i++ )
        if ( _coeffs[i] > 0 )
            allCoeffs -> push_back( -1.0 * _coeffs[i] );
        else
            allCoeffs -> push_back( _coeffs[i] );

    std::nth_element( allCoeffs->begin(), allCoeffs->begin()+n, allCoeffs->end() );
    float nth = -1.0 * (allCoeffs -> at(n));
    delete allCoeffs;

    return nth;
}

void
SamSlice::ReplaceCoeffs( float* replace )
{
    for( size_t i = 0; i < _rawlen; i++ )
        _coeffs[i] = replace[i];
}

void
SamSlice::Print10Coeffs()
{
    for( int i = 0; i < 10; i++ )
        cerr << "\t" << _coeffs[i] << endl;
}

void
SamSlice::Print10Raws()
{
    for( int i = 0; i < 10; i++ )
        cerr << "\t" << _raw[i] << endl;
}
