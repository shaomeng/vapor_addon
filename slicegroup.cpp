#include "vapor/SamSliceGroup3.h"

using namespace VAPoR;

SamSliceGroup3::SamSliceGroup3(string wavename, size_t rawlen, size_t nslices)
{
    _wavename = wavename;
    _rawlen = rawlen;
    _nslices = nslices;

    _buf = NULL;
    _L1d = NULL;
    _C1d = NULL;
    
    _mw = new MatWaveWavedec( _wavename );
    _nlevels1d = _mw -> wmaxlev( _nslices );
    _clen1d    = _mw -> coefflength( _nslices, _nlevels1d );
    assert( _clen1d = _nslices );
    if( _clen1d != _nslices ) {
        cerr << "SamSliceGroup3 error: clen1d != nslices " << endl;
        exit(1);
    }
}

SamSliceGroup3::~SamSliceGroup3()
{
    if( _buf )                      delete[] _buf;
    if( _L1d )                      delete[] _L1d; 
    if( _C1d )                      delete[] _C1d; 
    if( _mw )                       delete _mw;
}

void
SamSliceGroup3::UpdateRawPtr( vector< float* > &rawarr )
{
    assert( rawarr.size() == _nslices );
    if( _buf == NULL ) {
        _buf = new float[ _rawlen * _nslices ];
    }
    if( _buf == NULL ) {
        cerr << "SamSliceGroup3 memory allocation failed!" << endl;
        exit(1);
    }

    #pragma omp parallel for
    for( size_t i = 0; i < _rawlen; i++ )
        for( size_t j = 0; j < _nslices; j++ )
            _buf[i*_nslices + j] = rawarr[j][i];
}

void
SamSliceGroup3::Decompose( )
{
    if( _L1d == NULL )     _L1d = new size_t[ _nlevels1d + 2 ];
    _mw -> computeL( _nslices, _nlevels1d, _L1d );
    if( _C1d == NULL )      _C1d = new float[ _rawlen * _nslices ];

    #pragma omp parallel
    {
        MatWaveWavedec mv( _wavename );
        size_t l1d[ _nlevels1d+2 ];

        #pragma omp for
        for( size_t i = 0; i < _rawlen; i++ )
        {
            float* src = _buf + i*_nslices;
            float* dst = _C1d + i*_nslices;
            int rc = mv.wavedec( src, _nslices, _nlevels1d, dst, l1d );
            assert (rc >= 0 );
        }
    }

}

void
SamSliceGroup3::Reconstruct( int ratio )
{

    float nth = 0.0;
    if( ratio > 1 )
        nth = FindCoeffThreshold( ratio ); // use coeffs larger than nth.
    float nnth = -1.0 * nth;
    size_t totalC = _rawlen * _clen1d;


    #pragma omp parallel
    {
        MatWaveWavedec mv( _wavename );
        float src[ _nslices ];

        #pragma omp for
        for( size_t i = 0; i < _rawlen; i++ )
        {
            for( size_t j = 0; j < _nslices; j++ )
            {
                float c = _C1d[ i*_nslices + j ];
                if( c < nth && c > nnth )   src[j] = 0.0;
                else                        src[j] = c;
            }
            float* dst = _buf + i*_nslices;
            int rc = mv.waverec( src, _L1d, _nlevels1d, dst );
            assert (rc >= 0 );
        }
        
    }
}

void
SamSliceGroup3::FillReconstructedPtr( int idx, float* ptr )
{
    assert( ptr != NULL );
    #pragma omp parallel for
    for( size_t i = 0; i < _rawlen; i++ )
        ptr[i] = _buf[ i*_nslices + idx ];
}


float
SamSliceGroup3::FindCoeffThreshold( int ratio )
{
    size_t nCoeffs = _nslices * _rawlen;
    size_t n = nCoeffs / ratio - 1; // the nth largest, indexing from 1.


    vector<float> allCoeffs( nCoeffs, 0.0 );
    #pragma omp parallel for
    for( size_t i = 0; i < nCoeffs; i++ )
            if ( _C1d[i] > 0 )
                allCoeffs[i] = ( -1.0 * _C1d[i] ); 
            else
                allCoeffs[i] = _C1d[i]; 
    
    std::nth_element( allCoeffs.begin(), allCoeffs.begin()+n, allCoeffs.end() );
    float nth = -1.0 * allCoeffs[n];


    return nth;
}


/*
void
SamSliceGroup3::FreeReconstructed( int i )
{
    assert( i < _reconstructed.size() );
    if( _reconstructed[i] ) {
            delete[] _reconstructed[i];
            _reconstructed[i] = NULL;
    }
}

void
SamSliceGroup3::FreeRaw( int i )
{
    assert( i < _raw.size() );
    if( _raw[i] ) {
            delete[] _raw[i];
            _raw[i] = NULL;
    }
}

void
SamSliceGroup3::FreeCoeffs()
{
    for( size_t i = 0; i < _coeffs.size(); i++ )
        if( _coeffs[i] ) {
            delete[] _coeffs[i];
            _coeffs[i] = NULL;
        }
    if( _sigmapGroup ) {
        delete[] _sigmapGroup;
        _sigmapGroup = NULL;
    }
}

void
SamSliceGroup3::Print1DRaw()
{
    for( int i = 0; i < _nslices; i++ )
        cerr << "\t" << _raw[i][0] << endl;
}


void
SamSliceGroup3::Print1DCoeffs()
{
    for( int i = 0; i < _nslices; i++ )
        cerr << "\t" << _coeffs[0][i] << endl;
}
*/

