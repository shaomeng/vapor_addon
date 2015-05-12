#include "vapor/SamSliceGroup2.h"

using namespace VAPoR;

SamSliceGroup2::SamSliceGroup2(string wavename, size_t rawlen, size_t nslices)
{
    _wavename = wavename;
    _rawlen = rawlen;
    _nslices = nslices;

    _raw.resize( _nslices, NULL );
    _C1d.resize( _rawlen,  NULL );
    _reconstructed.resize( _nslices, NULL );
    _L1d = NULL;
    
    _mw = new MatWaveWavedec( _wavename );
    _nlevels1d = _mw -> wmaxlev( _nslices );
    _clen1d    = _mw -> coefflength( _nslices, _nlevels1d );
    assert( _clen1d = _nslices );
}

SamSliceGroup2::~SamSliceGroup2()
{
    for( size_t i = 0; i < _raw.size(); i++ )
        if( _raw[i] )           delete[] _raw[i];
    for( size_t i = 0; i < _C1d.size(); i++ )
        if( _C1d[i] )            delete[] _C1d[i];
    for( size_t i = 0; i < _reconstructed.size(); i++ )
        if( _reconstructed[i] )     delete[] _reconstructed[i];

    if( _L1d )                      delete[] _L1d; 
    if( _mw )                       delete _mw;
}

void
SamSliceGroup2::UpdateRawPtr( size_t i, float* ptr )
{
    if( _raw[i] ) {
        cerr << "WARNING! raw pointer of SliceGroup is not NULL when updating" << endl;
        delete[] _raw[i];
    }
    _raw[i] = ptr;
}

void
SamSliceGroup2::Decompose( )
{
    for( size_t i = 0; i < _raw.size(); i++ )
        assert( _raw[i] );
    for( int i = 0; i < _C1d.size(); i++ )
        if( _C1d[i] ) {
            delete[] _C1d[i];
            _C1d[i] = NULL;
        }
    if( _L1d == NULL )     _L1d = new size_t[ _nlevels1d + 2 ];
    _mw -> computeL( _nslices, _nlevels1d, _L1d );

/* 
 * Serial version, works great
 */
    int rc;
    float src[ _nslices ];
    #pragma unroll
    for( size_t i = 0; i < _rawlen; i++ )
    {
        for( size_t j = 0; j < _nslices; j++ )      src[j] = _raw[j][i];
        float* dst = new float[ _clen1d ];
        rc = _mw -> wavedec( src, _nslices, _nlevels1d, dst, _L1d );
        assert (rc >= 0);
        _C1d[i] = dst;
    }

/*
 * OpenMP version
 *
    #pragma omp parallel 
    {
        MatWaveWavedec mw( _wavename );
        int rc;
        #pragma omp for 
        for( size_t i = 0; i < _rawlen; i++ )
        {
            float src[ _nslices ]; 
            for( int j = 0; j < _nslices; j++ )         src[j] = _raw[j][i];
            size_t l1d[ _nlevels1d + 2 ];
            float* dst = new float[ _clen1d ];
            rc = mw.wavedec( src, _nslices, _nlevels1d, dst, l1d );
            assert (rc >= 0);
            _C1d[i] = dst;
    }
    }
*/
}

void
SamSliceGroup2::Reconstruct( int ratio )
{

    float nth = 0.0;
    if( ratio > 1 )
        nth = FindCoeffThreshold( ratio ); // use coeffs larger than nth.
    float nnth = -1.0 * nth;


    for( int i = 0; i < _nslices; i++ ) 
        if( _reconstructed[i] == NULL )
            _reconstructed[i] = new float[ _rawlen ];
    

    size_t inCount = 0;
    
/*
 * Serial Version, works great!
 */
    float src[ _clen1d ], dst[ _nslices ];  // they should have the same length
    int rc;
    #pragma unroll
    for( size_t i = 0; i < _rawlen; i++ )
    {    
        for( int j = 0; j < _clen1d; j++ ) {
            float c = _C1d[i][j];
            if( (c >= nth || c <= nnth) ) {
                src[j] = c;
                inCount++;
            }
            else          src[j] = 0.0;
        }

        rc = _mw -> waverec( src, _L1d, _nlevels1d, dst );
        assert( rc >= 0 );
        for( int j = 0; j < _nslices; j++ )
            _reconstructed[j][i] = dst[j];        
    }   

/*
 * OpenMP Version
 * The final result is different from the serial one.
 * Not sure where the problem is
 *
    #pragma omp parallel
    {
        int rc;
        float src[ _clen1d ], dst[ _nslices ];  // they should have the same length
        MatWaveWavedec mv( _wavename );
        #pragma omp for
        for( size_t i = 0; i < _rawlen; i++ )
        {
            for( int j = 0; j < _clen1d; j++ )
            {
                float c = _C1d[i][j];
                if( (c >= nth || c <= nnth) && inCount <= nc) {
                    src[j] = c;
                    #pragma omp atomic
                    inCount++;
                }
                else          src[j] = 0.0;
            }

            rc = mv.waverec( src, _L1d, _nlevels1d, dst );
            assert( rc >= 0 );
            for( int j = 0; j < _nslices; j++ )
                _reconstructed[j][i] = dst[j];        
        }
    }
*/

}

float*
SamSliceGroup2::GetReconstructedPtr( int i )
{
    if( i < _reconstructed.size() )
        return _reconstructed[i];
    else{
        cerr << "SamSliceGroup2::GetDecompressedCoeff(), idx too large! " << endl;
        return NULL;
    }
}


float
SamSliceGroup2::FindCoeffThreshold( int ratio )
{
    size_t nCoeffs = _nslices * _rawlen;
    size_t n = nCoeffs / ratio - 1; // the nth largest, indexing from 1.


    vector<float>* allCoeffs = new vector<float>( );
    for( size_t i = 0; i < _C1d.size(); i++ )
        for( int j = 0; j < _nslices; j++ ) 
            if ( _C1d[i][j] > 0 )
                allCoeffs -> push_back( -1.0 * _C1d[i][j] ); 
            else
                allCoeffs -> push_back( _C1d[i][j] ); 
    
    if( allCoeffs->size() != nCoeffs )
        cerr << "WARNING: coeffs number don't match in SamSliceGroup2::FindCoeffThreshold()!" << endl;

    std::nth_element( allCoeffs->begin(), allCoeffs->begin()+n, allCoeffs->end() );
    float nth = -1.0 * (allCoeffs -> at(n));

    delete allCoeffs;

    return nth;
}


/*
void
SamSliceGroup2::FreeReconstructed( int i )
{
    assert( i < _reconstructed.size() );
    if( _reconstructed[i] ) {
            delete[] _reconstructed[i];
            _reconstructed[i] = NULL;
    }
}

void
SamSliceGroup2::FreeRaw( int i )
{
    assert( i < _raw.size() );
    if( _raw[i] ) {
            delete[] _raw[i];
            _raw[i] = NULL;
    }
}

void
SamSliceGroup2::FreeCoeffs()
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
SamSliceGroup2::Print1DRaw()
{
    for( int i = 0; i < _nslices; i++ )
        cerr << "\t" << _raw[i][0] << endl;
}


void
SamSliceGroup2::Print1DCoeffs()
{
    for( int i = 0; i < _nslices; i++ )
        cerr << "\t" << _coeffs[0][i] << endl;
}
*/

