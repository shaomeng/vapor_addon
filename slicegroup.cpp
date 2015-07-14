#include "slicegroup.h"

using namespace VAPoR;

SliceGroup::SliceGroup(string wavename )
{
    _wavename = wavename;
    _rawlen = 0;
    _nslices = 0;

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
SliceGroup::AddSlice( Cube3D* slice )
{
    _sliceVec.push_back( slice );
}

void
SliceGroup::Initialize( )
{
    _nslices = _sliceVec.size();
    _sliceLen = _sliceVec[0] -> GetCoeffLen();
    assert( _buf == NULL );
    _buf = new float[ _nslices * _sliceLen ];
    assert( _buf != NULL );
    
//    #pragma omp parallel for
    for( size_t i = 0; i < _sliceLen; i++ )
        for( size_t j = 0; j < _nslices; j++ )
            _buf[i*_nslices + j] = _sliceVec[j] -> GetCoeff(i);
}


void
SliceGroup::Decompose( )
{
//    #pragma omp parallel
    {
        MatWaveWavedec mv( _wavename );
        size_t l1d[ _nlevels1d+2 ];

//        #pragma omp for
        for( size_t i = 0; i < _slicelen; i++ )
        {
            float* src = _buf + i*_nslices;
            float dst[ _nslices ];
            int rc = mv.wavedec( src, _nslices, _nlevels1d, dst, l1d );
            assert (rc >= 0 );
            memcpy( src, dst, sizeof(float) * _nslices );
        }
    }
}

void
SliceGroup::Reconstruct( int ratio )
{
    float L1d[ _nlevels1d + 2 ];
    _mw -> computeL( _nslices, _nlevels1d, L1d );
    float nth = FindCoeffThreshold( ratio ); // use coeffs larger than nth.
    float nnth = -1.0 * nth;
    size_t totalC = _sliceLen * _nslices;

//    #pragma omp parallel
    {
        MatWaveWavedec mv( _wavename );
        float src[ _nslices ];

//        #pragma omp for
        for( size_t i = 0; i < _sliceLen; i++ )
        {
            for( size_t j = 0; j < _nslices; j++ )
            {
                float c = _buf[ i*_nslices + j ];
                if( c < nth && c > nnth )   src[j] = 0.0;
                else                        src[j] = c;
            }
            float* dst = _buf + i*_nslices;
            int rc = mv.waverec( src, L1d, _nlevels1d, dst );
            assert (rc >= 0 );
        }
    }
}

void
SliceGroup::UpdateSlices( )
{
    for( size_t i = 0; i < _sliceLen; i++ )
        for( size_t j = 0; j < _nslices; j++ )
            _sliceVec[j] -> PutCoeff( i, _buf[ i*_nslices+j ] );
}


float
SliceGroup::FindCoeffThreshold( int ratio )
{
    if( ratio < 1 ) {
        cerr << "SliceGroup::FindCoeffThreshold( int ratio ): ratio < 1 " << endl;
        return 0.0;
    }
    else if (ratio == 1)
        return 0.0;

    size_t nCoeffs = _nslices * _sliceLen;
    size_t n = nCoeffs / ratio - 1; // the nth largest, indexing from 1.

    vector<float> allCoeffs( nCoeffs, 0.0 );
//    #pragma omp parallel for
    for( size_t i = 0; i < nCoeffs; i++ )
            if ( _buf[i] > 0 )
                allCoeffs[i] = ( -1.0 * _buf[i] ); 
            else
                allCoeffs[i] = _buf[i]; 
    
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

