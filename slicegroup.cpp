#include "slicegroup.h"

using namespace VAPoR;

SliceGroup::SliceGroup(string wavename )
{
    _wavename = wavename;
    _nslices = 0;
    _nlevels1d = 0;
    _sliceLen = 0;
    _buf = NULL;
    _mw = new MatWaveWavedec( _wavename );
//    _rms = 0.0;
//    _lmax = 0.0;
}

SliceGroup::~SliceGroup()
{
    if( _buf )                      delete[] _buf;
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
    _nlevels1d = _mw -> wmaxlev( _nslices );
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
        for( size_t i = 0; i < _sliceLen; i++ )
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
    size_t L1d[ _nlevels1d + 2 ];
    _mw -> computeL( _nslices, _nlevels1d, L1d );
    float nth = FindCoeffThreshold( ratio ); // use coeffs larger than nth.
    float nnth = -1.0 * nth;

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
 * Later I decide to put this function in the main file.
 * SliceGroup then only in charge of decompose and reconstruct.
 *
 *
void 
SliceGroup::Evaluate()
{
    double rmsArr[ _nslices ];
    double lmaxArr[ _nslices ];
    for( int i = 0; i < _nslices; i++ )
        _sliceVec[i]->Evaluate( rmsArr[i], lmaxArr[i] );

    double sum = 0.0;
    double c = 0.0;
    for( int i = 0; i < _nslices; i++ ) {
        double y = rmsArr[i] * rmsArr[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= double( _nslices );
    _rms = sqrt( sum );

    for( int i = 0; i < _nslices; i++ )
        if( lmaxArr[i] > _lmax )
            _lmax = lmaxArr[i];
}
 */
