#include "vapor/SamSliceGroup.h"

using namespace VAPoR;

SamSliceGroup::SamSliceGroup()
{
    _c1 = NULL;
}

void
SamSliceGroup::Setup( int nfiles, size_t rawlen, 
              const vector< float* > &raw )
{
    _nslices = nfiles;
    _rawlen = rawlen;
    _raw = raw;

    if( raw.size() != _nslices ){
        cerr << "ERROR: SamSliceGroup::SamSliceGroup() expecting " << _nslices
             << " coefficient arrays, but got " << raw.size() << endl;
        exit(1);
    }

    _coeffs.resize( rawlen, NULL);
    _reconstructed.resize( nfiles, NULL );
}

SamSliceGroup::~SamSliceGroup()
{
    for( int i = 0; i < _coeffs.size(); i++ )
        if( _coeffs[i] )      delete[] _coeffs[i];
    for( int i = 0; i < _reconstructed.size(); i++ )
        if( _reconstructed[i] )    delete[] _reconstructed[i];

    if( _c1 )                delete _c1; 
}

int
SamSliceGroup::Compress1( )
{
    for( int i = 0; i < _coeffs.size(); i++ )
        if( _coeffs[i] )      delete[] _coeffs[i];
    vector< size_t > dst_arr_len;
    dst_arr_len.push_back( _nslices );
    _sigmapGroup.clear();
    _sigmapGroup.resize( _rawlen );
    for( size_t i = 0; i < _rawlen; i++ )
        _sigmapGroup[i].resize( 1 );

    int rc;
    float src[ _nslices ];
    for( size_t i = 0; i < _rawlen; i++ ){

        for( int j = 0; j < _nslices; j++) 
            src[j] = _raw[j][i];
        
        float* dst = new float[ _nslices ];
        rc = _c1 -> Decompose( src, dst, dst_arr_len, _sigmapGroup[i] );
        if( rc != 0 ) {
            cerr << "Decompose() on time intervals has error: " << rc << endl;
            return rc;
        }
        _coeffs[i] = dst;
    }

    return 0;
}

int
SamSliceGroup::Decompress1( int ratio )
{
    size_t nCoeffs = _nslices * _rawlen;
//    if( nCoeffs % ratio != 0 )
//        cerr << "WARNING: requesting compression ratio is not divisible: " 
//             << ratio << endl;
    float nth = FindCoeffThreshold( ratio ); // use coeffs larger than nth.


    for( int i = 0; i < _nslices; i++ ) {
        if( _reconstructed[i] )     delete[] _reconstructed[i];
        _reconstructed[i] = new float[ _rawlen ];
    }

    float src[ _nslices], dst[ _nslices ];
    int rc;
    size_t inCount = 0;
    
    for( size_t i = 0; i < _rawlen; i++ ) {
        
        for( int j = 0; j < _nslices; j++ ) {
            float c = _coeffs[i][j];
            if( c >= nth || c <= -1.0 * nth ) {
                src[j] = c;
                inCount++;
            }
            else 
                src[j] = 0.0;
        }

        rc = _c1 -> Reconstruct( src, dst, _sigmapGroup[i], -1 );
        if( rc != 0 ) {
            cerr << "Reconstruct() error with code: " << rc << endl;
            return 1;
        }
        else{
            for( int j = 0; j < _nslices; j++ )
                _reconstructed[j][i] = dst[j];        
        }

    }   

    size_t n = nCoeffs / ratio;
//    if( n != inCount )
//        cerr << "WARNING: should use " << n << " coeffs, but actually used : "
//             << inCount << endl;
    
    return 0;
}

float*
SamSliceGroup::GetReconstructedPtr( int i )
{
    if( i < _reconstructed.size() )
        return _reconstructed[i];
    else{
        cerr << "SamSliceGroup::GetDecompressedCoeff(), idx too large! " << endl;
        return NULL;
    }
}


float
SamSliceGroup::FindCoeffThreshold( int ratio )
{
    size_t nCoeffs = _nslices * _rawlen;
    size_t n = nCoeffs / ratio - 1; // the nth largest, indexing from 1.


    vector<float>* allCoeffs = new vector<float>( );
    for( size_t i = 0; i < _coeffs.size(); i++ )
        for( int j = 0; j < _nslices; j++ ) 
            if ( _coeffs[i][j] > 0 )
                allCoeffs -> push_back( -1.0 * _coeffs[i][j] ); 
            else
                allCoeffs -> push_back( _coeffs[i][j] ); 
    
    if( allCoeffs->size() != nCoeffs )
        cerr << "WARNING: coeffs number don't match in SamSliceGroup::FindCoeffThreshold()!" << endl;

    std::nth_element( allCoeffs->begin(), allCoeffs->begin()+n, allCoeffs->end() );
    float nth = -1.0 * (allCoeffs -> at(n));

    delete allCoeffs;

    return nth;
}


int 
SamSliceGroup::SetupCompressor1( string wavename )
{
    vector< size_t > dims;
    dims.push_back( _nslices );
    
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

void
SamSliceGroup::FreeReconstructed( int i )
{
    assert( i < _reconstructed.size() );
    if( _reconstructed[i] ) {
            delete[] _reconstructed[i];
            _reconstructed[i] = NULL;
    }
}

void
SamSliceGroup::Print1DRaw()
{
    for( int i = 0; i < _nslices; i++ )
        cerr << "\t" << _raw[i][0] << endl;
}


void
SamSliceGroup::Print1DCoeffs()
{
    for( int i = 0; i < _nslices; i++ )
        cerr << "\t" << _coeffs[0][i] << endl;
}


