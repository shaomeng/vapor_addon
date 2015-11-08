
#include "vapor/Sam1D.h"

using namespace VAPoR;

Sam1D::Sam1D( int nfiles, size_t rawlen, 
              const vector< float* > &raw )
{
    _nfiles = nfiles;
    _rawlen = rawlen;
    _raw = raw;

    if( raw.size() != _nfiles ){
        cerr << "ERROR: Sam1D::Sam1D() expecting " << _nfiles
             << " coefficient arrays, but got " << raw.size() << endl;
        exit(1);
    }

    _c1 = NULL;
    _coeffs.resize( rawlen, NULL);
    _sigmaps.resize( rawlen, NULL );
    _reconstructed.resize( nfiles, NULL );
}

Sam1D::~Sam1D()
{
    for( int i = 0; i < _coeffs.size(); i++ )
        if( _coeffs[i] )      delete[] _coeffs[i];
    for( int i = 0; i < _sigmaps.size(); i++ )
        if( _sigmaps[i] )          delete _sigmaps[i];
    for( int i = 0; i < _reconstructed.size(); i++ )
        if( _reconstructed[i] )    delete[] _reconstructed[i];

    if( _c1 )                delete _c1; 
}


int 
Sam1D::SetupCompressor1( string wavename )
{
    vector< size_t > dims;
    dims.push_back( _nfiles );
    
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
Sam1D::Compress1( )
{

    for( int i = 0; i < _coeffs.size(); i++ )
        if( _coeffs[i] )      delete[] _coeffs[i];
    for( int i = 0; i < _sigmaps.size(); i++ )
        if( _sigmaps[i] )          delete _sigmaps[i];

    int rc;
    float src[ _nfiles ];
    for( size_t i = 0; i < _rawlen; i++ ){

        for( int j = 0; j < _nfiles; j++) {
            src[j] = _raw[j][i];
        }
        
        float* dst = new float[ _nfiles ];
        SignificanceMap* sigmap = new SignificanceMap();
        rc = _c1 -> Compress( src, dst, _nfiles, sigmap );
        if( rc != 0 ) {
            cerr << "compression on time intervals has error: " << rc << endl;
            return rc;
        }
      //else    cerr << "compression succeed!" << endl;
        
        _coeffs[i] = dst;
        _sigmaps[i] = sigmap;
    }

    return 0;
}

int
Sam1D::Decompress1( int ratio )
{
    float nth = FindCoeffThreshold( ratio ); // use coeffs larger than nth.

    for( int i = 0; i < _nfiles; i++ )
        _reconstructed[i] = new float[ _rawlen ];

    float src[ _nfiles], dst[ _nfiles ];
    int rc;
    size_t inCount = 0;
    
    for( size_t i = 0; i < _rawlen; i++ ) {
        
        for( int j = 0; j < _nfiles; j++ ) {
            float c = _coeffs[i][j];
            if( c > nth || c < -1.0 * nth ) {
                src[j] = c;
                inCount++;
            }
            else 
                src[j] = 0.0;
        }

        rc = _c1 -> Decompress( src, dst, _sigmaps[i] );

        if( rc != 0 ) {
            cerr << "Decompress() error with code: " << rc << endl;
            return 1;
        }
        else{
            for( int j = 0; j < _nfiles; j++ )
                _reconstructed[j][i] = dst[j];        
        }

    }   

    size_t nCoeffs = _nfiles * _rawlen;
    size_t n = nCoeffs / ratio;
//    if( n != inCount )
        cerr << "WARNING: should use " << n << " coeffs, but actually used : "
             << inCount << endl;
    
    return 0;
}

float*
Sam1D::GetDecompressedCoeff( size_t i )
{
    if( i < _reconstructed.size() )
        return _reconstructed[i];
    else{
        cerr << "Sam1D::GetDecompressedCoeff(), idx too large! " << endl;
        return NULL;
    }
}


float
Sam1D::FindCoeffThreshold( int ratio )
{
    if( ratio == 1 )    // use all coeffs
        return 0.0;

    size_t nCoeffs = _nfiles * _rawlen;
    if( nCoeffs % ratio != 0 )
        cerr << "WARNING: requesting compression ratio is not powers of two: " 
             << ratio << endl;
    size_t n = nCoeffs / ratio;


    vector<float>* allCoeffs = new vector<float>( );
    for( size_t i = 0; i < _coeffs.size(); i++ )
        for( int j = 0; j < _nfiles; j++ ) 
            if ( _coeffs[i][j] > 0 )
                allCoeffs -> push_back( -1.0 * _coeffs[i][j] ); 
            else
                allCoeffs -> push_back( _coeffs[i][j] ); 
    
    if( allCoeffs->size() != nCoeffs )
        cerr << "WARNING: coeffs number don't match in Sam1D::FindCoeffThreshold()!" << endl;


    std::nth_element( allCoeffs->begin(), allCoeffs->begin()+n, allCoeffs->end() );
    float mth = -1.0 * (allCoeffs -> at(n-1));
    float nth = -1.0 * (allCoeffs -> at(n));
    float oth = -1.0 * (allCoeffs -> at(n+1));

//std::sort( allCoeffs->begin(), allCoeffs->end() );
//cerr << "nth=" << allCoeffs->at(n);

cerr << "n-1= " << mth << ", n= " << nth << ", n+1= " << oth << endl;

    delete allCoeffs;

    return nth;
}


void
Sam1D::PrintRawThresholds( int ratio )
{
    if( _rawlen % ratio != 0 )
        cerr << "WARNING: requesting compression ratio is not powers of two: "
             << ratio << endl;
    size_t n = _rawlen / ratio;

    vector< float >* vec = new vector< float >( _rawlen, 0.0 );
    cerr << "The threshold coefficient magnitude at cratio " << ratio << ":1 are:" << endl;
    if( ratio == 1 )    // use all coeffs
        cerr << "\t" << 0.0 << endl;
    else{
        for( int i = 0; i < _nfiles; i++ ) {
            float* pt = _raw[i];
            for( size_t j = 0; j < _rawlen; j++ ) {
                if( pt[j] > 0 )
                    vec->at(j) = -1.0 * pt[j];
                else
                    vec->at(j) = pt[j];
            }
            
            std::nth_element( vec->begin(), vec->begin()+n, vec->end() );
            float nth = -1.0 * (vec->at(n));
            cerr << "\t" << nth << endl;
        }
    }

    delete vec;
}







