
#include "vapor/SamTI.h"

using namespace VAPoR;

SamTI::SamTI( int nfiles, size_t coeffArraySize, 
              const vector< float* > &arrays )
{
    _nfiles = nfiles;

    _coeffArraySize = coeffArraySize;

    _c1 = NULL;

    if( arrays.size() != _nfiles ){
        cerr << "ERROR: SamTI::SetupCoeffArrays() expecting " << _nfiles
             << " coefficient arrays, but got " << arrays.size() << endl;
        exit(1);
    }
    
    _coeffArrays = arrays;

    _nTCoeffArray = 0;
    _TCoeffSize = 0;
    _decompressedCoeffSize = 0;
}

SamTI::~SamTI()
{
    for( int i = 0; i < _TCoeffArrays.size(); i++ )
        if( _TCoeffArrays[i] )      delete[] _TCoeffArrays[i];
    for( int i = 0; i < _TSigmaps.size(); i++ )
        if( _TSigmaps[i] )          delete _TSigmaps[i];
    for( int i = 0; i < _decompressedCoeffs.size(); i++ )
        if( _decompressedCoeffs[i] )    delete _decompressedCoeffs[i];

    if( _c1 )                delete _c1; 
}


int 
SamTI::SetupTCompressor( string wavename )
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
SamTI::Compress1( size_t nTCoeffArray, size_t TCoeffSize )
{

    _nTCoeffArray = nTCoeffArray;
    _TCoeffSize   = TCoeffSize;

    for( int i = 0; i < _TCoeffArrays.size(); i++ )
        if( _TCoeffArrays[i] )      delete[] _TCoeffArrays[i];
    for( int i = 0; i < _TSigmaps.size(); i++ )
        if( _TSigmaps[i] )          delete _TSigmaps[i];

    int rc;
    float src[ _nfiles ];
    for( size_t i = 0; i < _nTCoeffArray; i++ ){

        for( int j = 0; j < _nfiles; j++) {
            src[j] = _coeffArrays[j][i];
        }
        
        float* dst = new float[ TCoeffSize ];
        SignificanceMap* sigmap = new SignificanceMap();
//cerr << "start compression" << endl;
        rc = _c1 -> Compress( src, dst, TCoeffSize, sigmap );
        if( rc != 0 ) {
            cerr << "compression on time intervals has error: " << rc << endl;
            return rc;
        }
//else    cerr << "compression succeed!" << endl;
        
        _TCoeffArrays.push_back( dst );
        _TSigmaps.push_back( sigmap );
    }

    return 0;
}

int
SamTI::Decompress1()
{
    _decompressedCoeffSize = _coeffArraySize;
    for( int i = 0; i < _nfiles; i++ )
        _decompressedCoeffs.push_back( new float[ _decompressedCoeffSize ] );

    float* src;
    float dst[ _nfiles ];
    int rc;
    
    for( size_t i = 0; i < _nTCoeffArray; i++ ) {
        src = _TCoeffArrays[i];
        rc = _c1 -> Decompress( src, dst, _TSigmaps[i] );
        if( rc != 0 ) {
            cerr << "Decompress() error with code: " << rc << endl;
            return 1;
        }
        else{
            for( int j = 0; j < _nfiles; j++ )
                _decompressedCoeffs[j][i] = dst[j];        
        }
    }   

    for( int i = 0; i < _nfiles; i++ )
        for( size_t j = _nTCoeffArray; j < _decompressedCoeffSize; j++ )
            _decompressedCoeffs[i][j] = 0.0;
    
    return 0;
}

float*
SamTI::GetDecompressedCoeff( size_t i )
{
    if( i < _decompressedCoeffs.size() )
        return _decompressedCoeffs[i];
    else{
        cerr << "SamTI::GetDecompressedCoeff(), idx too large! " << endl;
        return NULL;
    }
}











