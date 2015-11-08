#include "vapor/SamTS.h"

using namespace VAPoR;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

SamTS::SamTS( string filename ) {

    _filename = filename;
    _raw = NULL;
    _coeffs = NULL;
    _reconstructed = NULL;
    _coeffsInPosition = NULL;

    _reconstructed2 = NULL;

    FILE* f = fopen( filename.c_str(), "rb" );
    if( f != NULL ) {

        fseek( f, 0, SEEK_END );
        long size = ftell( f );
        _rawsize = size / 4;
        _raw = new float[ _rawsize ];

        fseek( f, 0, SEEK_SET );
        size_t rsize = fread( _raw, sizeof(float), _rawsize, f);
        fclose( f );
        if( rsize != _rawsize ) {
            cerr << "read size error: " << rsize << endl;
            exit(1);
        }
    
        _coeffs = new float[ _rawsize ];
    }
    else{
        cerr << "file open error: " << filename << endl;
        exit(1);
    }
}

SamTS::~SamTS(){
    if( _raw )          delete[] _raw;
    if( _coeffs )       delete[] _coeffs;
    if( _reconstructed) delete[] _reconstructed;
    if( _coeffsInPosition)  delete[] _coeffsInPosition;
    if( _reconstructed2 ) delete[] _reconstructed2;
}

void
SamTS::Print1stRaw(){
    if( _raw == NULL )
        cerr << "\t_raw == NULL! " << endl;
    else{
        cerr << "\t" << _raw[0] << endl; 
    }
}

void 
SamTS::SetupSigmap( size_t cvectorsize, int nsigmaps )
{
    _cvectorsize = cvectorsize;
    _coeffs = new float[ cvectorsize ];
    _originSigmaps.resize( nsigmaps );
}

void
SamTS::SetupReconstructedArray()
{
    _reconstructedSize = 0;
    for( int i = 0; i < _originSigmaps.size(); i++ )    
        _reconstructedSize += _originSigmaps[i].GetNumSignificant();
    if( _reconstructedSize != _rawsize )
        cerr << "Reconstructed size is " << _reconstructedSize
             << ", while raw size is " << _rawsize << endl;
    _reconstructed = new float[ _reconstructedSize ];
}

void
SamTS::PutCoeffsInPosition( )
{
    _coeffsInPosition = new float[ _cvectorsize ];
    for( int i = 0; i < _cvectorsize; i++ ) 
        _coeffsInPosition[i] = 0;

    size_t totalidx = 0;
    for( int i = 0; i < _originSigmaps.size(); i++ ) {
        for( int j = 0; j < _originSigmaps[i].GetNumSignificant(); j++ ) {
            size_t idx;
            _originSigmaps[i].GetCoordinates( j, &idx );
            _coeffsInPosition[ idx ] = _coeffs[ totalidx++ ];
        }
    }
}

float*
SamTS::GetReconstructed2Ptr( )
{
    if( _reconstructed2 == NULL )
        _reconstructed2 = new float[_rawsize];
    
    return _reconstructed2;
}
