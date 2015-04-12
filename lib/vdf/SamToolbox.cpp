
#include "vapor/SamToolbox.h"

using namespace VAPoR;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

int 
SamToolbox::SetupCompressor(
    Compressor* &_compressor, // output
    const vector<size_t> &bdims, // input
    const vector<size_t> &_cratios, // input
    vector< size_t > &_ncoeffs, // output: num wave coeff. at each compression level        
    size_t &_cvectorsize, // output: amount of space need for wavelet coefficients
    vector< SignificanceMap > &sigmaps )  // output: a signficance map for each compression level
{
    string wavename = "bior4.4";
    string boundary_mode = "symw";

    // compression ratios read in from the vdf file.
    _compressor = new Compressor( bdims, wavename, boundary_mode);

    if (Compressor::GetErrCode() != 0) {
        cerr << "Failed to allocate compressor" << endl;
        return(-1);
    }

    // Total number of wavelet coefficients in a forward transform
    size_t ntotal = _compressor->GetNumWaveCoeffs();

    size_t mincoeff = _compressor->GetMinCompression();
    assert(mincoeff >0);

    for (int i=0; i<_cratios.size(); i++) {
        size_t n = ntotal / _cratios[i];
        if (n < mincoeff) {
            cerr << "Invalid compression ratio for configuration : %" << _cratios[i];
            return(-1);
        }
    }

    _ncoeffs.clear();
    sigmaps.resize( _cratios.size() ); // a signficance map for each compression level
    size_t naccum = 0;

    for (int i = 0; i < _cratios.size(); i++) {
        if (_cratios[i] > ntotal) break;

        size_t n = (ntotal / _cratios[i]) - naccum;
        _ncoeffs.push_back(n);
        naccum += n;

        vector <size_t> sigdims;
        _compressor->GetSigMapShape(sigdims);
        sigmaps[i].Reshape(sigdims);
    }
    assert(naccum <= ntotal);

    _cvectorsize = naccum;

    return(0);   

}    

int 
SamToolbox::SetupNCoeffs(
    const Compressor* compressor,         // input
    const vector< size_t > &cratios,      // input:
    vector< size_t > &ncoeffs,            // output:
    size_t &cvectorsize )                 // output:
{
    size_t ntotal = compressor->GetNumWaveCoeffs();

    size_t mincoeff = compressor->GetMinCompression();
    assert(mincoeff >0);

    for (int i=0; i<cratios.size(); i++) {
        size_t n = ntotal / cratios[i];
        if (n < mincoeff) {
cerr << "mincoeff=" << mincoeff << ",  n=" << n << endl;
            cerr << "Invalid compression ratio for configuration : %" << cratios[i];
            return(-1);
        }
    }

    ncoeffs.clear();
    size_t naccum = 0;

    for (int i = 0; i < cratios.size(); i++) {
        if (cratios[i] > ntotal) break;

        size_t n = (ntotal / cratios[i]) - naccum;
        ncoeffs.push_back(n);
        naccum += n;
    }
    assert(naccum <= ntotal);

    cvectorsize = naccum;

    return(0);
}

float
SamToolbox::CompareArrays(
    const float* arr1,
    const float* arr2,
    size_t len,
    bool   print )
{
/*
    cerr << "comparing two arrays with length: " << len << endl;
    cerr << "The first ten elements are: " << endl;
    for( int i = 0; i < 10; i++ )
        cerr << "\t" << arr1[i] << "\t" << arr2[i] << endl; 
    cerr << "The last ten elements are: " << endl;
    for( int i = 10; i > 1; i-- )
        cerr << "\t" << arr1[len-i] << "\t" << arr2[len-i] << endl; 
*/

    float sum = 0.0;
    float c = 0.0;
    float max = 0.0;
    float v1, v2;
    float tmp;
    for( size_t i = 0; i < len; i++) {
        tmp = arr1[i] - arr2[i];
        if (tmp < 0)    tmp *= -1.0;
        float y = tmp * tmp - c;
        float t = sum + y;
        c = (t - sum) - y;
        sum = t;
        if (tmp > max){
            max = tmp;
            v1 = arr1[i];
            v2 = arr2[i];
        }
    }
    sum /= len * 1.0;
    sum = sqrt( sum );
    if( print )
        cerr << "\tRMS: " << sum << ", max error: " << max << endl;

    return sum;
}

int
SamToolbox::ReadBlock( 
    std::string filename, 
    size_t len, 
    float* ptr )
{
    FILE* f = fopen( filename.c_str(), "rb" );    
    if( f != NULL ) {
        fseek( f, 0, SEEK_SET );
        size_t rsize = fread( ptr, sizeof(float), len, f );
        fclose( f );
        if( rsize != len ) {
            cerr << "read size error: " << rsize << endl;
            return -1;
        }
        return 0;
    }
    
    cerr << "file open error: " << filename << endl;
    return -1;
}

void 
SamToolbox::PrintArray( const float* arr, size_t start, size_t end )
{
    cerr << "Printing elements from array: " << endl;
    for( size_t i = start; i < end; i++)
       cerr << "\t" << arr[i] << endl; 

}

void 
SamToolbox::CalcHistogram1( const vector< float* > &arrays, 
                           size_t arraylen, 
                           float epsilon,
                           vector< size_t > &histogram )
{
    int nfiles = arrays.size();
    histogram.clear();
    histogram.resize( nfiles + 1, 0 ); 

    for( size_t i = 0; i < arraylen; i++ ){
        int count = 0;
        for( int j = 0; j < nfiles; j++ ){
            // If the coeff at position i is non-zero
            if( arrays[j][i] > epsilon || arrays[j][i] < -1.0*epsilon )
                count++;
        }
        histogram[ count ]++;
    }

}


float 
SamToolbox::FindMax( const float* arr, size_t len ) {
    float max = 0;
    for( size_t i = 0; i < len; i++ )
        if( arr[i] > max )
            max = arr[i];
    return max;
}

float
SamToolbox::Findnth( const float* arr, size_t arrlen, size_t n)
{
    if( n == arrlen )
        return 0.0;

    vector< float >* vec = new vector< float >( arrlen, 0.0 );
    for( size_t i = 0; i < arrlen; i++ )
        if( arr[i] > 0 )
            vec -> at(i) = -1.0 * arr[i];
        else
            vec -> at(i) = arr[i];

    std::nth_element( vec->begin(), vec->begin()+n, vec->end() );
    float nth = -1.0 * (vec->at(n));
//    cerr << "\t" << nth << endl;

    delete vec;
    
    return nth;
}

float
SamToolbox::CalcRMS( const vector< float > &arr )
{
    float sum = 0.0;
    float c = 0.0;
    for( size_t i = 0; i < arr.size(); i++ ) {
        float y = arr[i] * arr[i] - c;
        float t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= 1.0 * arr.size();
    sum = sqrt( sum );    
    return sum;
}
