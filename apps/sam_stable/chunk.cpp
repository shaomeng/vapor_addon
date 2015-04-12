// This program evaluates a chunk of data (not necessarily to be cubic )

#include "vapor/SamToolbox.h"
#include "vapor/SamSlice.h"

#define NX 256
#define NY 256
#define NZ 256

using namespace VAPoR;


int main(int argc, char* argv[] ) {

    int ratio = 0;    
    string file = "/Users/samuel/Backyard/256cubes/e0chunks/18layer.float";
    if( argc == 2 )
         ratio = std::stoi( argv[1] );
    else if( argc == 3 ) {
        ratio = std::stoi( argv[1] );
        file = argv[2];
    }
    else {
        cerr << "Wrong number of parameters!" << endl;
        exit(1);
    }

    SamToolbox toolbox;
    SamSlice   slice( file );
    int rc;

    vector< size_t > dims;
    dims.push_back( NX );
    dims.push_back( NY );
    dims.push_back( NZ );
    rc = slice.SetupCompressor1( dims, "bior4.4" );
    if( rc != 0 )   cerr << "Setup compressor error: " << rc << endl;
    rc = slice.Compress1();
    if( rc != 0 )   cerr << "Compress error: " << rc << endl;
    rc = slice.Decompress1( ratio );
    if( rc != 0 )   cerr << "Decompress error: " << rc << endl;


    float* raw = slice.GetRawPtr();
    float* reconstructed = slice.GetReconstructedPtr();
    size_t size = slice.GetRawSize();
    toolbox.CompareArrays( raw, reconstructed, size, true );

// Test PutCoeffsInPosition() function
/*
    slice.PutCoeffsInPosition();

    rc = slice.Decompress1( ratio );
    if( rc != 0 )   cerr << "Decompress error: " << rc << endl;
    toolbox.CompareArrays( raw, reconstructed, size, true );
*/

}
