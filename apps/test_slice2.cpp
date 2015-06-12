#include "vapor/SamSlice2.h"
#include "vapor/SamToolbox.h"

#define NX 504
#define NY 504
#define NZ 504
#define NSLICES 9

using namespace VAPoR;

int main( int argc, char* argv[] )
{
    if( argc < 3 ) {
        cerr << "./a.out filepath ratio" << endl;
        exit(1);
    }

    string filename = argv[1];
    int ratio = stoi( argv[2] );
    SamToolbox toolbox;

    SamSlice2* slice = new SamSlice2( filename, "bior4.4", NX, NY, NZ );
    slice -> Decompose();
    slice -> Reconstruct( ratio );
    float* rawPtr = slice -> GetRawPtr();
    float* reconstructedPtr = slice -> GetReconstructedPtr();
    toolbox.CompareArrays( rawPtr, reconstructedPtr, NX*NY*NZ, true );

    delete slice;
}
