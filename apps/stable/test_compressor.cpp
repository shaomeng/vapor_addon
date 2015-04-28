#include "vapor/SamToolbox.h"
#include "vapor/SamSlice.h"


using namespace VAPoR;


int main(int argc, char* argv[] ) {

    int NX = 18;
    int NY = 18;
    if( argc == 2 )
        NX = std::stoi( argv[1] );

    vector< size_t > dims;
    dims.push_back( NX );
    dims.push_back( NY );

    float src[ NX*NY ], dst[ NX*NY ];

    cerr << "src array: " << endl;
    for( int i = 0; i < NX*NY; i++ ) {
        src[i] = sin( i*0.1 ) * 100;
//        cerr << "\t" << src[i] << endl;
    }

    
    Compressor c( dims, "bior4.4", "symw" ); 
//    Compressor c( dims, "bior3.3", "symh" ); 
//    Compressor c( dims, "bior1.1", "symh" ); 

    SignificanceMap map;
    c.Compress( src, dst, NX*NY, &map );

/*
    cerr << "dst array: " << endl;
    for( int i = 0; i < NX; i++ )
        cerr << "\t" << dst[i] << endl;
*/
    
    cerr << "Num of Significant: " << map.GetNumSignificant() << endl;
}
