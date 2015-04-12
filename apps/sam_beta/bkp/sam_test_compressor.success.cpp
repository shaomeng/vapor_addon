
#include "vapor/SamToolbox.h"

using namespace VAPoR;


int main(int argc, char* argv[] ) {

    if( argc != 2 ) {
        cerr << "Wrong number of parameters!" << endl;
        cerr << "Please specify the file to read!" << endl;
        exit(1);
    }

    SamToolbox sam;
    int rc;
    size_t onedim = 256;
    size_t totallen = onedim * onedim * onedim;
    float* src = new float[ totallen ];
//    for( size_t i = 0; i < totallen; i++ )
//        src[i] = float( rand() % 10000 ) / 100.0;
    string input = argv[1];
    sam.ReadBlock( input, totallen, src );
    vector<size_t> dims;    
    dims.push_back( onedim );
    dims.push_back( onedim );
    dims.push_back( onedim );
    vector <size_t> cratios;
    for( int i = 8; i >= 1; i /= 2 )
        cratios.push_back(i);
    Compressor* c1 = NULL;
    vector< size_t > ncoeffs;
    size_t cvectorsize;
    vector< SignificanceMap > sigmaps;
    sam.SetupCompressor( c1, dims, cratios, ncoeffs, cvectorsize, sigmaps );
    float* coeffs = new float[ cvectorsize ];

    cerr << "Number of collections of coeffs: " << ncoeffs.size() << endl;
    for ( int i = 0; i < ncoeffs.size(); i++ )
        cerr << "\t" << ncoeffs[i] << endl;
	cerr << "Number of Total coeffs: " << cvectorsize << endl;

    rc = c1 -> Decompose( src, coeffs, ncoeffs, sigmaps );
    if( rc != 0 ) {
        cerr << "Decompose error with error code: " << rc << endl;
    }

    float* rarr = new float[ totallen ];

// Test reconstructing using Lod
    vector< SignificanceMap > sigmaps2 = sigmaps;
    sigmaps2.pop_back();    // Use the next coarse level

    size_t cvectorsize2 = 0;
    for( int i = 0; i < sigmaps2.size(); i++ )
        cvectorsize2 += sigmaps2[i].GetNumSignificant();

// Test what is stored in sigmaps
/*
    SignificanceMap tmpmap = sigmaps[0];
    size_t maplen = tmpmap.GetNumSignificant();
    for( int i = 0; i < 10; i++ ) {
        size_t idx;
        tmpmap.GetCoordinates( maplen-10+i, &idx );
        cerr << idx << endl;
    }
*/

	cerr << "Number of coeffs for reconstruction: " << cvectorsize2 << endl;
    float* coeffs2 = new float[cvectorsize2];
    memcpy( coeffs2, coeffs, sizeof(float) * cvectorsize2 );

//    int rlevel = c1 -> GetNumLevels();
    rc = c1 -> Reconstruct( coeffs2, rarr, sigmaps2, -1);
    if( rc != 0 ){
        cerr << "Reconstruct error with error code: " << rc << endl;
    }
    
    sam.CompareArrays( src, rarr, totallen);    

    delete[] src;
    delete[] coeffs;
    delete[] coeffs2;
    delete[] rarr;
    delete c1;

}
