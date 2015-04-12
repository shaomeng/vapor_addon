// This program reads, stores, transforms, and compares 
// each slice separaely.

#include "vapor/SamToolbox.h"
#include "vapor/SamTS.h"

using namespace VAPoR;


int main(int argc, char* argv[] ) {

    string path;
    if( argc == 2 )
        path = argv[1];
    else if( argc == 1 )
        path = "/Users/samuel/Backyard/256cubes";
    else {
        cerr << "Wrong number of parameters!" << endl;
        cerr << "Please specify the path to the folder!" << endl;
        exit(1);
    }
    if( path.back() == '/' )
        path.pop_back();

    SamToolbox toolbox;

// Setup Compressor
    size_t onedim = 256;
    vector<size_t> dims;
    dims.push_back( onedim );
    dims.push_back( onedim );
    dims.push_back( onedim );
    Compressor* c1 = new Compressor( dims, "bior4.4", "symw" );

// Setup compression ratio, and number of coefficients
    vector <size_t> cratios;
    for( int i = 8; i >= 1; i /= 2 )
        cratios.push_back(i);
    vector< size_t > ncoeffs;
    size_t cvectorsize;
    toolbox.SetupNCoeffs( c1, cratios, ncoeffs, cvectorsize );

// Read in 8 time slices
    int nfiles = 8;
    string* filenames = new string[ nfiles ];
    for( int i = 0; i < nfiles; i++ ) {
        filenames[i] = path + "/e" + std::to_string(i) + ".float";
    }
    
    SamTS** slices = new SamTS* [ nfiles ];
    for( int i = 0; i < nfiles; i++ ) {
        slices[i] = new SamTS( filenames[i] );
        cerr << "finish reading slice " << i << endl;
    }

    int rc;

// Setup Sigmaps, and then DWT on each time slices
    for( int i = 0; i < nfiles; i++ ){
        slices[i] -> SetupSigmap( cvectorsize, cratios.size() );

        float* src = slices[i] -> GetRawPtr();
        float* dst = slices[i] -> GetCoeffsPtr();
        vector< SignificanceMap >* sigmaps = slices[i] -> GetSigMapsPtr();
        
        rc = c1 -> Decompose( src, dst, ncoeffs, *sigmaps);
        if( rc != 0 ) {
            cerr << "Decompose error with code: " << rc << endl;
        }
        else
            cerr << "finish decomposing slice: " << i << endl;
    }

// Reconstruct here
    for( int i = 0; i < nfiles; i++) {
        slices[i] -> SetupReconstructedArray();
        float* coeffs = slices[i] -> GetCoeffsPtr();
        float* reconstructed = slices[i] -> GetReconstructedPtr();
        vector< SignificanceMap >* sigmaps = slices[i] -> GetSigMapsPtr();

        // level of detail, used in prioritized coefficients
        // This should be between 1 and NumberOfBins.
        int lod = 1;
        vector< SignificanceMap > sigmaps2;
        SignificanceMap* copySigmaps[lod];
        for( int i = 0; i < lod; i++ ) {
            copySigmaps[i] = new SignificanceMap( (*sigmaps)[i] );
            sigmaps2.push_back( *copySigmaps[i] );
        }

        rc = c1 -> Reconstruct( coeffs, reconstructed, sigmaps2, -1 );
        if( rc != 0 )
            cerr << "Reconstruct error with code: " << rc << endl;
        else
            cerr << "finish reconstructing slice: " << i << endl;
        
        for( int i = 0; i < lod; i++ )
            if( copySigmaps[i] )        delete copySigmaps[i];
    }

// Collect stats for reconstructed arrays
    for( int i = 0; i < nfiles; i++ ) {
        float* raw = slices[i] -> GetRawPtr();
        float* reconstructed = slices[i] -> GetReconstructedPtr();
        size_t size = slices[i] -> GetReconstructedSize();

        toolbox.CompareArrays( raw, reconstructed, size );        
    }
    

/*
for( int i = 0; i < nfiles; i++ ){
    cerr << "printing file " << i << endl;
    for( int j = 0; j < slices[i]-> sigmaps.size(); j++ )
        cerr << "\t" << slices[i]-> sigmaps[j].GetNumSignificant() << endl;
}
*/

// Free memory
    if( c1 )                                delete c1;
    for( int i = 0; i < nfiles; i++ )
        if( slices[i] )                     delete slices[i];
    if( slices )                            delete[] slices;
    if( filenames )                         delete[] filenames;
}
