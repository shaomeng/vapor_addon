// This program calculates error introduced by recovering the array using partial coefficients, 
// e.g. the first half portion of the coefficients without prioritizing them.


#include "vapor/SamToolbox.h"
#include "vapor/SamTS.h"
#include "vapor/SamTI.h"

using namespace VAPoR;


int main(int argc, char* argv[] ) {

    int lod = 0;    
    string path = "/Users/samuel/Backyard/256cubes";
    if( argc == 2 )
        lod = std::stoi( argv[1] );
    else if( argc == 3 ) {
        lod = std::stoi( argv[1] );
        path = argv[2];
    }
    else {
        cerr << "Wrong number of parameters!" << endl;
        exit(1);
    }
    if( path.back() == '/' )
        path.pop_back();

    SamToolbox toolbox;

    size_t onedim = 256;
    size_t totallen = onedim * onedim * onedim;


// Read in nfiles time slices
    int nfiles = 16;
    string* filenames = new string[ nfiles ];
    for( int i = 0; i < nfiles; i++ ) {
        filenames[i] = path + "/e" + std::to_string(i) + ".float";
    }
    
    SamTS** slices = new SamTS* [ nfiles ];
    for( int i = 0; i < nfiles; i++ ) {
        slices[i] = new SamTS( filenames[i] );
//        cerr << "finish reading slice " << i << endl;
    }

    int rc;

    vector< float* > raw;
    for( int i = 0; i < nfiles; i++ ) {
        raw.push_back( slices[i] -> GetRawPtr() );
    }


// Perform DWT on time intervals 
    int TICoeffSize = nfiles / 2 + lod;    // how many coefficients to use on each time interval
    size_t nTICoeffArray = totallen;  

    SamTI* interval = NULL;
    interval = new SamTI( nfiles, totallen, raw );
    rc = interval -> SetupTCompressor( "bior4.4" );
    if( rc != 0 )
        cerr << "time interval compressor setup error: " << rc << endl;

    cerr << "compressing " << nTICoeffArray << " time intervals..." << endl;
    rc = interval -> Compress1( nTICoeffArray,  TICoeffSize );
    if( rc != 0 )
        cerr << "compression on time intervals error: " << rc << endl;
    else
        cerr << "finish compression on time intervals" << endl;
    

// Perform IDWT on time intervals to recover the coefficients
    rc = interval -> Decompress1();
    if( rc != 0 )
        cerr << "decompression on intervals error: " << rc << endl;
    else
        cerr << "finish decompression coefficients (time interval)" << endl;


    
// Collect stats for Decompressed arrays
    for( int i = 0; i < nfiles; i++ ) {
        float* raw = slices[i] -> GetRawPtr();
        float* decompressed = interval -> GetDecompressedCoeff(i);
        size_t size = totallen;

        toolbox.CompareArrays( raw, decompressed, size );        
    }


    




// Free memory
    if( interval )                          delete interval;
    for( int i = 0; i < nfiles; i++ )
        if( slices[i] )                     delete slices[i];
    if( slices )                            delete[] slices;
    if( filenames )                         delete[] filenames;
}
