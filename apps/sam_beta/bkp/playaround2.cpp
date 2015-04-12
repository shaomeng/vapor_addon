//* This program performs a whole set of tests:
//* 1) read in data
//* 2) perform DWT on spatial dimensions
//* 3) IDWT, and get the baseline results
//*
//* 4) put coefficients back to their original locations
//* 5) DWT on time intervals (short arrays)
//* 6) restore the array from time interval coeffs
//* 7) evaluate the time interval results
//* 
//* It essentially compares if the added time-dimension DWT helps...


#include "vapor/SamToolbox.h"
#include "vapor/SamTS.h"
#include "vapor/Sam1D.h"

using namespace VAPoR;


int main(int argc, char* argv[] ) {

    int lod = 0;    
    string path = "/Users/samuel/Backyard/256cubes/e0slices";
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
    if( path.back() != '/' )
        path.push_back('/');

    SamToolbox toolbox;

// Setup Compressor
    size_t onedim = 256;
    size_t totallen = onedim * onedim;
    vector<size_t> dims;
    dims.push_back( onedim );
    dims.push_back( onedim );
    Compressor* c1 = NULL;
    c1 = new Compressor( dims, "bior4.4", "symw" );

// Setup compression ratio, and number of coefficients
    vector <size_t> cratios;
    for( int i = 32; i >= 1; i /= 2 )
        cratios.push_back(i);   // lod:   0    1    2   3   4   5   
                                // ratio: 32:1 16:1 8:1 4:1 2:1 1:1
    vector< size_t > ncoeffs;
    size_t cvectorsize;
    toolbox.SetupNCoeffs( c1, cratios, ncoeffs, cvectorsize );
    // ( cvectorsize should equal to totallen, and should be the sum of nums in ncoeffs )

// Read in nfiles time slices
    int nfiles = 14;
    string* filenames = new string[ nfiles ];
    for( int i = 0; i < nfiles; i++ ) {
        filenames[i] = path + std::to_string(i);
    }
    
    SamTS** slices = new SamTS* [ nfiles ];
    for( int i = 0; i < nfiles; i++ ) {
        slices[i] = new SamTS( filenames[i] );
//        cerr << "finish reading slice " << i << endl;
    }

    int rc;

// Setup Sigmaps, and then DWT on each time slices
    for( int i = 0; i < nfiles; i++ ){
        slices[i] -> SetupSigmap( cvectorsize, cratios.size() );

        float* src = slices[i] -> GetRawPtr();
        float* dst = slices[i] -> GetCoeffsPtr();
        vector< SignificanceMap >* sigmaps = slices[i] -> GetOriginSigMapsPtr();
        
        rc = c1 -> Decompose( src, dst, ncoeffs, *sigmaps);
        if( rc != 0 ) {
            cerr << "Decompose error with code: " << rc << endl;
        }
//        else
//            cerr << "finish decomposing slice: " << i << endl;
    }


    int binToUse = lod+1;


// Reconstruct here
    for( int i = 0; i < nfiles; i++) {
        //slices[i] -> SetLOD( lod );
        slices[i] -> SetupReconstructedArray();
        float* coeffs = slices[i] -> GetCoeffsPtr();
        float* reconstructed = slices[i] -> GetReconstructedPtr();
        vector< SignificanceMap >* sigmaps = slices[i] -> GetOriginSigMapsPtr();

        // create new sigmaps used for partial reconstruction
        vector< SignificanceMap > sigmaps2;
        SignificanceMap* copySigmaps[ binToUse ];
        for( int i = 0; i <  binToUse ; i++ ) {
            copySigmaps[i] = new SignificanceMap( (*sigmaps)[i] );
            sigmaps2.push_back( *copySigmaps[i] );
        }

        rc = c1 -> Reconstruct( coeffs, reconstructed, sigmaps2, -1 );
        if( rc != 0 )
            cerr << "Reconstruct error with code: " << rc << endl;
//        else
//            cerr << "finish reconstructing slice: " << i << " in lod " << lod << endl;
        
        for( int i = 0; i <  binToUse ; i++ )
            if( copySigmaps[i] )        
                delete copySigmaps[i];
    }

// Collect stats for reconstructed arrays
    cerr << "Stats of reconstructed arrays at lod: " << lod << endl;
    int nprints;     // how many lines to print;
    nprints = 9 > nfiles? nfiles : 9;    // normally be nfiles;
    for( int i = 0; i < nprints; i++ ) {
        float* raw = slices[i] -> GetRawPtr();
        float* reconstructed = slices[i] -> GetReconstructedPtr();
        size_t size = totallen;

        toolbox.CompareArrays( raw, reconstructed, size );        
    }


// Put coefficients in right positions 
    for( int i = 0; i < nfiles; i++ ) {
        slices[i] -> PutCoeffsInPosition(  );
//        float* coeffsInPosition = slices[i] -> GetCoeffsInPosition();
//        toolbox.PrintArray( coeffsInPosition, 0, 20 );
    }
    vector< float* > coeffsInPosition;
    for( int i = 0; i < nfiles; i++ ) {
        coeffsInPosition.push_back( slices[i] -> GetCoeffsInPosition() );
    }


// Perform DWT on time intervals 
    int ratio;
    if( lod == 0 )          ratio = 32;
    else if( lod == 1 )     ratio = 16;
    else if( lod == 2 )     ratio = 8;
    else if( lod == 3 )     ratio = 4;
    else if( lod == 4 )     ratio = 2;
    else                    ratio = 1;


    Sam1D* sam1d = NULL;
    sam1d = new Sam1D( nfiles, totallen, coeffsInPosition );
    rc = sam1d -> SetupCompressor1( "bior4.4" );
    if( rc != 0 )
        cerr << "sam1d compressor setup error: " << rc << endl;
    rc = sam1d -> Compress1( );
    if( rc != 0 )
        cerr << "compression on time intervals error: " << rc << endl;
    else
        cerr << "finish compression on time intervals" << endl;
    

// Perform IDWT on time intervals to recover the coefficients
    rc = sam1d -> Decompress1( ratio );
    if( rc != 0 )
        cerr << "decompression on intervals error: " << rc << endl;
    else
        cerr << "finish decompression coefficients (time interval)" << endl;


// Perform IDWT on recovered coefficients to recover original array
    SignificanceMap* inOrderSigmap = NULL;
    inOrderSigmap = new SignificanceMap( dims );    
    for( size_t i = 0; i < totallen; i++ )
        inOrderSigmap -> Set (i );
    for( int i = 0; i < nfiles; i++ ) {
        float* src = sam1d -> GetDecompressedCoeff( i );
        float* dst = slices[i] -> GetReconstructed2Ptr();
        rc = c1 -> Decompress( src, dst, inOrderSigmap );
        if( rc != 0 )
            cerr << "error when decompressing file (time interval) " << i << endl;
//        else
//            cerr << "finish decompressing file (time interval) " << i << endl;
    }
    
// Collect stats for Decompressed arrays
    cerr << "Compression with separate last dimension:" << endl;
    for( int i = 0; i < nprints; i++ ) {
        float* raw = slices[i] -> GetRawPtr();
        float* decompressed = slices[i] -> GetReconstructed2Ptr();
        size_t size = totallen;

        toolbox.CompareArrays( raw, decompressed, size );        
    }


    

// Get histogram of how many positions have non-zero values
/*
    vector< size_t > histogram;
    float epsilon = 1e-6;
    toolbox.CalcHistogram( coeffsInPosition, totallen, epsilon, histogram );
    for( int i = 0; i < histogram.size(); i++ )
        cerr << "\t" << i << ":  " << histogram[i] << endl;
*/



// Free memory
    if( sam1d )                          delete sam1d;
    if( c1 )                                delete c1;
    for( int i = 0; i < nfiles; i++ )
        if( slices[i] )                     delete slices[i];
    if( slices )                            delete[] slices;
    if( filenames )                         delete[] filenames;
    if( inOrderSigmap )                     delete inOrderSigmap;
}
