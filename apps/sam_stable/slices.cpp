#include <string>

#include "vapor/SamToolbox.h"
#include "vapor/SamSlice.h"
#include "vapor/SamSliceGroup.h"


#define NX 256
#define NY 256
#define NZ 256
#define NSLICE 18

using namespace VAPoR;

int main(int argc, char* argv[] ) {

    int ratio = 0;    
    string path = "/Users/samuel/Backyard/256cubes/";
    if( argc == 2 )
        ratio = std::stoi( argv[1] );
    else {
        cerr << "Wrong number of parameters!" << endl;
        exit(1);
    }
    string* filenames = new string[ NSLICE ];
    for( long long i = 0; i < NSLICE; i++ ) 
        filenames[i] = path + "e" + std::to_string(i) + ".float";

    SamToolbox toolbox;
    int rc;

    vector< size_t > dims;
    dims.push_back( NX );
    dims.push_back( NY );
    if( NZ != 1 )   dims.push_back( NZ );

    vector< SamSlice > slices( NSLICE );
    for( int i = 0; i < NSLICE; i++ ) {
        slices[i].ReadFile( filenames[i] );
        rc = slices[i].SetupCompressor1( dims, "bior4.4" );
        if( rc != 0 ) cerr << "setup compressor error: " << rc << endl;
        rc = slices[i].Compress1();
        if( rc != 0 )   cerr << "Compress error: " << rc << endl;
        else            cerr << "3D compression success on slice " << i << endl;

        rc = slices[i].Decompress1( ratio );
        if( rc != 0 )   cerr << "Decompress error: " << rc << endl;
        float* raw = slices[i].GetRawPtr();
        float* reconstructed = slices[i].GetReconstructedPtr();
        toolbox.CompareArrays( raw, reconstructed, NX*NY*NZ, true );

    }


    vector< float* > arrs;
    for( int i = 0; i < NSLICE; i++ ) {
        slices[i].PutCoeffsInPosition();
        arrs.push_back( slices[i].GetCoeffsPtr() );
    }

    SamSliceGroup group;
    group.Setup( NSLICE, NX*NY*NZ, arrs );

// Sam
//cerr << "Sam: Printing first 10 coeffs after repositioning: " << endl;
//for( int i = 0; i < 5; i++ ) {
//    slices[i].Print10Coeffs();
//    cerr << endl;
//}

// Sam
//cerr << "Sam: printing the first 18 elements before 1D DWT: " << endl;
//group.Print1DRaw();

    rc = group.SetupCompressor1( "bior4.4" );
    if( rc != 0 )   cerr << "SamSliceGroup::SetupCompressor1() error: " << rc << endl;
    rc = group.Compress1();
    if( rc != 0 )   cerr << "SamSliceGroup::Compress1() error: " << rc << endl;

// Sam
//cerr << "Sam: printing the first 18 elements after DWT: " << endl;
//group.Print1DCoeffs();

    rc = group.Decompress1( ratio );
    if( rc != 0 )   cerr << "SamSliceGroup::Decompress1() error: " << rc << endl;

    vector< float > rms;
    for( int i = 0; i < NSLICE; i++ ) {
        slices[i].ReplaceCoeffs( group.GetReconstructedPtr(i) );
        group.FreeReconstructed( i );
        rc = slices[i].Decompress1( 1 );
        if( rc != 0 )   cerr << "Decompress1() error at slice: " << i << endl;
        float* raw = slices[i].GetRawPtr();
        float* reconstructed = slices[i].GetReconstructedPtr( );
        float f = toolbox.CompareArrays( raw, reconstructed, NX*NY, true );
        rms.push_back( f );
    }
    cerr << "Overall RMS: " << toolbox.CalcRMS( rms ) << endl;
    

    if( filenames )                         delete[] filenames;
}
