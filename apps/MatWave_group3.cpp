#include "vapor/SamSlice2.h"
#include "vapor/SamSliceGroup3.h"
#include "vapor/SamSliceGroup2.h"
#include "vapor/SamToolbox.h"
#include <sys/time.h>

#define NX 504
#define NY 504
#define NZ 504
#define NSLICES 9

using namespace VAPoR;

int main( int argc, char* argv[] )
{
    if( argc < 2 ) {
        cerr << "Please specify the compression ratios to test!" << endl;
        exit(1);
    }
    int ratio = stoi( argv[1] );

    struct timeval t_now;
    double timer1, timer2;

    string* filenames = new string[ NSLICES ];
    string path = "/work/02892/samuelli/maverick/15plume3d/float/plume_504cube/";
    for( long long i = 0; i < NSLICES; i++ )
        filenames[i] = path + to_string(i) + ".float";

/*
    vector< int > cratios;
    for( int i = 1; i < argc; i++ )
        cratios.push_back( stoi( argv[i] ) );
*/
    string wavename = "bior4.4";

    SamToolbox toolbox;
    SamSlice2** slices = new SamSlice2*[ NSLICES ];

    float rms[ NSLICES ];
    #pragma omp parallel for
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] = new SamSlice2( filenames[i], wavename, NX, NY, NZ );
        slices[i] -> Decompose();
        slices[i] -> Reconstruct( ratio );
        float* rawPtr = slices[i] -> GetRawPtr();
        float* reconstructedPtr = slices[i] -> GetReconstructedPtr();
        rms[i] = toolbox.CompareArrays( rawPtr, reconstructedPtr, NX*NY*NZ, false );
        slices[i] -> FreeReconstructed();   // make space for group operation
    }
    for( int i = 0; i < NSLICES; i++ )
        cerr << "Slice# " << i << ", RMS=" << rms[i] << endl;
    cerr << "\t==> 3D RMS: " << toolbox.CalcRMS( rms, NSLICES ) << endl;


    SamSliceGroup3* group = new SamSliceGroup3( "bior4.4", NX*NY*NZ, NSLICES );
    vector< float* > rawarr;
    for( int i = 0; i < NSLICES; i++ )
        rawarr.push_back( slices[i] -> GetCoeffsPtr() );
    group -> UpdateRawPtr( rawarr );

    gettimeofday(&t_now, NULL);
    timer1 = t_now.tv_sec + (t_now.tv_usec/1000000.0);
    group -> Decompose();
    gettimeofday(&t_now, NULL);
    timer2 = t_now.tv_sec + (t_now.tv_usec/1000000.0);
    cerr << "\tfinish decomposing in " << (timer2 - timer1) << " seconds, now reconstructing the group..." << endl;

    gettimeofday(&t_now, NULL);
    timer1 = t_now.tv_sec + (t_now.tv_usec/1000000.0);
    group -> Reconstruct( ratio );
    gettimeofday(&t_now, NULL);
    timer2 = t_now.tv_sec + (t_now.tv_usec/1000000.0);
    cerr << "\tfinish reconstructing in " << (timer2 - timer1) << " seconds, now reconstructing slices..." << endl;




/*
    SamSliceGroup2* group2 = new SamSliceGroup2( "bior4.4", NX*NY*NZ, NSLICES );
    for( int i = 0; i < NSLICES; i++ )
        group2 -> UpdateRawPtr( i, slices[i] -> HandoverCoeffs() );
    group2 -> Decompose();
    group2 -> Reconstruct( ratio );

    float* g2_buf = new float[ NX*NY*NZ ];
    g2_buf = group2 -> GetReconstructedPtr( NSLICES-1 );

    float* g3_buf = new float[ NX*NY*NZ ];
    group -> FillReconstructedPtr( NSLICES-1, g3_buf );

    for( size_t i = 0; i < NX*NY*NZ; i++ )
        if( g2_buf[i] != g3_buf[i] )
            cerr << i << "\t" << g2_buf[i] << "\t" << g3_buf[i] << endl;

    delete[]    g2_buf;
    delete[]    g3_buf;
*/







    #pragma omp parallel for
    for( int i = 0; i < NSLICES; i++ ){
        group -> FillReconstructedPtr( i, slices[i] -> GetCoeffsPtr() );
        slices[i] -> Reconstruct( 1 );   
    }

    for( int i = 0; i < NSLICES; i++ ) {
        float* rawPtr = slices[i] -> GetRawPtr();
        float* reconstructedPtr = slices[i] -> GetReconstructedPtr();
        assert( rawPtr != NULL );
        assert( reconstructedPtr != NULL );
        rms[i] = toolbox.CompareArrays( rawPtr, reconstructedPtr, NX*NY*NZ, false );
    }
    for( int i = 0; i < NSLICES; i++ )
        cerr << "Slice# " << i << ", RMS=" << rms[i] << endl;
    cerr << "\t==> 3D+1D RMS: " << toolbox.CalcRMS( rms, NSLICES ) << endl;

    



    if( group )                 delete group;
    for( int i = 0; i < NSLICES; i++ )
        if( slices[i] )         delete slices[i];    
    delete[]                    slices;
    if( filenames )             delete[] filenames;

}
