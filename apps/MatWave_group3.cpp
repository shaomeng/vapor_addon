#include "vapor/SamSlice2.h"
#include "vapor/SamSliceGroup3.h"
#include "vapor/SamToolbox.h"
#include <sys/time.h>

#define NX 504
#define NY 504
#define NZ 2048
#define NSLICES 10

using namespace VAPoR;

int main( int argc, char* argv[] )
{
    vector< int > cratios;
    if( argc < 4 ) {
        cerr << "usage: ./a.out UVW 0 ratio1 ratio2 ... " << endl;
        exit(1);
    }
    else
        for( int i = 3; i < argc; i++ )
            cratios.push_back( stoi( argv[i] ) );

    string varname  = argv[1];
    int startIdx    = stoi( argv[2] );
    string wavenameXYZ = "bior4.4";
    string wavenameT = "bior1.1";
    string* filenames = new string[ NSLICES ];
    struct timeval t_now;
    double timer1, timer2;

    
    string path = "/work/02892/samuelli/maverick/15plume3d/float/" + 
                  varname + "/";
    for( long long i = 0; i < NSLICES; i++ )
        filenames[i] = path + varname + "." + to_string(startIdx + i) + 
                        ".float";

    // Use replicated data files
    // string path = "/work/02892/samuelli/maverick/15plume3d/float/" + 
    //                 varname + "_replica/";

    // Use the 504*504*1024 files
    //string path = "/work/03546/tg828451/maverick/15plume3d_1024z/";
    //for( long long i = 0; i < NSLICES; i++ )
    //    filenames[i] = path + "U1024z." + to_string(startIdx + i) + ".float";



    SamToolbox toolbox;
    SamSlice2** slices = new SamSlice2*[ NSLICES ];
    SamErr errs[ cratios.size() ][ NSLICES ];

    #pragma omp parallel for
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] = new SamSlice2( filenames[i], wavenameXYZ, NX, NY, NZ );
        slices[i] -> Decompose();
    }

    cerr << "==> Variable " << varname << ", from time stamp " << startIdx+400 << endl;

    /* 
     * Test on different cratios, 3D compression
     */
/*
    #pragma omp parallel for
    for( int i = 0; i < NSLICES; i++ ) 
    {
        for( int j = 0; j < cratios.size(); j++ )
        {
            int ratio = cratios[j];
            slices[i] -> Reconstruct( ratio );
            float* rawPtr = slices[i] -> GetRawPtr();
            float* reconstructedPtr = slices[i] -> GetReconstructedPtr();
            errs[j][i] = toolbox.CompareArrays( rawPtr, reconstructedPtr, NX*NY*NZ, false );
        }
        slices[i] -> FreeReconstructed();   // make space for group operation
    }

    cerr << "==> 3D compression " << endl;
    cerr << "\tRatio,\t\tRMS,\t\t\tMAX" << endl;
    for( int i = 0; i < cratios.size(); i++ ){
        float max = toolbox.FindMax( errs[i], NSLICES );
        float rms = toolbox.CalcRMS( errs[i], NSLICES );
        cerr << "\t" << cratios[i] << ",\t\t" << rms << ",\t\t" << max << endl;
    }
*/


    /*
     * Setup Temporal Compression
     */ 
    SamSliceGroup3* group = new SamSliceGroup3( wavenameT, NX*NY*NZ, NSLICES );
    vector< float* > rawarr;
    for( int i = 0; i < NSLICES; i++ )
        rawarr.push_back( slices[i] -> GetCoeffsPtr() );
    group -> UpdateRawPtr( rawarr );
    for( int i = 0; i < NSLICES; i++ )
        slices[i] -> FreeCoeffs();


    //gettimeofday(&t_now, NULL);
    //timer1 = t_now.tv_sec + (t_now.tv_usec/1000000.0);
    group -> Decompose();
    //gettimeofday(&t_now, NULL);
    //timer2 = t_now.tv_sec + (t_now.tv_usec/1000000.0);
    //cerr << "\tfinish temporal decomposing in " << (timer2 - timer1) << endl;


    /*
     * Test on different cratios, 3D + 1D compression
     */
    for( int i = 0; i < cratios.size(); i++ )
    {
        group -> Reconstruct( cratios[i] );

        #pragma omp parallel for
        for( int j = 0; j < NSLICES; j++ )
        {
            group -> FillReconstructedPtr( j, slices[j] -> GetCoeffsPtr() );
            slices[j] -> Reconstruct( 1 );
            float* rawPtr           = slices[j] -> GetRawPtr();
            float* reconstructedPtr = slices[j] -> GetReconstructedPtr();
            errs[i][j] = toolbox.CompareArrays( rawPtr, reconstructedPtr, NX*NY*NZ, false );
        }
    }

    cerr << "==> 3D + 1D results:" << endl;
    cerr << "\tRatio,\t\tRMS,\t\t\tMAX" << endl;
    for( int i = 0; i < cratios.size(); i++ ){
        float max = toolbox.FindMax( errs[i], NSLICES );
        float rms = toolbox.CalcRMS( errs[i], NSLICES );
        cerr << "\t" << cratios[i] << ",\t\t" << rms << ",\t\t" << max << endl;
    }





    if( group )                 delete group;

    for( int i = 0; i < NSLICES; i++ )
        if( slices[i] )         delete slices[i];    
    delete[]                    slices;
    if( filenames )             delete[] filenames;

}
