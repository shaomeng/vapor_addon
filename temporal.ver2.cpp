/*
 * This is the main file for performing DWT, using cube3d and slicegroup.
 * This version compares the levels of DWTs applied on time dimension.
 *
 * Programmer: Samuel Li
 * Date: 7/21/2015
 *
 */

#include "cube3d.h"
#include "slicegroup.h"

#define NX 128
#define NY 128
#define NZ 128
#define NSLICES 19


double FindMax( const double* arr, size_t len ) {
    double max = 0;
    for( size_t i = 0; i < len; i++ )
        if( arr[i] > max )
            max = arr[i];
    return max;
}

double CalcRMS( const double* arr, size_t len)
{
    double sum = 0.0;
    double c = 0.0;
    for( size_t i = 0; i < len; i++ ) {
        double y = arr[i] * arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= double(len);
    sum = sqrt( sum );
    return sum;
}

int main( int argc, char* argv[] )
{
    if( argc != 4 ) {
        cerr << "usage: ./a.out varname startIdx ratio! " << endl;
        exit(1);
    }
    string varname  = argv[1];
    int startIdx    = atoi( argv[2] );
    int ratio       = atoi( argv[3] );
    string wavenameXYZ = "bior4.4";
    string wavenameT = "bior4.4";
    string* filenames = new string[ NSLICES ];
    double rmsArr[ NSLICES ];
    double lmaxArr[ NSLICES ];

    
    /* 
     * Initialize filenames
     */
    char namebuf[64];
    string prefix = "/home/samuel/Datasets/Ghost/250/" + varname + ".";
    for( int i = 0; i < NSLICES; i++ )
    {
        sprintf( namebuf, "%04d.float", i + startIdx );
        filenames[i] = prefix + namebuf;
        if( i == 0 )            printf( "Start sample idx = %s\n", namebuf );
        if( i == NSLICES-1 )    printf( " End  sample idx = %s\n", namebuf );
    }

    Cube3D** slices = new Cube3D*[ NSLICES ];
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] = new Cube3D( filenames[i], wavenameXYZ, NX, NY, NZ );
    }


    /*
     * Experiment 1: two levels of DWT 
     */

     /* 3D compression on indivisual files.  */ 
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] -> Decompose();
    }

    /* Then apply temporal compression.  */
    SliceGroup group( wavenameT );
    for( int i = 0; i < NSLICES; i++ )
        group.AddSlice( slices[i] );
    group.Initialize();
    group.Decompose( );
    group.Reconstruct( ratio );
    group.UpdateSlices();

    /* Reconstruct on individual slices, and then evaluate.  */
    for( int i = 0; i < NSLICES; i++ ){
        slices[i] -> Reconstruct( 1 );
        slices[i] -> Evaluate( rmsArr[i], lmaxArr[i] );
    }
    double lmax = FindMax( lmaxArr, NSLICES );
    double rms = CalcRMS( rmsArr, NSLICES );
    cout << "==> 3D+1D compression (two levels DWT) at ratio " << ratio << endl;
    cout << "\tRMS, MAX:  " << rms << ", " << lmax << endl;



    for( int i = 0; i < NSLICES; i++ )
    {
        rmsArr[i] = 0.0;
        lmaxArr[i] = 0.0;
    }


    /* 
     * Experiment 2: one level of DWT.
     */

    /* Step 1: one level of DWT * on first halves of the data */
    /* 3D compression on individual files */
    for( int i = 0; i < NSLICES / 2; i++ ) {
        slices[i] -> ReloadInputFile();
        slices[i] -> Decompose();
    }
    /* Then apply temporal compression */
    SliceGroup group1( wavenameT );
    for( int i = 0; i < NSLICES / 2; i++ )
        group1.AddSlice( slices[i] );
    group1.Initialize();
    group1.Decompose( );
    group1.Reconstruct( ratio );
    group1.UpdateSlices();
    /* Reconstruct individual slices */
    for( int i = 0; i < NSLICES / 2; i++ ){
        slices[i] -> Reconstruct( 1 );
        slices[i] -> Evaluate( rmsArr[i], lmaxArr[i] );
    }


    /* Step 2: one level of DWT on second halves of the data */
    /* 3D compression on individual files */
    for( int i = NSLICES / 2; i < NSLICES; i++ ) {
        slices[i] -> ReloadInputFile();
        slices[i] -> Decompose();
    }
    /* Then apply temporal compression */
    SliceGroup group2( wavenameT );
    for( int i = NSLICES / 2; i < NSLICES; i++ )
        group2.AddSlice( slices[i] );
    group2.Initialize();
    group2.Decompose( );
    group2.Reconstruct( ratio );
    group2.UpdateSlices();
    /* Reconstruct individual slices */
    for( int i = NSLICES / 2; i < NSLICES; i++ ){
        slices[i] -> Reconstruct( 1 );
        slices[i] -> Evaluate( rmsArr[i], lmaxArr[i] );
    }


    /*
     * Evaluate the results from two halves
     */
    lmax = FindMax( lmaxArr, NSLICES );
    rms = CalcRMS( rmsArr, NSLICES );
    cout << "==> 3D+1D compression (two one-level DWTs) at ratio " << ratio << endl;
    cout << "\tRMS, MAX:  " << rms << ", " << lmax << endl;
        
    
    for( int i = 0; i < NSLICES; i++ )
        if( slices[i] )         delete slices[i];    
    delete[]                    slices;
    if( filenames )             delete[] filenames;

}
