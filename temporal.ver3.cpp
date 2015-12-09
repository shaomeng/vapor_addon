/*
 * This is the main file for performing DWT, using cube3d and slicegroup.
 * This version compares evaluates how the sparsity of time slices affect accuracy.
 *
 * Programmer: Samuel Li
 * Date: 7/24/2015
 *
 */

#include "cube3d.h"
#include "slicegroup.h"
#include <sys/time.h>

#define NX 64
#define NY 64
#define NZ 64
#define NSLICES 20


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
    struct timeval t_now;
    double timer1, timer2;
    double rmsArr[ NSLICES ];
    double lmaxArr[ NSLICES ];

    
    /* 
     * Test on the small test data set.
     */
    string prefix = "/home/samuel/Git/utilities/small_test_datasets/plume_64/" 
                    + varname + ".";
    for( int i = 0; i < NSLICES; i++ )
        filenames[i] = prefix + to_string(startIdx + i) + ".float";

    Cube3D** slices = new Cube3D*[ NSLICES ];
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] = new Cube3D( filenames[i], wavenameXYZ, NX, NY, NZ );
    }


    /* 
     * 1st experiment: 2 levels of temporal DWT 
     */
    /* 3D compression on indivisual files. */ 
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] -> Decompose();
    }
    /* Then apply temporal compression. */
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
    cerr << "==> 3D+1D compression (two levels DWT) at ratio " << ratio << endl;
    cerr << "\tRMS, MAX:  " << rms << ", " << lmax << endl;


    /* Initialization */
    for( int i = 0; i < NSLICES; i++ )
    {
        rmsArr[i] = 0.0;
        lmaxArr[i] = 0.0;
    }


    /*  
     * 2nd experiment: one level of DWT 
     * on a sparser version of the data set.
     * More specifically, use every other time slice, 
     * but less aggresive compression.
     */
    /* 3D compression on individual files */
    for( int i = 0; i < NSLICES; i += 2 ) {
        slices[i] -> ReloadInputFile();
        slices[i] -> Decompose();
    }
    /* Then apply temporal compression */
    SliceGroup group1( wavenameT );
    for( int i = 0; i < NSLICES; i += 2 )
        group1.AddSlice( slices[i] );
    group1.Initialize();
    group1.Decompose( );
    group1.Reconstruct( ratio * 2 );    // less aggresive compression here
    group1.UpdateSlices();
    /* Reconstruct individual slices */
    for( int i = 0; i < NSLICES; i += 2 ){
        slices[i] -> Reconstruct( 1 );
        slices[i] -> Evaluate( rmsArr[i], lmaxArr[i] );
    }
    /* put the results into new arrays */
    double  rmsArr1[ NSLICES / 2 ];
    double lmaxArr1[ NSLICES / 2 ];
    for( int i = 0; i < NSLICES; i += 2 )
    {
        rmsArr1[  i/2 ] = rmsArr[ i ];
        lmaxArr1[ i/2 ] = lmaxArr[ i ];
    }


    /*
     * Evaluate the results from sparser slices.
     */
    lmax = FindMax( lmaxArr1, NSLICES/2 );
    rms = CalcRMS( rmsArr1, NSLICES/2 );
    cerr << "==> 3D+1D compression (on sparser slices) at ratio " << ratio << endl;
    cerr << "\tRMS, MAX:  " << rms << ", " << lmax << endl;
        
    
    for( int i = 0; i < NSLICES; i++ )
        if( slices[i] )         delete slices[i];    
    delete[]                    slices;
    if( filenames )             delete[] filenames;

}
