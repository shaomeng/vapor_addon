/*
 * This is the main file for performing DWT, using cube3d and slicegroup.
 * This version compares 3D vs (3D + 1D)
 *
 * Programmer: Samuel Li
 * Date: 7/17/2015
 *
 */

#include "cube3d.h"
#include "slicegroup.h"

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

    

    /* 
     * Initialize filenames.
     */
    char namebuf[64];
    string prefix = "/flash_buffer/Sam/HD_128/" + varname + ".";
	//string prefix = "/home/users/samuelli/Git/utilities/small_test_datasets/plume_64/" + varname + ".";
    for( int i = 0; i < NSLICES; i++ )
    {
        sprintf( namebuf, "%04d.out", i + startIdx );
        filenames[i] = prefix + namebuf;
        if( i == 0 )            printf( "Start sample idx = %d\n", i + startIdx );
        if( i == NSLICES-1 )    printf( "End  sample idx = %d\n", i + startIdx );
    }

    Cube3D** slices = new Cube3D*[ NSLICES ];
    for( int i = 0; i < NSLICES; i++ ) {
//        slices[i] = new Cube3D( filenames[i], wavenameXYZ, NX, NY, NZ );
        slices[i] = new Cube3D( filenames[i], wavenameXYZ, NX, NY, NZ,
								128, 128, 128,
								0, 64, 64, 128, 0, 64 );
    }


    /* 
     * 3D compression on indivisual files, and evaluate.
     */
    double rmsArr[ NSLICES ];
    double lmaxArr[ NSLICES ];
    for( int i = 0; i < NSLICES; i++ ) 
    {
        slices[i] -> Decompose();
        slices[i] -> Reconstruct( ratio );
        slices[i] -> Evaluate( rmsArr[i], lmaxArr[i] );
    }
    double lmax = FindMax( lmaxArr, NSLICES );
    double rms = CalcRMS( rmsArr, NSLICES );
    cout  << "==> 3D compression at ratio " << ratio << endl;
    cout  << "\tRMS, MAX:  " << rms << ", " << lmax << endl;

	/*
	 * individual errors
  	 */
	for( int i = 0; i < NSLICES; i++ )
		cout << "RMS, LMAX: " << rmsArr[i] << "\t" << lmaxArr[i] << endl;


    /*
     * 3D compression on indivisual files.
     */ 
    for( int i = 0; i < NSLICES; i++ ) {
        slices[i] -> ReloadInputFile();
        slices[i] -> Decompose();
    }
    /*
     * Then apply temporal compression.
     */
    SliceGroup group( wavenameT );
    for( int i = 0; i < NSLICES; i++ )
        group.AddSlice( slices[i] );
    group.Initialize();
    group.Decompose( );
    group.Reconstruct( ratio );
    group.UpdateSlices();
    /*
     * Reconstruct on individual slices, and then evaluate.
     */
    for( int i = 0; i < NSLICES; i++ ){
        slices[i] -> Reconstruct( 1 );
        slices[i] -> Evaluate( rmsArr[i], lmaxArr[i] );
    }
    lmax = FindMax( lmaxArr, NSLICES );
    rms = CalcRMS( rmsArr, NSLICES );
    cout  << "==> 3D+1D compression at ratio " << ratio << endl;
    cout  << "\tRMS, MAX:  " << rms << ", " << lmax << endl;

	/*
	 * individual errors
  	 */
	for( int i = 0; i < NSLICES; i++ )
		cout << "RMS, LMAX: " << rmsArr[i] << "\t" << lmaxArr[i] << endl;
        
    
    for( int i = 0; i < NSLICES; i++ )
        if( slices[i] )         delete slices[i];    
    delete[]                    slices;
    if( filenames )             delete[] filenames;

}
