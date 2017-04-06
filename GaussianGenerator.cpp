#include "cstdlib"
#include "cstdio"
#include "cmath"

#define FLOAT float

void WriteGaussian2D( long int dimX, long int dimY, char* filename )
{
  // 2D gaussian parameters
  double amp     = 6400.0;
  double x0      = (double)dimX / 2.0;  // center
  double y0      = (double)dimY / 2.0;  // center
  double sigmaX  = (double)dimX / 4.0;  // spread
  double sigmaY  = (double)dimY / 4.0;  // spread
  double sigmaX2 = sigmaX * sigmaX * 2.0;
  double sigmaY2 = sigmaY * sigmaY * 2.0;
  const double EulerConstant = exp(1.0);
  
  FLOAT* buffer = new FLOAT[dimX * dimY];

  for( long int i = 0; i < dimX * dimY; i++ )
  {
    long int x = i % dimX;
    long int y = i / dimX;
    double power = (x-x0) * (x-x0) / sigmaX2 + (y-y0) * (y-y0) / sigmaY2;
    double result = pow( EulerConstant, power * -1.0 );
    buffer[i] = (FLOAT)( result * amp );
  }

  FILE* outFile = fopen( filename, "wb" );
  fwrite( buffer, sizeof(FLOAT), dimX * dimY, outFile );
  fclose( outFile );

  delete[] buffer;

}

void WriteGaussian3D( long int dimX, long int dimY, long int dimZ, char* filename )
{
  FLOAT amp = 10.0;
  FLOAT sigmaX  = (FLOAT)dimX / 4.0;  
  FLOAT sigmaY  = (FLOAT)dimY / 4.0;  
  FLOAT sigmaZ  = (FLOAT)dimZ / 4.0;  
  FLOAT sigmaX2 = sigmaX * sigmaX * 2.0;
  FLOAT sigmaY2 = sigmaY * sigmaY * 2.0;
  FLOAT sigmaZ2 = sigmaZ * sigmaZ * 2.0;
  FLOAT EulerConstant = exp(1.0);

  long int totalLen = dimX * dimY * dimZ;
  FLOAT* buffer = new FLOAT[ totalLen ];
  for( long int i = 0; i < totalLen; i++ )
  {
    long int z = i / (dimX * dimY);
    long int y = (i - z * dimX * dimY) / dimX;
    long int x = i % dimX;
    FLOAT   xf = (FLOAT)x - (FLOAT)dimX / 2;
    FLOAT   yf = (FLOAT)y - (FLOAT)dimY / 2;
    FLOAT   zf = (FLOAT)z - (FLOAT)dimZ / 2;
    buffer[i]  = amp * pow( EulerConstant, -1.0 * (xf*xf/sigmaX2 + yf*yf/sigmaY2 + zf*zf/sigmaZ2) );
  }  

  FILE* outFile = fopen( filename, "wb" );
  fwrite( buffer, sizeof(FLOAT), dimX * dimY * dimZ, outFile );
  fclose( outFile );

  delete[] buffer;
}

void WriteIncrement( long int totalLen, char* filename )
{
  FLOAT* buffer = new FLOAT[ totalLen ];
  for( long int i = 0; i < totalLen; i++ )
    buffer[i] = (FLOAT)i; 
  FILE* outFile = fopen( filename, "wb" );
  fwrite( buffer, sizeof(FLOAT), totalLen, outFile );
  fclose( outFile );
  delete[] buffer;
}

int main(int argc, char** argv)
{
  if( argc != 5 )
  {
    printf("Usage: ./a.out  dimX  dimY  dimZ  output_filename\n");
    return 1;
  }
  long int dimX   = atol( argv[1] ); 
  long int dimY   = atol( argv[2] ); 
  long int dimZ   = atol( argv[3] ); 
  char* filename  = argv[4];

  WriteGaussian3D( dimX, dimY, dimZ, filename );
  //WriteIncrement( dimX * dimY * dimZ, filename );

  return 0;
}
