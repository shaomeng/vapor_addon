#include "cstdlib"
#include "cstdio"
#include "cmath"

#define FLOAT double

int main(int argc, char** argv)
{
  if( argc != 4 )
  {
    printf("Usage: ./a.out dimX dimY output_filename\n");
    return 1;
  }
  long int dimX   = atol( argv[1] ); 
  long int dimY   = atol( argv[2] ); 
  char* filename  = argv[3];

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

  return 0;
}
