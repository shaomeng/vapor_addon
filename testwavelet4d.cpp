#include "wavelet4d.h"
#include <string>

using std::string;

int main(int argc, char* argv[] )
{
	int startIdx = 8502;
	int increment = 2;
	int total_steps = 500;
	int each_group = 18;
	int cratio = 8;
	string name = "vortmag";
	char d = 'A';
	if( argc == 2 )
		cratio = atoi( argv[1] );
	if( argc == 3 )
	{
		cratio = atoi( argv[1] );
		name = argv[2];
	}
	if( argc == 4 )
	{
		cratio = atoi( argv[1] );
		name = argv[2];
		d = argv[3][0];
	}
	Wavelet4D wav( 490, 490, 280, each_group);
	wav.SetCRatio( cratio );
	string path = "/home/users/samuelli/Datasets/Orf/" + name + "_cropped";
	//string pathCompare = "/home/users/samuelli/Datasets/HD512_500/enstrophy";
	wav.SetPath( path );

	for( int i = 0; i < total_steps / each_group; i++ )
	{
		wav.GenerateFilenames( name, i * increment * each_group + startIdx );
		//wav.GenerateBkpFilenames( pathCompare, name, i * each_group + startIdx );
		wav.PrintFilenames();
		cout << "start Idx: " << i * each_group + startIdx << endl;
		wav.StartMonitor();
		if( argc < 4 )
		{
			wav.Output3DReconstruct();
			wav.Output4DReconstruct();
		}
		else 
		{
			if	   ( d == 'X' )		wav.XDimParallelExec();
			else if( d == 'Y' ) 	wav.YDimParallelExec();
			else if( d == 'Z' ) 	wav.ZDimParallelExec();
			else if( d == 'T' )		wav.TimeDimParallelExec();
			else					wav.ParallelExec();
		}
	}

}
