#include "wavelet4d.h"

#include <string>

using std::string;

int main(int argc, char* argv[] )
{
	int startIdx = 2;
	int total_steps = 198;
	int each_group = 9;
	int cratio = 8;
	string name = "viscosity";
	if( argc == 2 )
		cratio = atoi( argv[1] );
	if( argc == 3 )
	{
		cratio = atoi( argv[1] );
		name = argv[2];
	}
	Wavelet4D wav( 96, 96, 96, each_group);
	wav.SetCRatio( cratio );
	string path = "/home/users/samuelli/Datasets/cloverleaf/viscosity";
	string pathCompare = "/home/users/samuelli/Datasets/Hurricane_w_special_values";
	wav.SetPath( path );


	for( int i = 0; i < total_steps / each_group; i++ )
	{
		wav.GenerateFilenames( name, i * each_group + startIdx );
		wav.GenerateBkpFilenames( pathCompare, name, i * each_group + startIdx );
		//wav.PrintFilenames();
		cout << "start Idx: " << i * each_group + startIdx << endl;
		wav.StartMonitor();
		wav.ParallelExec();
	}
	

}
