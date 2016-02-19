#include "wavelet4d.h"

#include <string>

using std::string;

int main(int argc, char* argv[] )
{
	int startIdx = 1;
	int total_steps = 90;
	int each_group = 18;
	int cratio = 8;
	string name = "xd";
	if( argc == 2 )
		cratio = atoi( argv[1] );
	if( argc == 3 )
	{
		cratio = atoi( argv[1] );
		name = argv[2];
	}
	Wavelet4D wav( 63, 63, 63, each_group);
	wav.SetCRatio( cratio );
	string path = "/home/users/samuelli/Datasets/lulesh_field";
	wav.SetPath( path );

	for( int i = 0; i < total_steps / each_group; i++ )
	{
		wav.GenerateFilenames( name, i * each_group + startIdx );
		//wav.PrintFilenames();
		cout << "start Idx: " << i * each_group + startIdx << endl;
		wav.StartMonitor();
		wav.ParallelExec();
	}
	

}
