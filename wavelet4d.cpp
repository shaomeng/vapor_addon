#include "wavelet4d.h"


using std::string;
using std::cout;
using std::cerr;
using std::endl;


Wavelet4D::Wavelet4D( long NX, long NY, long NZ, long NT )
{

	_BLOCKDIM = 64;
	_wavename = "bior4.4";
	_cratio = 1;

	_NX = NX;
	_NY = NY;
	_NZ = NZ;
	_NT = NT;

	_block_indices = NULL;
	_path.clear();
	CalcBlockIndices();
}

Wavelet4D::~Wavelet4D()
{
	if( _block_indices ) {
		delete[] _block_indices;
		_block_indices = NULL;
	}
}

void 
Wavelet4D::SetPath( const std::string &path )
{
	_path.clear();
	_path = path;
}

void
Wavelet4D::GenerateFilenames( const string &name, long idx)
{
	assert( !_path.empty());
	char buf[256];
	_filenames.clear();

	for( long i = 0; i < _NT; i++ ) {
		sprintf( buf, ".%04ld.out", i+idx);
		string f = _path + "/" + name + buf;
		_filenames.push_back( f );
	}	
}

void
Wavelet4D::CalcBlockIndices()
{
	assert( _NX % _BLOCKDIM == 0 );
	assert( _NY % _BLOCKDIM == 0 );
	assert( _NZ % _BLOCKDIM == 0 );
	
	long x_block = _NX / _BLOCKDIM;				// number of blocks in X
    long y_block = _NY / _BLOCKDIM;				// number of blocks in X
    long z_block = _NZ / _BLOCKDIM;				// number of blocks in X
    _BLOCKNUM = x_block * y_block * z_block;   

	_block_indices = new long[ _BLOCKNUM * 6 ];
	long counter = 0;

	/*
	 * 6 indices for a cell are stored together.
	 */
	for( long k = 0; k < z_block; k++ )
		for( long j = 0; j < y_block; j++ )
			for( long i = 0; i < x_block; i++ )
			{
                _block_indices[ counter++ ] = _BLOCKDIM *  i;		// startX
                _block_indices[ counter++ ] = _BLOCKDIM * (i+1);	// endX
                _block_indices[ counter++ ] = _BLOCKDIM *  j;		// startY
                _block_indices[ counter++ ] = _BLOCKDIM * (j+1);	// endY
                _block_indices[ counter++ ] = _BLOCKDIM *  k;		// startZ
                _block_indices[ counter++ ] = _BLOCKDIM * (k+1);	// endZ
			}
}

void
Wavelet4D::PrintBlockIndices()
{
	if( _block_indices == NULL )
		cerr << "_block_indices == NULL " << endl;
	else
		for( long i = 0; i < _BLOCKNUM; i++ )
			printf("%lu->%lu, %lu->%lu, %lu->%lu\n", 
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ], 
					_block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ], 
					_block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );
		
}

void
Wavelet4D::PrintFilenames()
{
	for( unsigned int i = 0; i < _filenames.size(); i++ )
		cout << _filenames[i] << endl;
}



int
Wavelet4D::ParallelExec()
{
#ifdef EVALUATE
	/*
  	 * RMS and LMAX with same block index but different _NT are stored together
	 */
	double* rms3d = new double[ _BLOCKNUM * _NT ];
	double* lmax3d = new double[ _BLOCKNUM * _NT ];
	double* rms4d = new double[ _BLOCKNUM * _NT ];
	double* lmax4d = new double[ _BLOCKNUM * _NT ];
#endif

	/*
	 * Each thread takes care of one index from _BLOCKNUM.
	 */
	int nthreadSys = omp_get_num_threads();
	int nthreads = (nthreadSys < _BLOCKNUM)? nthreadSys : _BLOCKNUM;
	omp_set_num_threads( nthreads );
	#pragma omp parallel for schedule( dynamic )
	for( long i = 0; i < _BLOCKNUM; i++ )
	{
		Cube3D** slices = new Cube3D*[ _NT ];
		SliceGroup* group = new SliceGroup( _wavename );

		for( long t = 0; t < _NT; t++ )
		{
			slices[t] = new Cube3D( _filenames[t], _wavename, 
					_BLOCKDIM, _BLOCKDIM, _BLOCKDIM, _NX, _NY, _NZ,
					_block_indices[ 6*i ],   _block_indices[ 6*i+1 ],
                    _block_indices[ 6*i+2 ], _block_indices[ 6*i+3 ],
                    _block_indices[ 6*i+4 ], _block_indices[ 6*i+5 ] );
			slices[t] -> Decompose();

#ifdef EVALUATE
			/*
			 * individual slice reconstruction and evaluation
			 */
			slices[t] -> Reconstruct (_cratio);
			slices[t] -> Evaluate( rms3d[ i*_NT + t ], lmax3d[ i*_NT + t ] );
			slices[t] -> ReloadInputFile();
			slices[t] -> Decompose();
#endif
			group -> AddSlice( slices[t] );
		}							
		/*
		 * Temporal compression
		 */
		group -> Initialize();
		group -> Decompose();		// All coefficients stored now!

		/*
		 * Filename for output coefficients
 		 */
		string coeff_name = _filenames[0] + ".block" + to_string(i) + ".coeff";
		group -> OutputFile( coeff_name, _cratio );

#ifdef EVALUATE
		group -> Reconstruct( _cratio );
		group -> UpdateSlices();
		for( long t = 0; t < _NT; t++ ) {
			slices[t] -> Reconstruct(1);
			slices[t] -> Evaluate( rms4d[ i*_NT + t ], lmax4d[ i*_NT + t ] );
		}
#endif

		if( group )				delete group;
		for( long t = 0; t < _NT; t++ )
			if( slices[t] ) 	delete slices[t];
		if( slices ) 			delete[] slices;
	}

#ifdef EVALUATE
	printf("3D DWT RMS = %e, LMAX = %e\n", FindRMS( rms3d, _BLOCKNUM * _NT ),
										   FindMax( lmax3d, _BLOCKNUM * _NT));
	printf("4D DWT RMS = %e, LMAX = %e\n", FindRMS( rms4d, _BLOCKNUM * _NT ),
										   FindMax( lmax4d, _BLOCKNUM * _NT));
	delete[] rms3d;
	delete[] lmax3d;
	delete[] rms4d;
	delete[] lmax4d;
#endif

	return 0;

}

double 
Wavelet4D::FindMax( const double* arr, long len ) {
    double max = 0;
    for( long i = 0; i < len; i++ )
        if( arr[i] > max )
            max = arr[i];
    return max;
}

double 
Wavelet4D::FindRMS( const double* arr, long len)
{
    double sum = 0.0;
    double c = 0.0;
    for( long i = 0; i < len; i++ ) {
        double y = arr[i] * arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum /= double(len);
    sum = sqrt( sum );
    return sum;
}

int main(int argc, char* argv[] )
{
	int cratio = 1;
	if( argc == 2 )
		cratio = atoi( argv[1] );
	Wavelet4D wav( 128, 128, 128, 5);
	wav.SetCRatio( cratio );
	string path = "/home/users/samuelli/Datasets/HD_128";
	string name = "vz";
	wav.SetPath( path );
	wav.GenerateFilenames( name, 1 );
	
	//wav.PrintFilenames();
	wav.StartMonitor();
	//wav.ParallelExec();

}

#ifdef INOTIFY
bool
Wavelet4D::FinishWriteAll(struct inotify_event *event)
{
	string last = _filenames.back();
	string current = _path + "/" + event->name;
	printf("%s\n", event->name );
	if( (event->mask & IN_CLOSE_WRITE) && (current.compare(last)==0) ) 
	{
		/* inspect size */
		long size_should = _NX * _NY * _NZ * 4;
		struct stat buffer;
		int status = stat(current.c_str(), &buffer);
		
		if( status==0 && size_should==(long)buffer.st_size )
		{
			/* make sure all files exist and have the correct size */
			for( unsigned int i = 0; i < _filenames.size()-1; i++ )
			{
				status = stat(_filenames[i].c_str(), &buffer);
				if( status!=0 || size_should!=(long)buffer.st_size )
				{
					cerr << "Error in file: " << _filenames[i] << endl;
					return false;
				}
			}
			return true;
		}
	}
	return false;
}

void
Wavelet4D::StartMonitor()
{
	assert( !_path.empty() );
	
    int inotifyFd, wd;
    char buf[BUF_LEN]; 
    long numRead;
    char *p;
    struct inotify_event *event;


    inotifyFd = inotify_init();                 /* Create inotify instance */
    if (inotifyFd == -1)
        perror("inotify_init() failed\n");

    /* Add a watch for _path: only watch CLOSE_WRITE */
	wd = inotify_add_watch( inotifyFd, _path.c_str(), IN_CLOSE_WRITE ); 
	if (wd == -1)
	{
		perror("inotify_add_watch failed: ");
		perror(_path.c_str());
	}
	printf("Watching %s using wd %d\n", _path.c_str(), wd);

	bool stop = false;
    while(!stop) {                                  /* Read events forever */
        numRead = read(inotifyFd, buf, BUF_LEN);
		if( numRead <= 0 )
			perror("inotify read return value error!");

        /* Process all of the events in buffer returned by read() */
        for (p = buf; p < buf + numRead; ) {
            event = (struct inotify_event *) p;
			if( FinishWriteAll(event) )
			{
				printf("finish write all files in this group!\n");
				stop = true;
			}
			if( event->len > 0 && (strcmp( event->name, "stop") == 0) )
				stop = true;

            p += sizeof(struct inotify_event) + event->len;
        }

    }
	inotify_rm_watch( inotifyFd, wd );
	close( inotifyFd );		/* close file discriptor */
}
#endif
