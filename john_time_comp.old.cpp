//
// Compress and decompress a 1D array
//
void test1d(std::string wavename, const float *srcarr, float *dstarr, float cratio ) {

	for (size_t i=0; i<NX*NY*NZ; i++) dstarr[i] = 0.0;

	VAPoR::MatWaveWavedec mw(wavename);

	// Compute # of wavelet decompositions. Same # along each 
	// dimension
	//
	size_t nlevels = mw.wmaxlev( NX*NY*NZ );

	size_t clen = mw.coefflength(NX* NY* NZ, nlevels);
	float *C = new float[clen];
	assert(C);

	size_t L[ nlevels + 2];

	// Compute the bookkeeping vector, L.
	//
	// N.B. L will be recomputed by wavedec(), but not
	// waverec();
	//
	mw.computeL(NX* NY* NZ, nlevels, L);

	//size_t startCoeffIdx = 0;
	//if (KeepAppCoeff) startCoeffIdx = L[0]*L[1]*L[2];
	int rc = mw.wavedec(srcarr, NX* NY* NZ, nlevels, C, L);
	assert (rc>=0);

  for( size_t i = 0; i < NX* NY* NZ; i++ )
    cout << C[i] << endl;

  size_t NumCoeff = (float) (NX*NY*NZ) / cratio;

	//
    // sort the coefficients  and find the threshold for culling
	// coefficients.
    //
   	vector <float> sortedC; 
    for (size_t i=0; i<clen; i++)  sortedC.push_back(C[i]);
    sort(sortedC.begin(), sortedC.end(), my_compare<float>);

	cout << "sortedC.size() = " << sortedC.size() << endl;

	size_t ti = NumCoeff - 1;

	float threshold  = sortedC[ti];

	// Zero out wavelet coefficients that are smaller than threshold
	//
	cout << "ti " << ti << endl;
	cout << "threshold  " << threshold  << endl;
    size_t should_squash = clen - clen / cratio;    // should squash this many
	size_t squashed = 0;
	for (size_t i=0; i<clen; i++) {
		if (fabs(C[i]) <= fabs(threshold) && squashed <= should_squash) {
			squashed += 1;
			C[i] = 0;
		}
	}
	cout << "squashed " << squashed << endl;
	
	// Inverse DWT
	//
	mw.waverec(C,L,nlevels,dstarr);	

    delete[] C;
}



// Test 2d dwt, and print out the coefficients
template< typename T >
void test2d(std::string wavename, const T *srcarr, T *dstarr, float cratio,
            size_t dimX, size_t dimY) 
{

	for (size_t i=0; i<dimX * dimY; i++) dstarr[i] = 0.0;

	// 2D transform first
	//
	VAPoR::MatWaveWavedec mw(wavename);

	size_t nlevels2d = min(mw.wmaxlev(dimX), mw.wmaxlev(dimY));

	size_t clen2d = mw.coefflength2(dimX, dimY, nlevels2d);
	assert (clen2d == dimX*dimY);

	// space for 2D wavelet coefficients
	//
	T* C2d = new T[clen2d];

	size_t L2d[(6*nlevels2d)+4];

	// Compute the bookkeeping vector, L2d.
	//
	mw.computeL2(dimX, dimY, nlevels2d, L2d);

	// perform 2D transforms, one for each of the NZ slices
	//
  int rc = mw.wavedec2(srcarr, dimX, dimY, nlevels2d, C2d, L2d);
  assert (rc>=0);

  // threshold coefficients
  //
  if( cratio > 1 )
  {
    vector<T> sortedC( clen2d, 0 );
    for( size_t i = 0; i < clen2d; i++ )
      sortedC[i] = C2d[i];
    int n = clen2d / cratio;
    //nth_element( sortedC.begin(), sortedC.begin() + n, sortedC.end(), my_compare<T> );
    sort(sortedC.begin(), sortedC.end(), my_compare<T>);
    T threshold = sortedC[n];
    if( threshold < 0 )
      threshold *= -1.0 ;
    for( size_t i = 0; i < clen2d; i++ )
      if( fabs( C2d[i] ) < threshold )
        C2d[i] = 0.0;
  }

  // reconstruct
  // 
  rc = mw.waverec2( C2d, L2d, nlevels2d, dstarr );
  assert( rc >= 0 );

  delete[] C2d;
}



// 
// Compress and decompress a 3D volume using a 2D wavelet transform
// transform along X and Y axes, followed by a 1D transform along
// the Z axis.
//
void test2dp1d(std::string wavename, const float *srcarr, float *dstarr, float cratio) 
{

	for (size_t i=0; i<NX*NY*NZ; i++) dstarr[i] = 0.0;

	// 
	// 2D transform first
	//
	VAPoR::MatWaveWavedec mw(wavename);

	size_t nlevels2d = min(mw.wmaxlev(NX), mw.wmaxlev(NY));

	size_t clen2d = mw.coefflength2(NX, NY, nlevels2d);
	assert (clen2d == NX*NY);

	// space for 2D wavelet coefficients
	//
	vector <float *> C2d;
	for (int i=0; i<NZ; i++) C2d.push_back(new float[clen2d]);

	size_t L2d[(6*nlevels2d)+4];

	// Compute the bookkeeping vector, L2d.
	//
	mw.computeL2(NX, NY, nlevels2d, L2d);

	// perform 2D transforms, one for each of the NZ slices
	//
	for (int i=0; i<NZ; i++) {
		int rc = mw.wavedec2(srcarr+i*NX*NY, NX, NY, nlevels2d, C2d[i], L2d);
		assert (rc>=0);
	}

	//
	// Now do 1D forward DWT
	//

	size_t nlevels1d = mw.wmaxlev(NZ);

	size_t clen1d = mw.coefflength(NZ, nlevels1d);
	assert(clen1d == NZ);

	// Space for coefficients resulting from NX*NY 1D transforms
	//
	vector <float *> buf1d, C1d;
	for (int i=0; i<NX*NY; i++) {
		C1d.push_back(new float[clen1d]);
		buf1d.push_back(new float[clen1d]);
	}

	size_t L1d[(nlevels1d)+2];
	assert(L1d);

	mw.computeL(NZ, nlevels1d, L1d);

	// Copy 2D coefficients to buf1d, tranpose array in process
	//
	for (int i=0; i<NX*NY; i++) {
	for (int j=0; j<NZ; j++) {
		buf1d[i][j] = C2d[j][i];
	}
	}


	// Forward 1D DWTs
	//
	for (int i=0; i<NX*NY; i++) {
		int rc = mw.wavedec(buf1d[i], NZ, nlevels1d, C1d[i], L1d);
		assert (rc>=0);
	}

	// Copy coefficients to single vector for easy sorting
	//
	vector <float> sortedC;
	for (int i=0; i<NX*NY; i++) {
		for (size_t j=0; j<clen1d; j++)  {
			sortedC.push_back(C1d[i][j]);
		}
	}

	cout << "sortedC.size() = " << sortedC.size() << endl;

	  //
    // sort the coefficients 
    //
    sort(sortedC.begin(), sortedC.end(), my_compare<float>);

    size_t NumCoeff = (float) (NX*NY*NZ) / cratio;
	size_t ti = NumCoeff - 1;

	double threshold  = sortedC[ti];

	// Kill coefficients that are smaller than threshold
	//
	cout << "ti " << ti << endl;
	cout << "threshold  " << threshold  << endl;
	size_t squashed = 0;
	for (int i=0; i<NX*NY; i++) {
	for (size_t j=0; j<clen1d; j++) {
		if (fabs(C1d[i][j]) <= fabs(threshold )) {
			squashed += 1;
			C1d[i][j] = 0;
		}
	}
	}
	cout << "squashed " << squashed << endl;

	//
	// Now perform IDWTs. First 1D, than 2D transforms
	//

	for (int i=0; i<NX*NY; i++) {
		int rc = mw.waverec(C1d[i], L1d, nlevels1d, buf1d[i]);
		assert (rc>=0);
	}

	// copy and tranpose coefficients
	//
	for (int i=0; i<NX*NY; i++) {
	for (int j=0; j<NZ; j++) {
		C2d[j][i] = buf1d[i][j];
	}
	}

	for (int i=0; i<NZ; i++) {
		int rc = mw.waverec2(C2d[i], L2d, nlevels2d, dstarr+i*NX*NY);
		assert (rc>=0);
	}
	
    for( size_t i = 0; i < C2d.size(); i++ )
        if( C2d[i] )    delete[] C2d[i];
    for( size_t i = 0; i < C1d.size(); i++ )
        if( C1d[i] )    delete[] C1d[i];
    for( size_t i = 0; i < buf1d.size(); i++ )
        if( buf1d[i] )    delete[] buf1d[i];
}
	
template< typename T >
void test2D_openmp( std::string wname, const T* srcarr, T* dstarr, float cratio, 
                    size_t dimX, size_t dimY, size_t nblocksX, size_t nblocksY )
{
  assert( dimX % nblocksX == 0 );
  assert( dimY % nblocksY == 0 );
  size_t blockDimX  = dimX / nblocksX;
  size_t blockDimY  = dimY / nblocksY;
  size_t nblocks    = nblocksX * nblocksY;
  
  // Calculate indices of each block
  //
  size_t* indices = new size_t[nblocks * 2];
  for(   size_t y = 0; y < nblocksY; y++ )
    for( size_t x = 0; x < nblocksX; x++ )
    {
      size_t idx = y * nblocksX + x;
      indices[idx * 2 + 0] = blockDimX * x;       // startX
      indices[idx * 2 + 1] = blockDimY * y;       // startY
    }

#ifdef OPENMP
  #pragma omp parallel for
#endif
  for( size_t i = 0; i < nblocks; i++ )
  {
#ifdef OPENMP
    if( i == 0 )
      cout << "number of threads: " <<  omp_get_num_threads() << endl;
#endif
    
    T* inBuf  = new T[ blockDimX * blockDimY ];
    T* outBuf = new T[ blockDimX * blockDimY ];
    RectangleCopy( srcarr, dimX, dimY,
                   inBuf, blockDimX, blockDimY,
                   indices[i*2 + 0], indices[i*2 + 1],
                   0, 0,
                   blockDimX, blockDimY ); 
    test2d( wname, inBuf, outBuf, cratio, blockDimX, blockDimY );
    RectangleCopy( outBuf, blockDimX, blockDimY,
                   dstarr, dimX, dimY,
                   0, 0,
                   indices[i*2 + 0], indices[i*2 + 1],
                   blockDimX, blockDimY );
    delete[] inBuf;
    delete[] outBuf;
  }
}
