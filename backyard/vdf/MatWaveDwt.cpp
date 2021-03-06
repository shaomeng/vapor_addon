#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vapor/MatWaveDwt.h>
#ifdef WIN32
#include <float.h>
#define isfinite _finite
#endif
#include <sys/time.h>

using namespace VAPoR;

namespace {

#ifdef DEBUG
template <class T>
void printmatrix1d(const char *msg, const T *mat, size_t nx) {


	fprintf(stdout, "%s\n", msg);
		for (size_t i=0; i<nx; i++) {
			fprintf(stdout, "%-10f", (float) mat[i]);
		}
		fprintf(stdout, "\n");
}

template <class T>
void printmatrix2d(const char *msg, const T *mat, size_t nx, size_t ny) {

	fprintf(stdout, "%s\n", msg);
	for (size_t j=0; j<ny; j++) {
		fprintf(stdout, "%-4d", j);
		for (size_t i=0; i<nx; i++) {
			fprintf(stdout, "%-10f", (float) mat[j*nx+i]);
		}
		fprintf(stdout, "\n");
	}
}

template <class T>
void printmatrix3d(const char *msg, const T *mat, size_t nx, size_t ny, size_t nz) {

	fprintf(stdout, "%s\n", msg);
	
	for (size_t k=0; k<nz; k++) {
		fprintf(stdout, "\n");
		for (size_t j=0; j<ny; j++) {
			fprintf(stdout, "%-4d", k);
			for (size_t i=0; i<nx; i++) {
				fprintf(stdout, "%-10f", (float) mat[k*nx*ny + j*nx+i]);
			}
			fprintf(stdout, "\n");
		}
	}
}

#else

#define printmatrix1d(a,b,c)
#define printmatrix2d(a,b,c,d)
#define printmatrix3d(a,b,c,d,e)

#endif


double *buf_alloc(
	double *buf, size_t *buf_size, size_t len
) {
	if (*buf_size < len) {
		if (buf) delete [] buf;
		buf = new double[len];
		*buf_size = len;
	}
	if (buf == NULL) {
		MatWaveDwt::SetErrMsg(
			"Memory allocation of %lu bytes failed", 
			(size_t) len * sizeof(double)
		);
		return(NULL);
	}
	return(buf);
}


/*-------------------------------------------
 * Signal Extending
 *-----------------------------------------*/


template <class T>
int wextend_1D_center (
	const T *sigIn, size_t sigInLen,
	double *sigOut, size_t addLen,
	MatWaveBase::dwtmode_t leftExtMethod,
	MatWaveBase::dwtmode_t rightExtMethod,
	bool invalid_float_abort

) {
  int count = 0;

  for (count = 0; count < addLen; count++)
    {
      sigOut[count] = 0;
      sigOut[count + sigInLen + addLen] = 0;
    }

  for (count = 0; count < sigInLen; count++)
    {
      if (! isfinite((double) sigIn[count])) {
        if (invalid_float_abort) {
          MatWaveDwt::SetErrMsg(
            "Invalid floating point value : %lf",
            (double) sigIn[count]
          );
          return(-1);
        }
        sigOut[count + addLen] = 0.0;
      }
      else {
        sigOut[count + addLen] = sigIn[count];
      }
    }
  if (! addLen) return(0);

  switch (leftExtMethod) {
  case MatWaveBase::ZPD: break;
  case MatWaveBase::SYMH:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count] = sigIn[addLen - count - 1];
	}
      break;
    }
  case MatWaveBase::SYMW:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count] = sigIn[addLen - count];
	}
      break;
    }
  case MatWaveBase::ASYMH:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count] = sigIn[addLen - count - 1] * (-1);
	}
      break;
    }
  case MatWaveBase::ASYMW:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count] = sigIn[addLen - count] * (-1);
	}
      break;
    }
  case MatWaveBase::SP0:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count] = sigIn[0];
	}
      break;
    }
  case MatWaveBase::SP1:
    {
      for (count = (addLen - 1); count >= 0; count--)
	{
		sigOut[count] = sigIn[0]-(sigIn[1]-sigIn[0])*(addLen-count);
	}
      break;
    }
  case MatWaveBase::PPD:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count] = sigIn[sigInLen - addLen + count];
	}
      break;
    }
  case MatWaveBase::PER:
    {
      if (sigInLen%2 == 0)
	{
	  for (count = 0; count < addLen; count++)
	    {
	      sigOut[count] = sigIn[sigInLen - addLen + count];
	    }
	}
      else
	{
	  sigOut[addLen-1] = sigIn[sigInLen-1];
	  addLen--;
	  for (count = 0; count < addLen; count++)
	    {
	      sigOut[count] = sigIn[sigInLen - addLen + count];
	    }
	}
      break;
    }
  default: break;
  }

  switch (rightExtMethod) {
  case MatWaveBase::ZPD: break;
  case MatWaveBase::SYMH:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count + sigInLen + addLen] =
	    sigIn[sigInLen - count - 1];
	}
      break;
    }
  case MatWaveBase::SYMW:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count + sigInLen + addLen] =
	    sigIn[sigInLen - count - 2];
	}
      break;
    }
  case MatWaveBase::ASYMH:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count + sigInLen + addLen] =
	    sigIn[sigInLen - count - 1] * (-1);
	}
      break;
    }
  case MatWaveBase::ASYMW:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count + sigInLen + addLen] =
	    sigIn[sigInLen - count - 2] * (-1);
	}
      break;
    }
  case MatWaveBase::SP0:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count + sigInLen + addLen] = 
	    sigIn[sigInLen - 1];
	}
      break;
    }
  case MatWaveBase::SP1:
    {
      for (count = (addLen - 1); count >= 0; count--)
	{
		sigOut[sigInLen + 2 * addLen - count - 1] = 
			sigIn[sigInLen - 1] - (sigIn[sigInLen-2] - sigIn[sigInLen-1])*(addLen-count);
	}
      break;
    }
  case MatWaveBase::PPD:
    {
      for (count = 0; count < addLen; count++)
	{
	  sigOut[count + sigInLen + addLen] = sigIn[count];
	}
      break;
    }
  case MatWaveBase::PER:
    {
      if (sigInLen%2 == 0)
	{
	  for (count = 0; count < addLen; count++)
	    {
	      sigOut[count + sigInLen + addLen] = sigIn[count];
	    }
	}
      else
	{
	  sigOut[addLen+sigInLen] = sigIn[sigInLen-1];
	  addLen--;
	  for (count = 0; count < addLen; count++)
	    {
	      sigOut[count + sigInLen + addLen+2] = sigIn[count];
	    }
	}
      break;
    }
  default: break;
  }

  return(0);
}

//
// Perform single-level, 1D forward wavelet transform 
// (convolution + downsampling)
//
//
// The number of samples computed for both cA and cD is: sigInLen / 2
//
// If oddlow is true the odd indexed low pass samples are computed (the first
// input sample is ignored), else the even samples are computed. This
// parameter provides control over the centering of the filter. Similar
// for the oddhigh parameter.
//
// sigIn must contain sigInLen + filterLen + 1 samples if oddlow or oddhigh
// is true, otherwise sigInLen + filterLen samples are required
//
// See G. Strang and T. Nguyen, "Wavelets and Filter Banks", chap 8, finite
// length filters
//
void
forward_xform (
	const double *sigIn, size_t sigInLen, 
	const double *low_filter, const double *high_filter, 
	int filterLen, double *cA, double *cD, bool oddlow, bool oddhigh
) {
//	assert(sigInLen > filterLen);

	size_t xlstart = oddlow ? 1 : 0;
	size_t xl;
	size_t xhstart = oddhigh ? 1 : 0;
	size_t xh;

	for (size_t yi = 0; yi < sigInLen; yi += 2) {
		cA[yi>>1] = cD[yi>>1] = 0.0;    // Sam: cA[ yi/2 ] = cD[ yi/2 ] = 0.0;

		xl = xlstart;
		xh = xhstart;

		for (int k = filterLen - 1; k >= 0; k--) {
			cA[yi>>1] += low_filter[k] * sigIn[xl];
			cD[yi>>1] += high_filter[k] * sigIn[xh];
			xl++;
			xh++;
		}
		xlstart+=2;
		xhstart+=2;
	}

	return;
}

void
inverse_xform_even (
	const double *cA, const double *cD, size_t sigInLen, 
	const double *low_filter, const double *high_filter, 
	int filterLen, double *sigOut, bool matlab
) {
	size_t xi; // input and out signal indecies
	int k; // filter index

	assert((filterLen % 2) == 0);

	for (size_t yi = 0; yi < 2*sigInLen; yi++ ) {
		sigOut[yi] = 0.0;

		if (matlab  || (filterLen>>1)%2) { // odd length half filter
			xi = yi >> 1;
			if (yi % 2) {
				k =  filterLen - 1;
			} else {
				k =  filterLen - 2;
			}
		} else {
			xi = (yi+1) >> 1;
			if (yi % 2) {
				k = filterLen - 2;
			} else {
				k = filterLen - 1;
			}
		}

		for (; k >= 0; k-=2) {
			sigOut[yi] += (low_filter[k] * cA[xi]) + (high_filter[k] * cD[xi]);
			xi++;
		}
	}

	return;
}

//
// Inverse transform for odd length, symmetric filters. In this case
// it is assumed that cA coefficients come from even indexed samples
// and cD coefficients come from odd indexed samples.
//
// See G. Strang and T. Nguyen, "Wavelets and Filter Banks", 
// chap 8, finite length filters
//
void
inverse_xform_odd (
	const double *cA, const double *cD, size_t sigInLen, 
	const double *low_filter, const double *high_filter, 
	int filterLen, double *sigOut
) {
	size_t xi; // input and out signal indecies
	int k; // filter index

	assert((filterLen % 2) == 1);

	for (size_t yi = 0; yi < 2*sigInLen; yi++ ) {
		sigOut[yi] = 0.0;

		xi = (yi+1) >> 1;
		if (yi % 2) {
			k = filterLen - 2;
		} else {
			k = filterLen - 1;
		}
		for (; k >= 0; k-=2) {
			sigOut[yi] += (low_filter[k] * cA[xi]);
			xi++;
		}

		xi = (yi) >> 1;
		if (yi % 2) {
			k = filterLen - 1;
		} else {
			k = filterLen - 2;
		}
		for (; k >= 0; k-=2) {
			sigOut[yi] += (high_filter[k] * cD[xi]);
			xi++;
		}

	}

	return;
}


#define Minimum(a,b) ((a<b)?a:b)
#define BlockSize 32

  //
  // blocked submatrix transpose suitable for multithreading
  //   *a : pointer to input matrix
  //   *b : pointer to output matrix
  //    p1,p2: starting index of submatrix (row,col)
  //    m1,m2: size of submatrix (row,col)
  //    s1,s2: size of entire matrix (row,col)
  //
  
template <class T, class U>
void transpose(const T *a,U *b,size_t p1,size_t p2,size_t m1,size_t m2,size_t s1,size_t s2)
  {
    size_t I1,I2;
    size_t i1,i2;
    size_t q,r;
    register double c0;
    const size_t block=BlockSize;
    for(I2=p2;I2<p2+m2;I2+=block)
      for(I1=p1;I1<p1+m1;I1+=block)
	for(i2=I2;i2<Minimum(I2+block,p2+m2);i2++)
	  for(i1=I1;i1<Minimum(I1+block,p1+m1);i1++)
	    {
	      q=i2*s1+i1;
	      r=i1*s2+i2;
	      c0=a[q];
	      b[r]=c0;
	    }
  }

  // specialization for Real -> Complex
  // note the S1 matrix dimension is for the Real matrix
  // and the size of the Complex output is then s2 x (S1/2+1)

  
  //
  // blocked matrix transpose single threaded
  //   *a : pointer to input matrix
  //   *b : pointer to output matrix
  //    s1,s2: size of entire matrix (row,col)
  //
  
template <class T, class U>
void transpose(const T *a,U *b,size_t s1,size_t s2)
  {
    transpose(a,b,0,0,s1,s2,s1,s2);
  }


};


// Sam sets timers
void MatWaveDwt::ResetFXformTimer()
{
    _fxformTimer.tv_sec = 0;
    _fxformTimer.tv_usec = 0;
}
void MatWaveDwt::ResetIXformTimer()
{
    _ixformTimer.tv_sec = 0;
    _ixformTimer.tv_usec = 0;
}
struct timeval MatWaveDwt::GetFXformTimer()
{
    return _fxformTimer;
}
struct timeval MatWaveDwt::GetIXformTimer()
{
    return _ixformTimer;
}
void MatWaveDwt::AddFXformTime( struct timeval &start, struct timeval &end )
{
    long int sec = end.tv_sec - start.tv_sec;
    long int usec = end.tv_usec - start.tv_usec;
    if( usec < 0 ) {
        usec += 1000000;
        sec--;
    }
    _fxformTimer.tv_sec += sec;
    _fxformTimer.tv_usec += usec;
}
void MatWaveDwt::AddIXformTime( struct timeval &start, struct timeval &end )
{
    long int sec = end.tv_sec - start.tv_sec;
    long int usec = end.tv_usec - start.tv_usec;
    if( usec < 0 ) {
        usec += 1000000;
        sec--;
    }
    _ixformTimer.tv_sec += sec;
    _ixformTimer.tv_usec += usec;
}
// Sam test CUDA
void MatWaveDwt::TestCuda() {}


MatWaveDwt::MatWaveDwt(
	const string &wname, const string &mode
) : MatWaveBase(wname, mode) {

	_dwt1dBufSize = 0;
	_dwt1dBuf = NULL;

	_dwt2dBufSize = 0;
	_dwt2dBuf = NULL;

	_dwt3dBuf1Size = 0;
	_dwt3dBuf1 = NULL;
	_dwt3dBuf2Size = 0;
	_dwt3dBuf2 = NULL;
} 

MatWaveDwt::MatWaveDwt(const string &wname ) : MatWaveBase(wname) {

	_dwt1dBufSize = 0;
	_dwt1dBuf = NULL;

	_dwt2dBufSize = 0;
	_dwt2dBuf = NULL;

	_dwt3dBuf1Size = 0;
	_dwt3dBuf1 = NULL;
	_dwt3dBuf2Size = 0;
	_dwt3dBuf2 = NULL;
} 

MatWaveDwt::~MatWaveDwt() {
	if (_dwt1dBuf) delete [] _dwt1dBuf;

	if (_dwt2dBuf) delete [] _dwt2dBuf;

	if (_dwt3dBuf1) delete [] _dwt3dBuf1;
	if (_dwt3dBuf2) delete [] _dwt3dBuf2;
}

template <class T, class U>
int dwt_template(
	MatWaveDwt *dwt,
	const T *sigIn, size_t sigInLen, const WaveFiltBase *wf,
	MatWaveBase::dwtmode_t mode,
	U *cA, U *cD, size_t L[3], double **buf, size_t *bufsize
) {

	if (dwt->wmaxlev(sigInLen) < 1) {
		MatWaveDwt::SetErrMsg("Can't transform signal of length : %d", sigInLen);
		return(-1);
	}
	if (mode == MatWaveBase::PER) {
		MatWaveDwt::SetErrMsg("Invalid boundary extension mode: %d", mode);
		return(-1);
	}

	L[0] = dwt->approxlength(sigInLen);
	L[1] = dwt->detaillength(sigInLen);
	L[2] = sigInLen;


	int filterLen = wf->GetLength();

	//
	// See if we can do symmetric convolution
	//
	bool do_sym_conv = false;
	if (wf->issymmetric()) {
		if (
			(mode == MatWaveBase::SYMW && (filterLen % 2)) ||
			(mode == MatWaveBase::SYMH && (! (filterLen % 2)))
		)  {
		
			do_sym_conv = true;
		}
	}

	//cout << "filter length " << filterLen << endl;
	//printmatrix1d("dwt: low pass decomp filter", wf->GetLowDecomFilCoef(), filterLen);
	//printmatrix1d("dwt: high pass decomp filter", wf->GetHighDecomFilCoef(), filterLen);
	//cout << endl;
	printmatrix1d("dwt: input signal", sigIn,sigInLen);

	// length of signal after boundary extension. We extend both
	// left and right boundary by the width of the filter. 
	//
	size_t sigConvolvedLen =  L[0] + L[1];
	size_t extendLen;

	bool oddlow = true;
	bool oddhigh = true;
	if (filterLen % 2) oddlow = false;
	if (do_sym_conv) {
		extendLen = filterLen>>1;
		if (sigInLen % 2) sigConvolvedLen += 1;
	}
	else {
		extendLen = filterLen-1;
	}
	size_t sigExtendedLen = sigInLen + (2*extendLen);

	*buf = buf_alloc(
		*buf, bufsize, sigExtendedLen + sigConvolvedLen
	);
	if (! buf) return(-1);

	double *sigExtended = *buf;
	double  *sigConvolved = sigExtended + sigExtendedLen;

	// Signal boundary extension
	//
	int rc = wextend_1D_center(
		sigIn, sigInLen, sigExtended, 
		extendLen, mode, mode, dwt->InvalidFloatAbortOnOff()
	);
	if (rc<0) return(-1);
	printmatrix1d("dwt: extended signal", sigExtended, sigExtendedLen);

    // Sam adds timing
    // <==============>
    struct timeval sam_start, sam_end;
    gettimeofday( &sam_start, NULL );
	forward_xform(
		sigExtended, L[0]+L[1], wf->GetLowDecomFilCoef(), 
		wf->GetHighDecomFilCoef(), filterLen, sigConvolved, sigConvolved+L[0],
		oddlow, oddhigh
	);
    gettimeofday( &sam_end, NULL );
    dwt->AddFXformTime( sam_start, sam_end );
    // <==============>

	for (size_t i=0; i<L[0]; i++) {
		cA[i] = sigConvolved[i];
	}
	printmatrix1d("dwt: convolved lowpass signal", cA, L[0]);

	for (size_t i=0; i<L[1]; i++) {
		cD[i] = sigConvolved[i+L[0]];
	}
	printmatrix1d("dwt: convolved high signal", cD, L[1]);

	return(0);

}


int MatWaveDwt::dwt(
	const double *sigIn, size_t sigInLen, double *C, size_t L[3]
) {
	double *cA = C;
	double *cD = C + approxlength(sigInLen);

	return(dwt_template(this,
		sigIn, sigInLen, wavelet(), dwtmodeenum(), cA, cD, L,
		&_dwt1dBuf,  &_dwt1dBufSize
	));
}

int MatWaveDwt::dwt(
	const float *sigIn, size_t sigInLen, float *C, size_t L[3]
) {
	float *cA = C;
	float *cD = C + approxlength(sigInLen);
	return(dwt_template(this,
		sigIn, sigInLen, wavelet(), dwtmodeenum(), cA, cD, L,
		&_dwt1dBuf,  &_dwt1dBufSize
	));
}

int MatWaveDwt::dwt(
	const double *sigIn, size_t sigInLen, double *cA, double *cD, size_t L[3]
) {
	return(dwt_template(this,
		sigIn, sigInLen, wavelet(), dwtmodeenum(), cA, cD, L,
		&_dwt1dBuf,  &_dwt1dBufSize
	));
}

int MatWaveDwt::dwt(
	const float *sigIn, size_t sigInLen, float *cA, float *cD, size_t L[3]
) {
	return(dwt_template(this,
		sigIn, sigInLen, wavelet(), dwtmodeenum(), cA, cD, L,
		&_dwt1dBuf,  &_dwt1dBufSize
	));
}


template <class T, class U>
int idwt_template(
	MatWaveDwt *dwt,
	const T *cA, const T *cD, const size_t L[3], const WaveFiltBase *wf,
	MatWaveBase::dwtmode_t mode, U *sigOut,
	double **buf, size_t *bufsize
) {
	if (mode == MatWaveBase::PER) {
		MatWaveDwt::SetErrMsg("Invalid boundary extension mode: %d", mode);
		return(-1);
	}

	int filterLen = wf->GetLength();

	bool do_sym_conv = false;
	MatWaveBase::dwtmode_t cALeftMode = mode;
	MatWaveBase::dwtmode_t cARightMode = mode;
	MatWaveBase::dwtmode_t cDLeftMode = mode;
	MatWaveBase::dwtmode_t cDRightMode = mode;
	if (wf->issymmetric()) {
		if (
			(mode == MatWaveBase::SYMW && (filterLen % 2)) ||
			(mode == MatWaveBase::SYMH && (! (filterLen % 2)))
		)  {

			if (mode == MatWaveBase::SYMH) {
				cDLeftMode = MatWaveBase::ASYMH;
				if (L[2]%2) {
					cARightMode = MatWaveBase::SYMW;
					cDRightMode = MatWaveBase::ASYMW;
				}
				else {
					cDRightMode = MatWaveBase::ASYMH;
				}
			}
			else {
				cDLeftMode = MatWaveBase::SYMH;
				if (L[2]%2) {
					cARightMode = MatWaveBase::SYMW;
					cDRightMode = MatWaveBase::SYMH;
				}
				else {
					cARightMode = MatWaveBase::SYMH;
				}
			}
		
			do_sym_conv = true;
		}
	}

	size_t cATempLen, cDTempLen, reconTempLen;

	size_t extendLen = 0;
	size_t cDPadLen = 0;
	if (do_sym_conv) {
		extendLen = filterLen>>2;
		if ((L[0] > L[1]) && (mode == MatWaveBase::SYMH)) cDPadLen = L[0];

		cATempLen = L[0] + (2*extendLen);

		// cD length must be same as cA (it's not for odd length signals)
		//
		cDTempLen = cATempLen;
	} else {
		cATempLen = L[0];
		cDTempLen = L[1];
	}
	reconTempLen = L[2];
	if (reconTempLen % 2) reconTempLen++;

	*buf = buf_alloc(
		*buf, bufsize, 
		cATempLen + cDTempLen + reconTempLen + cDPadLen
	);
	if (! buf) return(-1);


	double *cATemp = *buf;
	double *cDTemp = cATemp + cATempLen;
	double *reconTemp = cDTemp + cDTempLen;
	double *cDPad = reconTemp + reconTempLen;

	//printmatrix1d("idwt: low pass reconstruct filter", wf->GetLowReconFilCoef(), filterLen);
	//printmatrix1d("idwt: high pass reconstruct filter", wf->GetHighReconFilCoef(), filterLen);
	//cout << endl;

	// For symmetric filters we need to add the boundary coefficients 
	//
	if (do_sym_conv) {
		int rc = wextend_1D_center(
			cA, L[0], cATemp, 
			extendLen, cALeftMode, cARightMode, dwt->InvalidFloatAbortOnOff()
		);
		if (rc<0) return(-1);

		// For odd length signals we need to add back the missing final cD 
		// coefficient before signal extension.
		// See G. Strang and T. Nguyen, "Wavelets and Filter Banks", 
		// chap 8, finite length filters
		// 
		if (cDPadLen) {
			for (size_t i=0; i<L[1]; i++) cDPad[i] = cD[i];
			cDPad[L[1]] = 0.0;

			rc = wextend_1D_center(
				cDPad, L[0], cDTemp, extendLen, cDLeftMode, cDRightMode, 
				dwt->InvalidFloatAbortOnOff()
			);
			if (rc<0) return(-1);
		}
		else {
			rc = wextend_1D_center(
				cD, L[1], cDTemp, extendLen, cDLeftMode, cDRightMode, 
				dwt->InvalidFloatAbortOnOff()
			);
			if (rc<0) return(-1);
		}
	}
	else {
		for (size_t i=0; i<L[0]; i++) {
      		if (! isfinite((double) cA[i])) {
				if (dwt->InvalidFloatAbortOnOff()) {
					MatWaveDwt::SetErrMsg(
						"Invalid floating point value : %lf", (double) cA[i]
					);
					return(-1);
				}
				cATemp[i] = 0.0;
			}
			else {
				cATemp[i] = (double) cA[i];
			}
		}

		for (size_t i=0; i<L[1]; i++) {
      		if (! isfinite((double) cD[i])) {
				if (dwt->InvalidFloatAbortOnOff()) {
					MatWaveDwt::SetErrMsg(
						"Invalid floating point value : %lf", (double) cD[i]
					);
					return(-1);
				}
				cDTemp[i] = 0.0;
			}
			else {
				cDTemp[i] = (double) cD[i];
			}
		}
		
	}

	printmatrix1d("idwt: extended cA signal", cATemp, cATempLen);
	printmatrix1d("idwt: extended cD signal", cDTemp, cDTempLen);

    // Sam adds timing
    // <==============>
    struct timeval sam_start, sam_end;
    gettimeofday( &sam_start, NULL );
	if (filterLen % 2) {
		
		inverse_xform_odd (
			cATemp, cDTemp, L[0], wf->GetLowReconFilCoef(), 
			wf->GetHighReconFilCoef(), filterLen,
			reconTemp 
		);
	}
	else {
		inverse_xform_even (
			cATemp, cDTemp, L[0], wf->GetLowReconFilCoef(), 
			wf->GetHighReconFilCoef(), filterLen,
			reconTemp, ! do_sym_conv
		);
	}
    gettimeofday( &sam_end, NULL );
    dwt->AddFXformTime( sam_start, sam_end );
    // <==============>

	for(size_t count=0;count<L[2];count++) {
		sigOut[count] = (U) reconTemp[count];
	}
	printmatrix1d("idwt: reconstructed signal", sigOut, L[2]);

	return(0);
}

int MatWaveDwt::idwt(
	const double *C, const size_t L[3], double *sigOut
) {
	const double *cA = C;
	const double *cD = C + L[0];
	return idwt_template(
		this, cA, cD, L, wavelet(), dwtmodeenum(), sigOut, 
		&_dwt1dBuf,  &_dwt1dBufSize
	);
} 

int MatWaveDwt::idwt(
	const float *C, const size_t L[3], float *sigOut
) {
	const float *cA = C;
	const float *cD = C + L[0];
	return idwt_template(
		this, cA, cD, L, wavelet(), dwtmodeenum(), sigOut, 
		&_dwt1dBuf,  &_dwt1dBufSize
	);
} 

int MatWaveDwt::idwt(
	const double *cA, const double *cD, const size_t L[3], double *sigOut
) {
	return idwt_template(
		this, cA, cD, L, wavelet(), dwtmodeenum(), sigOut, 
		&_dwt1dBuf,  &_dwt1dBufSize
	);
} 

int MatWaveDwt::idwt(
	const float *cA, const float *cD, const size_t L[3], float *sigOut
) {
	return idwt_template(
		this, cA, cD, L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt1dBuf,  &_dwt1dBufSize
	);
} 

template <class T, class U>
int dwt2d_template(
	MatWaveDwt *dwt,
	const T *sigIn, size_t sigInX, size_t sigInY, const WaveFiltBase *wf,
    MatWaveBase::dwtmode_t mode, U *cA,  U *cDh, U *cDv, U *cDd, size_t L[10],
	double **buf2d, size_t *buf2dsize, double **buf1d, size_t *buf1dsize
) {

	//
	// Store dimensions of cA, cDh, cDv, cDd, and original signal
	// in book keeping vector, L
	//
	// N.B.
	//	L[0] == L[2]
	//	L[1] == L[5]
	//	L[3] == L[7]
	//	L[4] == L[6]
	//
	//
	//      ____L[0]_______L[4]____
	//      |          |          |
	// L[1] |  cA      |  cDv     | L[5]
	//      |  (LL)    |  (HL)    |
	//      |          |          |
	//      |---------------------|
	//      |          |          |
	//      |  cDh     |  cDd     | L[7]
	// L[3] |  (LH)    |  (HH)    |
	//      |          |          |
	//      |__________|__________|
	//         L[2]       L[6]


	L[0] = dwt->approxlength(sigInX); L[1] = dwt->approxlength(sigInY); // cA
	L[2] = dwt->approxlength(sigInX); L[3] = dwt->detaillength(sigInY); // cDh
	L[4] = dwt->detaillength(sigInX); L[5] = dwt->approxlength(sigInY); // cDv
	L[6] = dwt->detaillength(sigInX); L[7] = dwt->detaillength(sigInY); // cDd
	L[8] = sigInX; L[9] = sigInY;

	// First: transform rows
	//
	size_t passXLen = (L[0] + L[4]) * sigInY;
	size_t transposeLen = max(L[0],L[4]) * sigInY;
	size_t passYLen = max(L[0],L[4]) * (L[1] + L[3]);
	
	*buf2d = buf_alloc(*buf2d, buf2dsize, passXLen + transposeLen + passYLen);
	if (! buf2d) return(-1);

	double *cAXbuf = *buf2d;
	double *cDXbuf = cAXbuf + (L[0] * sigInY);

	double *buftranspose = cAXbuf + passXLen;

	double *cAYbuf = buftranspose + transposeLen;
	double *cDYbuf = cAYbuf + (max(L[0],L[4]) * L[1]);

	int rc;
	for (size_t y = 0; y<sigInY; y++) {
		size_t xL[3];
		const T *row = &sigIn[sigInX*y];
		double *cAptr = &cAXbuf[L[0]*y];
		double *cDptr = &cDXbuf[L[4]*y];

		rc = dwt_template(
			dwt, row, sigInX, wf, mode, cAptr, cDptr, xL,
			buf1d, buf1dsize
		);
		if (rc < 0) return(-1);
	}

	// Second: transform columns. First approximation coefficients, then
	// detail coefficients
	//

	transpose(cAXbuf, buftranspose, L[0], sigInY);

	for (size_t y = 0; y<L[0]; y++) {
		size_t yL[3];
		const double *row = &buftranspose[sigInY*y];
		double *cAptr = &cAYbuf[L[1]*y];
		double *cDptr = &cDYbuf[L[3]*y];

		rc = dwt_template(
			dwt, row, sigInY, wf, mode, cAptr, cDptr, yL,
			buf1d, buf1dsize
		);
		if (rc < 0) return(-1);
	}

	transpose(cAYbuf, cA, L[1], L[0]);
	transpose(cDYbuf, cDh, L[3], L[2]);

printmatrix2d("cA", cA, L[0], L[1]);
printmatrix2d("cDh", cDh, L[2], L[3]);

	// Now detail coefficients
	//
	//
	transpose(cDXbuf, buftranspose, L[4], sigInY);

	for (size_t y = 0; y<L[4]; y++) {
		size_t yL[3];
		const double *row = &buftranspose[sigInY*y];
		double *cAptr = &cAYbuf[L[1]*y];
		double *cDptr = &cDYbuf[L[3]*y];

		rc = dwt_template(
			dwt, row, sigInY, wf, mode, cAptr, cDptr, yL,
			buf1d, buf1dsize
		);
		if (rc < 0) return(-1);
	}
	transpose(cAYbuf, cDv, L[5], L[4]);
	transpose(cDYbuf, cDd, L[7], L[6]);

printmatrix2d("cDv", cDv, L[4], L[5]);
printmatrix2d("cDd", cDd, L[6], L[7]);

	return(0);
}


int MatWaveDwt::dwt2d(
	const double *sigIn, size_t sigInX, size_t sigInY, double *C, size_t L[10]
) {
	double *cA = C;
	double *cDh = cA + (approxlength(sigInX) * approxlength(sigInY));
	double *cDv = cDh + (approxlength(sigInX) * detaillength(sigInY));
	double *cDd = cDv + (detaillength(sigInX) * approxlength(sigInY));

	return dwt2d_template(
		this, sigIn, sigInX, sigInY, wavelet(), dwtmodeenum(),
		cA, cDh, cDv, cDd, L,
		&_dwt2dBuf,  &_dwt2dBufSize, 
		&_dwt1dBuf,  &_dwt1dBufSize
	);
}

int MatWaveDwt::dwt2d(
	const float *sigIn, size_t sigInX, size_t sigInY, float *C, size_t L[10]
) {
	float *cA = C;
	float *cDh = cA + (approxlength(sigInX) * approxlength(sigInY));
	float *cDv = cDh + (approxlength(sigInX) * detaillength(sigInY));
	float *cDd = cDv + (detaillength(sigInX) * approxlength(sigInY));

	return dwt2d_template(
		this, sigIn, sigInX, sigInY, wavelet(), dwtmodeenum(),
		cA, cDh, cDv, cDd, L,
		&_dwt2dBuf,  &_dwt2dBufSize, &_dwt1dBuf,  &_dwt1dBufSize
	);
}

int MatWaveDwt::dwt2d(
	const double *sigIn, size_t sigInX, size_t sigInY, 
	double *cA, double *cDh, double *cDv, double *cDd, size_t L[10]
) {
	return dwt2d_template(
		this, sigIn, sigInX, sigInY, wavelet(), dwtmodeenum(),
		cA, cDh, cDv, cDd, L,
		&_dwt2dBuf,  &_dwt2dBufSize, &_dwt1dBuf,  &_dwt1dBufSize
	);
}

int MatWaveDwt::dwt2d(
	const float *sigIn, size_t sigInX, size_t sigInY, 
	float *cA, float *cDh, float *cDv, float *cDd, size_t L[10]
) {
	return dwt2d_template(
		this, sigIn, sigInX, sigInY, wavelet(), dwtmodeenum(),
		cA, cDh, cDv, cDd, L,
		&_dwt2dBuf,  &_dwt2dBufSize, &_dwt1dBuf,  &_dwt1dBufSize
	);
}


template <class T, class U>
int idwt2d_template(
	MatWaveDwt *dwt,
	const T *cA, const T *cDh, const T *cDv, const T *cDd,
	const size_t L[10], const WaveFiltBase *wf,
	MatWaveBase::dwtmode_t mode, U *sigOut, 
	double **buf2d, size_t *buf2dsize,
	double **buf1d, size_t *buf1dsize
) {

    size_t passYLen = max(L[0],L[4]) * (L[1] + L[3]);
    size_t transposeLen = max(L[0],L[4]) * L[9];
    size_t passXLen = (L[0] + L[4]) * L[9];

	*buf2d = buf_alloc(
		*buf2d, buf2dsize, passYLen + transposeLen + passXLen
	);
	if (! buf2d) return(-1);

	double *cAYbuf = *buf2d;
	double *cDYbuf = cAYbuf + (max(L[0],L[4]) * L[1]);

	double *buftranspose = cAYbuf + passYLen;

    double *cAXbuf = buftranspose + transposeLen;
    double *cDXbuf = cAXbuf + (L[0] * L[9]);

	// First: transform columns. First detail coefficients, then
	// approximation coefficients
	//

	// cDv and cDd detail coefficients
	//

	transpose(cDv, cAYbuf, L[4], L[1]);
	transpose(cDd, cDYbuf, L[4], L[3]);
	int rc; 
	for (size_t y = 0; y<L[4]; y++) {
		size_t yL[3] = {L[1], L[3], L[9]};
		const double *cAptr = &cAYbuf[L[1]*y];
		const double *cDptr = &cDYbuf[L[3]*y];
		double *row = &buftranspose[L[9]*y];

		rc = idwt_template(
			dwt, cAptr, cDptr, yL, wf, mode, row,
			buf1d, buf1dsize
		);
		if (rc < 0) return (-1);
	}
	transpose(buftranspose, cDXbuf, L[9], L[4]);
	//printmatrix2d("cDXbuf", cDXbuf, L[4], L[9]);


	// cA approximation and cDh detail coefficients
	//

	transpose(cA, cAYbuf, L[0], L[1]);
	transpose(cDh, cDYbuf, L[0], L[3]);
	for (size_t y = 0; y<L[0]; y++) {
		size_t yL[3] = {L[1], L[3], L[9]};
		const double *cAptr = &cAYbuf[L[1]*y];
		const double *cDptr = &cDYbuf[L[3]*y];
		double *row = &buftranspose[L[9]*y];

		rc = idwt_template(
			dwt, cAptr, cDptr, yL, wf, mode, row,
			buf1d, buf1dsize
		);
		if (rc < 0) return (-1);
	}
	transpose(buftranspose, cAXbuf, L[9], L[0]);

	//
	//  Second: tranform rows
	//
	for (size_t y = 0; y<L[9]; y++) {
		size_t xL[3] = {L[0], L[4], L[8]};
		const double *cAptr = &cAXbuf[L[0]*y];
		const double *cDptr = &cDXbuf[L[4]*y];
		U *row = &sigOut[L[8]*y];

		rc = idwt_template(
			dwt, cAptr, cDptr, xL, wf, mode, row,
			buf1d, buf1dsize
		);
		if (rc < 0) return (-1);
	}


	return(0);
}

int MatWaveDwt::idwt2d(
	const double *C, const size_t L[10], double *sigOut
) {
	const double *cA = C;
	const double *cDh = cA + (L[0] * L[1]);
	const double *cDv = cDh + (L[2] * L[3]);
	const double *cDd = cDv + (L[4] * L[5]);

	return idwt2d_template(
		this, cA, cDh, cDv, cDd, L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
}


int MatWaveDwt::idwt2d(
	const float *C, const size_t L[10], float *sigOut
) {
	const float *cA = C;
	const float *cDh = cA + (L[0] * L[1]);
	const float *cDv = cDh + (L[2] * L[3]);
	const float *cDd = cDv + (L[4] * L[5]);

	return idwt2d_template(
		this, cA, cDh, cDv, cDd, L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
}

int MatWaveDwt::idwt2d(
	const double *cA, const double *cDh, const double *cDv, const double *cDd,
	const size_t L[10], double *sigOut
) {
	return idwt2d_template(
		this, cA, cDh, cDv, cDd, L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
}

int MatWaveDwt::idwt2d(
	const float *cA, const float *cDh, const float *cDv, const float *cDd,
	const size_t L[10], float *sigOut
) {
	return idwt2d_template(
		this, cA, cDh, cDv, cDd, L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
}

template <class T>
int _dwtz_template(
	MatWaveDwt *dwt,
	double *sigIn, size_t sigInX, size_t sigInY, size_t sigInZ,
	const WaveFiltBase *wf, 
	MatWaveBase::dwtmode_t mode,
	T *cA, T *cD, 
	double **buf3d, size_t *buf3dsize,
	double **buf1d, size_t *buf1dsize
) {
	size_t cALen = dwt->approxlength(sigInZ); 
	size_t cDLen = dwt->detaillength(sigInZ);

	size_t sigInLen = sigInX * sigInY * sigInZ;
	size_t sigOutLen = sigInX * sigInY * (cALen + cDLen);

	*buf3d = buf_alloc(*buf3d, buf3dsize, sigInLen + sigOutLen);
	if (! buf3d) return(-1);

	double *sigtranspose = *buf3d;
	double *cAbuf = sigtranspose + sigInLen;
	double *cDbuf = cAbuf + (sigInX * sigInY * cALen);

	// XZ plane transpose
	//
	transpose(sigIn, sigtranspose, sigInX*sigInY, sigInZ);

	int rc; 
	for (size_t z = 0; z<sigInX; z++) {
		for (size_t y = 0; y<sigInY; y++) {
			size_t L[3];
			const double *row = &sigtranspose[sigInY*sigInZ*z + y*sigInZ];
			double *cAptr = &cAbuf[sigInY*cALen*z + y*cALen];
			double *cDptr = &cDbuf[sigInY*cDLen*z + y*cDLen];


			rc = dwt_template(
				dwt, row, sigInZ, wf, mode, cAptr, cDptr, 
				L, buf1d, buf1dsize
			); 
			if (rc < 0) return (-1);
		}
	}


	//
	// now inverse transpose to restore original memory order
	//
	transpose(cAbuf, cA, cALen, sigInX * sigInY);
	transpose(cDbuf, cD, cDLen, sigInX * sigInY);
	return(0);
}

template <class T>
int dwt3d_template(
	MatWaveDwt *dwt,
	const T *sigIn, size_t sigInX, size_t sigInY, size_t sigInZ,
	const WaveFiltBase *wf, MatWaveBase::dwtmode_t mode,
	T *C, size_t L[27],
	double **buf3d1, size_t *buf3d1size,
	double **buf3d2, size_t *buf3d2size,
	double **buf2d, size_t *buf2dsize,
	double **buf1d, size_t *buf1dsize
) {

	//
	// Store dimensions of  LLL, LLH, LHL, LHH, HLL,
    // HLH, HHL, HHH
	// in book keeping vector, L
	//

	//                  L[3]        L[15]
	//               -----------------------
	//              /          /          /|
	//        L[5] /          /          / |
	//            /  LLH     /  HLH     /  |
	//           /          /          /   | L[16]
	//          -----------------------    |
	//         /          /          /|    |
	//   L[2] /          /          / |   /|
	//       /          /          /  |  / |
	//      /___L[0]___/___L[12]__/   | /  | L[22]
	//      |          |          |   |/   | 
	// L[1] |          |          |   /HHH /
	//      |   LLL    |   HLL    |  /|   /
	//      |          |          | / |  / L[23]
	//      |---------------------|/  | /
	//      |          |          |   |/
	//      |          |          |   /
	// L[7] |   LHL    |   HHL    |  /
	//      |          |          | / L[20]
	//      |__________|__________|/
	//          L[6]       L[18]

	// LLL
	//
	L[0] = dwt->approxlength(sigInX); 
    L[1] = dwt->approxlength(sigInY); 
    L[2] = dwt->approxlength(sigInZ); 

	// LLH
	//
	L[3] = dwt->approxlength(sigInX); 
    L[4] = dwt->approxlength(sigInY); 
    L[5] = dwt->detaillength(sigInZ); 

	// LHL
	//
	L[6] = dwt->approxlength(sigInX); 
    L[7] = dwt->detaillength(sigInY); 
    L[8] = dwt->approxlength(sigInZ); 

	// LHH
	//
	L[9] = dwt->approxlength(sigInX); 
    L[10] = dwt->detaillength(sigInY); 
    L[11] = dwt->detaillength(sigInZ); 

	// HLL
	//
	L[12] = dwt->detaillength(sigInX); 
    L[13] = dwt->approxlength(sigInY); 
    L[14] = dwt->approxlength(sigInZ); 

	// HLH
	//
	L[15] = dwt->detaillength(sigInX); 
    L[16] = dwt->approxlength(sigInY); 
    L[17] = dwt->detaillength(sigInZ); 

	// HHL
	//
	L[18] = dwt->detaillength(sigInX); 
    L[19] = dwt->detaillength(sigInY); 
    L[20] = dwt->approxlength(sigInZ); 

	// HHH
	//
	L[21] = dwt->detaillength(sigInX); 
    L[22] = dwt->detaillength(sigInY); 
    L[23] = dwt->detaillength(sigInZ); 

	L[24] = sigInX; L[25] = sigInY; L[26] = sigInZ;


    T *cLLL = C;
    T *cLLH = cLLL + L[0]*L[1]*L[2];
    T *cLHL = cLLH + L[3]*L[4]*L[5];
    T *cLHH = cLHL + L[6]*L[7]*L[8];
    T *cHLL = cLHH + L[9]*L[10]*L[11];
    T *cHLH = cHLL + L[12]*L[13]*L[14];
    T *cHHL = cHLH + L[15]*L[16]*L[17];
    T *cHHH = cHHL + L[18]*L[19]*L[20];

	// First: transform XY planes
	//
	size_t passXYLen = (L[0]+L[12]) * (L[1]+L[7])  * sigInZ;
	
	*buf3d1 = buf_alloc(*buf3d1, buf3d1size, passXYLen);
	if (! buf3d1) return(-1);

	double *xyC = *buf3d1;
	double *cAXYbuf = xyC;
	double *cDhXYbuf = cAXYbuf + (L[0]*L[1]*sigInZ);
	double *cDvXYbuf = cDhXYbuf + (L[6]*L[7]*sigInZ);
	double *cDdXYbuf = cDvXYbuf + (L[12]*L[13]*sigInZ);

	int rc;
	for (size_t z = 0; z<sigInZ; z++) {
		size_t xyL[10];
		const T *plane = &sigIn[sigInX*sigInY*z];
		double *cAptr = &cAXYbuf[L[0]*L[1]*z];
		double *cDhptr = &cDhXYbuf[L[6]*L[7]*z];
		double *cDvptr = &cDvXYbuf[L[12]*L[13]*z];
		double *cDdptr = &cDdXYbuf[L[18]*L[19]*z];

		rc = dwt2d_template(
			dwt, plane, sigInX, sigInY, wf, mode, 
			cAptr, cDhptr, cDvptr, cDdptr, xyL,
			buf2d, buf2dsize, buf1d, buf1dsize
		);
		if (rc < 0) return(rc);
	}

	// Second: transform along Z
	//

	rc = _dwtz_template(
		dwt, cAXYbuf, L[0], L[1], sigInZ, wf, mode, cLLL, cLLH, 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(rc);
	printmatrix3d("cLLL", cLLL, L[0], L[1], L[2]);
	printmatrix3d("cLLH", cLLH, L[3], L[4], L[5]);
	rc = _dwtz_template(
		dwt, cDhXYbuf, L[6], L[7], sigInZ, wf, mode, cLHL, cLHH, 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(rc);
	printmatrix3d("cLHL", cLHL, L[6], L[7], L[8]);
	printmatrix3d("cLHH", cLHH, L[9], L[10], L[11]);
	rc = _dwtz_template(
		dwt, cDvXYbuf, L[12], L[13], sigInZ, wf, mode, cHLL, cHLH, 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(rc);
	printmatrix3d("cHLL", cHLL, L[12], L[13], L[14]);
	printmatrix3d("cHLH", cHLH, L[15], L[16], L[17]);
	rc = _dwtz_template(
		dwt, cDdXYbuf, L[18], L[19], sigInZ, wf, mode, cHHL, cHHH, 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(rc);
	printmatrix3d("cHHL", cHHL, L[18], L[19], L[20]);
	printmatrix3d("cHHH", cHHH, L[21], L[22], L[23]);

	return(0);
}


int MatWaveDwt::dwt3d(
	const double *sigIn, size_t sigInX, size_t sigInY, size_t sigInZ,
	double *C, size_t L[27]
) {

	return dwt3d_template(
		this, sigIn, sigInX, sigInY, sigInZ, wavelet(), dwtmodeenum(), C, L,
		&_dwt3dBuf1, &_dwt3dBuf1Size, &_dwt3dBuf2, &_dwt3dBuf2Size,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
}

int MatWaveDwt::dwt3d(
	const float *sigIn, size_t sigInX, size_t sigInY, size_t sigInZ,
	float *C, size_t L[27]
) {
	return dwt3d_template(
		this, sigIn, sigInX, sigInY, sigInZ, wavelet(), dwtmodeenum(), C, L,
		&_dwt3dBuf1, &_dwt3dBuf1Size, &_dwt3dBuf2, &_dwt3dBuf2Size,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
}

template <class T>
int _idwtz_template(
	MatWaveDwt *dwt,
	T *cA, T *cD, size_t sigInX, size_t sigInY, size_t cALen, size_t cDLen,
	const WaveFiltBase *wf, MatWaveBase::dwtmode_t mode, double *sigOut,
	size_t sigOutZ,
	double **buf3d, size_t *buf3dsize,
	double **buf1d, size_t *buf1dsize
) {


//	size_t sigInLen = 2 * sigInX * sigInY * (cALen + cDLen);
	size_t sigInLen = sigInX * sigInY * (cALen + cDLen);
	size_t sigOutLen = sigInX * sigInY * sigOutZ;
	*buf3d = buf_alloc(
		*buf3d, buf3dsize, sigInLen + sigOutLen
	);
	if (! buf3d) return(-1);

	double *cAtranspose = *buf3d;
	double *cDtranspose = *buf3d + (sigInX * sigInY * cALen);
	double *sigtranspose = cAtranspose + sigInLen;

	//
	// transpose to Z memory order
	//
	transpose(cA, cAtranspose, sigInX*sigInY, cALen);
	transpose(cD, cDtranspose, sigInX*sigInY, cDLen);

	int rc;
	for (size_t z = 0; z<sigInX; z++) {
		for (size_t y = 0; y<sigInY; y++) {
			size_t zL[3] = {cALen, cDLen, sigOutZ};
			const double *cAptr = &cAtranspose[sigInY*cALen*z + y*cALen];
			const double *cDptr = &cDtranspose[sigInY*cDLen*z + y*cDLen];
			double *row = &sigtranspose[sigInY*sigOutZ*z + y*sigOutZ];

			rc = idwt_template(
				dwt, cAptr, cDptr, zL, wf, mode, row,
				buf1d, buf1dsize
			); 
			if (rc < 0) return(-1);
		}
	}

	// Inverse transpose back to X memory order
	//
	transpose(sigtranspose, sigOut, sigOutZ, sigInY*sigInX);

	return(0);
}

template <class T, class U>
int idwt3d_template(
	MatWaveDwt *dwt,
	const T *cLLL, const T *cLLH, const T *cLHL, const T *cLHH,
	const T *cHLL, const T *cHLH, const T *cHHL, const T *cHHH,
	const size_t L[27], const WaveFiltBase *wf,
	MatWaveBase::dwtmode_t mode, U *sigOut, 
	double **buf3d1, size_t *buf3d1size,
	double **buf3d2, size_t *buf3d2size,
	double **buf2d, size_t *buf2dsize,
	double **buf1d, size_t *buf1dsize
) {

	size_t passXYLen = (L[0]+L[12]) * (L[1]+L[7])  * L[26];

	*buf3d1 = buf_alloc(*buf3d1, buf3d1size, passXYLen);
	if (! buf3d1) return(-1);

	double *xyC = *buf3d1;
	double *cAXYbuf = xyC;
	double *cDhXYbuf = cAXYbuf + (L[0]*L[1]*L[26]);
	double *cDvXYbuf = cDhXYbuf + (L[6]*L[7]*L[26]);
	double *cDdXYbuf = cDvXYbuf + (L[12]*L[13]*L[26]);


	// First: inverse transform along Z
	//
	int rc;
	rc = _idwtz_template(
		dwt, cLLL, cLLH, L[0], L[1], L[2], L[5], wf, mode,
		cAXYbuf, L[26], 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(-1);
	rc = _idwtz_template(
		dwt, cLHL, cLHH, L[6], L[7], L[8], L[11], wf, mode,
		cDhXYbuf, L[26], 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(-1);
	rc = _idwtz_template(
		dwt, cHLL, cHLH, L[12], L[13], L[14], L[17], wf, mode,
		cDvXYbuf, L[26], 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(-1);
	rc = _idwtz_template(
		dwt, cHHL, cHHH, L[18], L[19], L[20], L[23], wf, mode,
		cDdXYbuf, L[26], 
		buf3d2, buf3d2size, buf1d, buf1dsize
	);
	if (rc < 0) return(-1);

	// Second: inverse transform XY planes
	//


	size_t xyL[10] = {L[0],L[1],L[6],L[7],L[12],L[13],L[18],L[19],L[24],L[25]};
	for (size_t z = 0; z<L[26]; z++) {

		const double *cAptr = &cAXYbuf[L[0]*L[1]*z];
		const double *cDhptr = &cDhXYbuf[L[6]*L[7]*z];
		const double *cDvptr = &cDvXYbuf[L[12]*L[13]*z];
		const double *cDdptr = &cDdXYbuf[L[18]*L[19]*z];
		U *plane = &sigOut[L[24]*L[25]*z];

		rc = idwt2d_template(
			dwt, cAptr, cDhptr, cDvptr, cDdptr, xyL, wf, mode,
			plane, buf2d, buf2dsize, buf1d, buf1dsize
		); 
		if (rc < 0) return(-1);
	}
	return(0);

}

int MatWaveDwt::idwt3d(
	const double *C, const size_t L[27], double *sigOut
) {
	const double *cLLL = C;
	const double *cLLH = cLLL + L[0]*L[1]*L[2];
	const double *cLHL = cLLH + L[3]*L[4]*L[5];
	const double *cLHH = cLHL + L[6]*L[7]*L[8];
	const double *cHLL = cLHH + L[9]*L[10]*L[11];
	const double *cHLH = cHLL + L[12]*L[13]*L[14];
	const double *cHHL = cHLH + L[15]*L[16]*L[17];
	const double *cHHH = cHHL + L[18]*L[19]*L[20];

	return idwt3d_template(
		this, cLLL, cLLH, cLHL, cLHH, cHLL, cHLH, cHHL, cHHH,
		L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt3dBuf1, &_dwt3dBuf1Size, &_dwt3dBuf2, &_dwt3dBuf2Size,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
} 

int MatWaveDwt::idwt3d(
	const float *C, const size_t L[27], float *sigOut
) {

	const float *cLLL = C;
	const float *cLLH = cLLL + L[0]*L[1]*L[2];
	const float *cLHL = cLLH + L[3]*L[4]*L[5];
	const float *cLHH = cLHL + L[6]*L[7]*L[8];
	const float *cHLL = cLHH + L[9]*L[10]*L[11];
	const float *cHLH = cHLL + L[12]*L[13]*L[14];
	const float *cHHL = cHLH + L[15]*L[16]*L[17];
	const float *cHHH = cHHL + L[18]*L[19]*L[20];

	return idwt3d_template(
		this, cLLL, cLLH, cLHL, cLHH, cHLL, cHLH, cHHL, cHHH,
		L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt3dBuf1, &_dwt3dBuf1Size, &_dwt3dBuf2, &_dwt3dBuf2Size,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
} 

int MatWaveDwt::idwt3d(
	const double *cLLL, const double *cLLH, const double *cLHL,
	const double *cLHH,
	const double *cHLL, const double *cHLH, const double *cHHL,
	const double *cHHH,
	const size_t L[27], double *sigOut
) {

	return idwt3d_template(
		this, cLLL, cLLH, cLHL, cLHH, cHLL, cHLH, cHHL, cHHH,
		L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt3dBuf1, &_dwt3dBuf1Size, &_dwt3dBuf2, &_dwt3dBuf2Size,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
} 

int MatWaveDwt::idwt3d(
	const float *cLLL, const float *cLLH, const float *cLHL, const float *cLHH,
	const float *cHLL, const float *cHLH, const float *cHHL, const float *cHHH,
	const size_t L[27], float *sigOut
) {

	return idwt3d_template(
		this, cLLL, cLLH, cLHL, cLHH, cHLL, cHLH, cHHL, cHHH,
		L, wavelet(), dwtmodeenum(), sigOut,
		&_dwt3dBuf1, &_dwt3dBuf1Size, &_dwt3dBuf2, &_dwt3dBuf2Size,
		&_dwt2dBuf, &_dwt2dBufSize, &_dwt1dBuf, &_dwt1dBufSize
	);
} 
