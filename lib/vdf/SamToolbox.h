//
// A class contains simplified functions for Sam 
// to research on wavelet compression applied in time domain. 
// 

#ifndef _SamToolbox_h_
#define _SamToolbox_h_

#include <iostream>
#include <cstdio>
#include "vapor/Compressor.h"

namespace VAPoR{

struct SamErr{
    double rms = 0.0;
    double max = 0.0;
};

class SamToolbox { 
public:
    // Constructor.
    // Doesn't do anything yet.
//    SamToolbox();
//    ~SamToolbox();

    // Setup right parameters for a constructor.
    int SetupCompressor(Compressor* &_compressor, // input
        const vector<size_t> &bdims, // input
        const vector<size_t> &_cratios, // input
        vector< size_t > &_ncoeffs, // output: num wave coeff. at each compression level
        size_t &_cvectorsize, // output: amount of space need for wavelet coefficients
        vector< SignificanceMap > &sigmaps // output: a signficance map for each compression level
    );

    // Setup coefficient number for each level, based on the compression ratios
    int SetupNCoeffs( const Compressor* compressor,
                      const vector< size_t > &cratios,
                      vector< size_t > &ncoeffs,
                      size_t &cvectorsize );

    // Compare two arrays with stats on errors.
    // 'print' controls if it prints the results out
    // It always returns the RMS value.
    SamErr CompareArrays( const float* arr1, const float* arr2, size_t len, bool print );
    SamErr CompareArrays( const double* arr1, const double* arr2, size_t len, bool print );

    // Read a file of size: len
    // It requires the file name, length, and an allocated space to put read file
    int ReadBlock( std::string filename, size_t len, float* ptr );

    // Output a certain range of an array
    void PrintArray( const float* arr, size_t start, size_t end );

    // A special version of histogram that basically provides the idea of 
    // how many values in an array are non-zero.
    // Calculate the histogram of n of arrays.
    // Input arrays must have the same length: arraylen.
    // The threshold value epsilon defines how small we start to omit values (consider as zero)
    // Output hisgram shows on each position of the array (len in total),
    // how many have n non-zero values, (n-1) non-zero values, (n-2) non-zero values, etc
    void CalcHistogram1( const vector< float* > &arrays,   // Input
                        size_t arraylen,            // Input
                        float epsilon,              // Input
                        vector< size_t > &histogram ); // Output

    float FindMax( const float* arr, size_t len );
    double FindMax( const SamErr* arr, size_t len );

    float Findnth( const float* arr, size_t len, size_t n);

    float CalcRMS( const vector< float > &arr );
    float CalcRMS( const float* arr, size_t len );
    double CalcRMS( const SamErr* arr, size_t len );

};

}

#endif
