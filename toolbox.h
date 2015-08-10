/*
 * This class, Toolbox, is forked from the SamToolbox class.
 * It simplifies SamToolbox with less, but all useful, function provided.
 *
 * Programmer: Samuel Li
 * Date: 7/13/2015
 */

#ifndef _TOOLBOX_
#define _TOOLBOX_

#include <iostream>
#include <cstdio>

namespace VAPoR{

struct Error{
    double rms = 0.0;
    double max = 0.0;
};

class Toolbox { 
public:

    // Compare two arrays and return error evaluations
    // 'print' controls if it prints out the results 
    template <typename T>
    Error CompareArrays( const T* arr1, const T* arr2, size_t len, bool print );


    // Output a certain range of an array
    template <typename T>
    void PrintArray( const T* arr, size_t start, size_t end );


    float Findnth( const float* arr, size_t len, size_t n);

    float CalcRMS( const vector< float > &arr );
    float CalcRMS( const float* arr, size_t len );
    double CalcRMS( const SamErr* arr, size_t len );

};

}

#endif
