// pyardrone.i
%module submodular
%{
#define SWIG_FILE_WITH_INIT
#include "submodular.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%include "typemaps.i"
%include "std_vector.i"
%include "std_string.i"

%include exception.i
%exception {
	try {
		$action
	} catch(std::exception& e) {
		// const_cast for Perl
		SWIG_exception(SWIG_RuntimeError, const_cast<char *>(e.what()));
	} catch(...) {
		SWIG_exception(SWIG_RuntimeError, "Unknown error");
	}
}

namespace std {
   %template(vectori) vector<int>;
   %template(vectord) vector<double>;
};

%apply (float* INPLACE_ARRAY1, int DIM1) {(float* data, int n)};
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* source, int n1)}
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* sink, int n2)};
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* I, int ni)};
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* J, int nj)};
%apply (unsigned char* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned char* buffer, int xsize, int ysize, int zsize)};
%apply (signed char* INPLACE_ARRAY1, int DIM1) {(signed char* x, int nx)};

%include "submodular.h"

%init %{
	// Init code
%}