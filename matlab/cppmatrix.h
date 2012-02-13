//
// CPPMATRIX
//
// C++ wrapper for Matlab matrices 
//
// Petter Strandmark 2010
//


#ifndef CPPMATRIX_MATLAB_HEADER
#define CPPMATRIX_MATLAB_HEADER

#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "mex.h"
#include "mexutils.h"


// Extremely annoying macros needs to be undefed in
// order to make std::min and std::max work.
#ifdef min 
#undef min
#endif
#ifdef max 
#undef max
#endif


namespace {
	template<typename T>
	mxClassID typeToID()
	{
		mexErrMsgTxt("Unknown type!");
		return mxUINT8_CLASS;
	}

	template<>
	mxClassID typeToID<double>()
	{
		return mxDOUBLE_CLASS;
	}
	template<>
	mxClassID typeToID<float>()
	{
		return mxSINGLE_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned char>()
	{
		return mxUINT8_CLASS;
	}
	template<>
	mxClassID typeToID<signed char>()
	{
		return mxINT8_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned int>()
	{
		ASSERT(sizeof(unsigned int)==4);
		return mxUINT32_CLASS;
	}
	template<>
	mxClassID typeToID<signed int>()
	{
		ASSERT(sizeof(signed int)==4);
		return mxINT32_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned short>()
	{
		ASSERT(sizeof(unsigned short)==4);
		return mxUINT16_CLASS;
	}
	template<>
	mxClassID typeToID<signed short>()
	{
		ASSERT(sizeof(signed short)==4);
		return mxINT16_CLASS;
	}
}



template<typename T>
class matrix
{
public:
    
	T* data;
    mxArray* array;
	mwSize M,N,O,P;
	
	matrix(const mxArray* array)
	{
		this->array = (mxArray*)array; //Hack to allow const mxArrays
		ASSERT(!mxIsSparse(array));		
		ASSERT(mxGetClassID(array) == typeToID<T>());

		int ndim = mxGetNumberOfDimensions(array);
		ASSERT(ndim<=4);

		if (ndim <= 2) {
			M = mxGetM(array);
			N = mxGetN(array);	
			O = 1;
			P = 1;
		}
		else if (ndim==3) {
			const mwSize* dims = mxGetDimensions(array);
			M = dims[0];
			N = dims[1];
			O = dims[2];
			P = 1;
		}
		else {
			const mwSize* dims = mxGetDimensions(array);
			M = dims[0];
			N = dims[1];
			O = dims[2];
			P = dims[3];
		}
		data = (T*)mxGetPr(array);
		
		shouldDestroy = false;
	}

	matrix(mwSize M, mwSize N=1, mwSize O=1, mwSize P=1)
	{
		this->M = M;
		this->N = N;
		this->O = O;
		this->P = P;
		if (O==1 && P==1) {
			mwSize size[] = {M,N};
			this->array = mxCreateNumericArray(2,size,typeToID<T>(),mxREAL);
		}
		else if (P==1) {
			mwSize size[] = {M,N,O};
			this->array = mxCreateNumericArray(3,size,typeToID<T>(),mxREAL);
		}
		else {
			mwSize size[] = {M,N,O,P};
			this->array = mxCreateNumericArray(4,size,typeToID<T>(),mxREAL);
		}
		data = (T*)mxGetPr(array);
		
		shouldDestroy = true;
	}
	
	matrix(const matrix& org) 
	{
		*this = org;
	}
	
	matrix()
	{
		M = N = O = P = 0;
		data = 0;
		array = 0;
	}
	
	~matrix()
	{
		if (shouldDestroy) {
			mxDestroyArray(array);
		}
	}
	
	mwSize numel() 
	{
		return M*N*O*P;
	}
	
	void operator=(const matrix& org) 
	{
		if (org.shouldDestroy) {
			throw std::runtime_error("matrix() : Cannot copy a managed matrix");
		}
		shouldDestroy = false;
		M = org.M;
		N = org.N;
		O = org.O;
		P = org.P;
		data = org.data;
		array = org.array;
	}
	
	T& operator[](mwSize i)
	{
		ASSERT(i>=0 && i<numel());
		return data[i];
	}
	T& operator()(mwSize i)
	{
		return operator[](i);
	}
	T& operator()(mwSize i, mwSize j)
	{
		return operator[](i + M*j);
	}
	T& operator()(mwSize i, mwSize j, mwSize k)
	{
		return operator[](i + M*j + M*N*k);
	}
	T& operator()(mwSize i, mwSize j, mwSize k, mwSize l)
	{
		return operator[](i + M*j + M*N*k + M*N*O*l);
	}
	
	operator mxArray*()
	{
		shouldDestroy = false;
		return array;
	}
    
    T min() 
    {
        return *std::min_element(data, data+numel());
    }
    
    T max()
    {
        return *std::max_element(data, data+numel());
    }
	
private:
	
	bool shouldDestroy;

};






#endif
