//
// Petter Strandmark 2012
//
// Using an interface similar to vgg_qpbo
// (available in Oliver Woodford's imrender toolbox)
//

#include <cstdlib>
#include <stdexcept>
#include <iostream>
using namespace std;

#include "mex.h"
#include "cppmatrix.h"
#include "mexutils.h"

#include <PseudoBoolean.h>

template<typename real>
void grd(int			nlhs, 		/* number of expected outputs */
	mxArray		*plhs[],	/* mxArray output pointer array */
	int			nrhs, 		/* number of inputs */
	const mxArray	*prhs[]		/* mxArray input pointer array */)
{
	using namespace Petter;
	ASSERT_STR(nrhs >= 7, "Need arguments");

	//First param
	int curarg = 0;
	matrix<real> UE(prhs[curarg++]);
	matrix<unsigned int>  PI(prhs[curarg++]);
	matrix<real> PE(prhs[curarg++]);
	matrix<unsigned int>  TI(prhs[curarg++]);
	matrix<real> TE(prhs[curarg++]);
	matrix<unsigned int>  QI(prhs[curarg++]);
	matrix<real> QE(prhs[curarg++]);

	//Additional params
	MexParams params(nrhs-curarg, prhs+curarg);

	//Get some params
	string method_str = params.get<string>("method","GRD"); 
	string generators = params.get<string>("generators","../bin/generators/generators.txt");
	//double d = params.get<double>("d",42); 
	//vector<double> vec = params.get< vector<double> >("vec");

	size_t n = UE.N;
	ASSERT(UE.M == 2);
	ASSERT(PI.numel() == 0 || PI.M == 2 || PI.M == 3);
	ASSERT(PI.numel() == 0 || PE.M == 4);
	ASSERT(TI.numel() == 0 || TI.M == 3 || TI.M == 4);
	ASSERT(TI.numel() == 0 || TE.M == 8);
	ASSERT(QI.numel() == 0 || QI.M == 4 || QI.M == 5);
	ASSERT(QI.numel() == 0 || QE.M == 16);


	PseudoBoolean<real> f;

	//
	// Unary terms 
	//
	for (size_t i=0; i<n; ++i) {
		f.add_clique(i, UE(0,i), UE(1,i));
	}

	//
	// Quadratic terms
	//
	for (size_t t=0; t<PI.N; ++t) {
		unsigned int i   = PI(0,t)-1; 
		unsigned int j   = PI(1,t)-1; 
		unsigned int ind;
		if (PI.M == 3) {
			ind = PI(2,t)-1;
		}
		else {
			ind = t;
		}
		ASSERT( i < n && j < n );
		ASSERT(ind < PE.N);
		f.add_clique(i,j, PE(0,ind), PE(1,ind), PE(2,ind), PE(3,ind));
		//mexPrintf("(%d,%d) : %g %g %g %g\n",i,j, PE(0,ind), PE(1,ind), PE(2,ind), PE(3,ind));
	}

	//
	// Ternary terms
	//
	for (size_t t=0; t<TI.N; ++t) {
		unsigned int i   = TI(0,t)-1; 
		unsigned int j   = TI(1,t)-1; 
		unsigned int k   = TI(2,t)-1; 
		unsigned int ind;
		if (TI.M == 4) {
			ind = TI(3,t)-1;
		}
		else {
			ind = t;
		}
		ASSERT( i < n && j < n && k < n);
		ASSERT(ind < TE.N);
		f.add_clique(i,j,k, TE(0,ind), TE(1,ind), TE(2,ind), TE(3,ind), 
		                    TE(4,ind), TE(5,ind), TE(6,ind), TE(7,ind));
	}

	//
	// Quartic terms
	//
	for (size_t t=0; t<QI.N; ++t) {
		unsigned int i   = QI(0,t)-1; 
		unsigned int j   = QI(1,t)-1; 
		unsigned int k   = QI(2,t)-1; 
		unsigned int l   = QI(3,t)-1; 
		unsigned int ind;
		if (QI.M == 5) {
			ind = QI(4,t)-1;
		}
		else {
			ind = t;
		}
		ASSERT( i < n && j < n && k < n && l < n );
		ASSERT(ind < QE.N);
		f.add_clique(i,j,k,l, QE(0,ind), QE(1,ind), QE(2,ind), QE(3,ind), 
		                      QE(4,ind), QE(5,ind), QE(6,ind), QE(7,ind),
		                      QE(8,ind), QE(9,ind), QE(10,ind), QE(11,ind), 
		                      QE(12,ind), QE(13,ind), QE(14,ind), QE(15,ind));
	}

	Method method = HOCR;
	if (method_str == "GRD") {
		method = GRD;
	}
	else if (method_str == "HOCR") {
		method = HOCR;
	}
	else if (method_str == "GRD-heur") {
		method = GRD_heur;
	}
	else if (method_str == "Fix") {
		method = Fix;
	}
	else if (method_str == "GRD-gen") {
		method = GRD_gen;
	}
	else {
		throw runtime_error("Unknown optimization method");
	}

	vector<label> x(n,-1);
	matrix<real> energy(1);
	energy(0) = f.minimize(x, method, generators.c_str() );

	matrix<int> lab(n);
	for (size_t i=0; i<n; ++i) {
		lab[i] = x.at(i); 
	}

	if (nlhs >= 1) {
		plhs[0] = lab;
	}
	if (nlhs >= 2) {
		plhs[1] = energy;
	}
}


void mexFunctionReal(int nlhs, mxArray *plhs[],int nrhs, const mxArray	*prhs[])
{
	ASSERT_STR(nrhs >= 1, "Need arguments");

	switch (mxGetClassID(prhs[0])) {
	case mxINT32_CLASS:
		grd<int>(nlhs, plhs, nrhs, prhs);
		break;
	case mxDOUBLE_CLASS:
		grd<double>(nlhs, plhs, nrhs, prhs);
		break;
	default:
		throw runtime_error("Inputs are of an unsupported type");
		break;
	}
}

void mexFunction(int		nlhs, 		/* number of expected outputs */
	mxArray	*plhs[],	/* mxArray output pointer array */
	int		nrhs, 		/* number of inputs */
	const mxArray	*prhs[]		/* mxArray input pointer array */)
{
	try {
		mexFunctionReal(nlhs,plhs,nrhs,prhs);
	}
	catch (bad_alloc& ) {
		mexErrMsgTxt("Out of memory");
	}
	catch (exception& e) {
		mexErrMsgTxt(e.what());
	}
	catch (...) {
		mexErrMsgTxt("Unknown exception");
	}
}
