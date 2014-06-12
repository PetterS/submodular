//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Some tests to verify the functionality of PseudoBoolean and 
// SymmetricPseudoBoolean
//

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <new>
#include <cstdlib>
#include <ctime>
#include <map>
#include <list>
#include <random>
#include <functional>
using namespace std;


#include "Petter-Color.h"
#include "PseudoBoolean.h"
using namespace Petter;

#include "HOCR.h"
#include "graph.h"

#include "Minimizer.h"

// Random number generator
namespace {
	mt19937 engine(unsigned(time(0)&0xffffffff));
}

template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}


template<typename MINIMIZER>
void create_function(MINIMIZER& minimizer)
{
	minimizer.AddUnaryTerm(0, 23,-64);
	minimizer.AddUnaryTerm(1, 80,-34);
	minimizer.AddUnaryTerm(2, -23,-64);

	minimizer.AddPairwiseTerm(0,1, 23,64,-20,-303);
	minimizer.AddPairwiseTerm(0,2, -23,60,-12,-302);
	minimizer.AddPairwiseTerm(1,2, -20,-21,-22,-301);

	int ind[] = {0,1,2};
	int E[] = {34,21,-56,23,1,2,-15,1};
	minimizer.AddHigherTerm(3,ind,E);
}


void test_minimize()
{
	PBF<int,3> hocr;
	Petter::Minimizer<int> minimizer;

	create_function(hocr);
	create_function(minimizer);

	//
	// Minimize HOCR
	//
	int graph_constant = hocr.cnst();
	Graph<int,int,int> graph(hocr.maxID() + 1, hocr.size());
	graph.add_node(hocr.maxID() + 1);

	PBF<int,2>::VID vars[2];
	int c;
	int size;
	hocr.startEnum();
	while (hocr.get(size, vars, c))
	{
		if (size == 1) {
			// c*x(i)
			if (c >= 0) {
				graph.add_tweights(vars[0], c, 0);
			}
			else {
				graph_constant += c;
				graph.add_tweights(vars[0], 0,  -c);
			}
		}
		else {
			// c*x(i)*x(j)

			// c has no be non-positive for this function to be 
			// submodular
			ASSERT_STR(c <= 0, "Function not submodular");

			graph_constant += c;
			graph.add_tweights(vars[1], 0,  -c);

			graph.add_edge(vars[0],vars[1], -c, 0);
		}
	}

	int hocr_energy = graph_constant + graph.maxflow();

	//
	// Minimize using Minimizer
	//

	int minimizer_energy = minimizer.minimize();

	if (hocr_energy != minimizer_energy) {
		statusFailed();
		cerr << "HOCR      : " << hocr_energy << endl;
		cerr << "Minimizer : " << minimizer_energy << endl;

		throw runtime_error("Minimizer test failed");
	}
}


template<typename real>
void test_pseudoboolean()
{
	///////////////
	// Reduction //
	///////////////
	{
		Petter::PseudoBoolean<real> pb;
		
		int nTerms = 100;
		int nVars = 100;
		
		uniform_int_distribution<int> distribution(-100, 100);
		auto random_value = bind(distribution, engine);

		map< quad, bool> exists;
		for (int t=0; t < nTerms; ++t) {
			int i,j,k,l;
			do {
				i = rand()%nVars;
				j = rand()%nVars;
				k = rand()%nVars;
				l = rand()%nVars;
			} while (i>=j || i>=k || j>=k || k>=l || exists[ make_quad(i,j,k,l) ] );

			exists[ make_pair(make_pair(i,j) , make_pair(k,l)) ] = true;

			if (i<0 || i>=j || j>=k || k>=nVars) {
				cout << "(" << i << "," << j << "," << k << ")" << endl;
				throw runtime_error("ijk failed");
			}
			
			pb.add_clique(i,j,k,l, random_value(),random_value(),random_value(),random_value(),
				random_value(),random_value(),random_value(),random_value(),
				random_value(),random_value(),random_value(),random_value(),
				random_value(),random_value(),random_value(),random_value());
		}


		//Random reduction
		int i1 = 1;
		int i2 = 10;
		int i3 = 30;
		int i4 = 40;
		int i5 = 50;
		ASSERT(nVars > 50);
		int x1 = rand()%2;
		int x2 = rand()%2;
		int x3 = rand()%2;
		int x4 = rand()%2;
		int x5 = rand()%2;

		vector<label> x(nVars, -1);
		x[i1]=x1;
		x[i2]=x2;
		x[i3]=x3;
		x[i4]=x4;
		x[i5]=x5;
		Petter::PseudoBoolean<real> pb2 = pb;

		pb2.reduce(x);

		// Test that the two functions are equal
		// for the rest of the variables
		for (int iter=0; iter<=100; ++iter) {
			for (int i=0;i<nVars;++i) {
				x[i] = rand()%2;
			}
			real f2 = pb2.eval(x);

			// Compute the original value
			x[i1]=x1;
			x[i2]=x2;
			x[i3]=x3;
			x[i4]=x4;
			x[i5]=x5;
			real f1 = pb.eval(x);	
			if ( absolute(f1-f2) > 1e-6 ) {
				cout << "f1 = " << f1 << endl;
				cout << "f2 = " << f2 << endl;
				throw runtime_error("f1(x) =/= f2(x)");
			}
		}
	}

	
}

// Instantiate
template void test_pseudoboolean<double>();
template void test_pseudoboolean<int>();

