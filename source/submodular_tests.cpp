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

	/////////////
	// Order-2 //
	/////////////
	{
		Petter::PseudoBoolean<real> pb;

		real E00 = 1000-rand()%2001;
		real E01 = 1000-rand()%2001;
		real E10 = 1000-rand()%2001;
		real E11 = 1000-rand()%2001;

		pb.add_clique(0,1, E00, E01, E10, E11);
		vector<label> x(2);

		x[0]=0; x[1]=0;
		ASSERT( absolute( pb.eval(x) - E00 ) < 1e-10 );
		x[0]=0; x[1]=1;
		ASSERT( absolute( pb.eval(x) - E01 ) < 1e-10 );
		x[0]=1; x[1]=0;
		ASSERT( absolute( pb.eval(x) - E10 ) < 1e-10 );
		x[0]=1; x[1]=1;
		ASSERT( absolute( pb.eval(x) - E11 ) < 1e-10 );
	}

	/////////////
	// Order-3 //
	/////////////
	{
		Petter::PseudoBoolean<real> pb;

		vector<real> E(8,0);
		for (int i=0;i<8;++i) {
			E[i] = 1000-rand()%2001;
		}

		pb.add_clique(0,1,2, E);
		pb.add_clique(0,1,2, E);

		vector<label> x(3);
		for (x[0]=0;x[0]<=1;++x[0]) {
		for (x[1]=0;x[1]<=1;++x[1]) {
		for (x[2]=0;x[2]<=1;++x[2]) {
			real f = pb.eval(x);
			real e = 2*E[4*x[0] + 2*x[1] + x[2]];
			if ( absolute(f-e) > 1e-10 ) {
				throw runtime_error("f =/= E for order 3");
			}
		}}}
	}

	/////////////
	// Order-4 //
	/////////////
	{
		Petter::PseudoBoolean<real> pb;

		vector<real> E(16,0);
		for (int i=0;i<16;++i) {
			E[i] = 1000-rand()%2001;
		}

		pb.add_clique(0,1,2,3, E);
		pb.add_clique(0,1,2,3, E);

		vector<label> x(4);
		for (x[0]=0;x[0]<=1;++x[0]) {
		for (x[1]=0;x[1]<=1;++x[1]) {
		for (x[2]=0;x[2]<=1;++x[2]) {
		for (x[3]=0;x[3]<=1;++x[3]) {
			real f = pb.eval(x);
			real e = 2*E[8*x[0] + 4*x[1] + 2*x[2] + x[3]];
			if ( absolute(f-e) > 1e-10 ) {
				throw runtime_error("f =/= E for order 4");
			}
		}}}}
	}




	////////////////
	// Relaxation //
	////////////////
	{
		Petter::SymmetricPseudoBoolean<real> spbf;
		spbf.test_clear();
		spbf.test_clear();

		Petter::PseudoBoolean<real> pb;

		vector<real> E(8,0);
		for (int i=0;i<8;++i) {
			E[i] = 1000-rand()%2001;
		}
		pb.add_clique(0,1,2, E);
		
		for (int i=0;i<8;++i) {
			E[i] = 1000-rand()%2001;
		}
		pb.add_clique(1,2,3, E);

		spbf.create_lp(pb);
		spbf.test();

		vector<label> x(4,0), y(4,0);
		spbf.eval(x,y);

		int labeled;
		pb.minimize_reduction(x, labeled);
		
		//
		// Verify that f(x) = g(x,bar(x))
		//
		for (x[0]=0;x[0]<=1;++x[0]) {
		for (x[1]=0;x[1]<=1;++x[1]) {
		for (x[2]=0;x[2]<=1;++x[2]) {
		for (x[3]=0;x[3]<=1;++x[3]) {
			real f = pb.eval(x);
			vector<label> y(6);
			y[0] = 1-x[0];
			y[1] = 1-x[1];
			y[2] = 1-x[2];
			y[3] = 1-x[3];
			real g = spbf.eval(x,y);
			if ( absolute(f-g) > 1e-10 ) {
				cout << endl;
				for (int i=0;i<4;++i) {
					cout << int(x[i]);
				}
				cout << endl;
				throw runtime_error("f(x) =/= g(x,bar(x)) for order 3");
			}
		}}}}
	}


	////////////////
	// Relaxation //
	////////////////
	{
		Petter::SymmetricPseudoBoolean<real> spbf;
		spbf.test_clear();
		spbf.test_clear();

		Petter::PseudoBoolean<real> pb;

		vector<real> E(16,0);
		for (int i=0;i<16;++i) {
			E[i] = 1000-rand()%2001;
		}
		pb.add_clique(0,1,2,3, E);
		
		for (int i=0;i<16;++i) {
			E[i] = 1000-rand()%2001;
		}
		pb.add_clique(2,3,4,5, E);

		spbf.create_lp(pb);
		spbf.test();

		vector<label> x(6,0), y(6,0);
		spbf.eval(x,y);

		int labeled;
		pb.minimize_reduction(x, labeled);
		
		//
		// Verify that f(x) = g(x,bar(x))
		//
		for (x[0]=0;x[0]<=1;++x[0]) {
		for (x[1]=0;x[1]<=1;++x[1]) {
		for (x[2]=0;x[2]<=1;++x[2]) {
		for (x[3]=0;x[3]<=1;++x[3]) {
		for (x[4]=0;x[4]<=1;++x[4]) {
		for (x[5]=0;x[5]<=1;++x[5]) {
			real f = pb.eval(x);
			vector<label> y(6);
			y[0] = 1-x[0];
			y[1] = 1-x[1];
			y[2] = 1-x[2];
			y[3] = 1-x[3];
			y[4] = 1-x[4];
			y[5] = 1-x[5];
			real g = spbf.eval(x,y);
			if ( absolute(f-g) > 1e-10 ) {
				cout << endl;
				for (int i=0;i<6;++i) {
					cout << int(x[i]);
				}
				cout << endl;
				throw runtime_error("f(x) =/= g(x,bar(x)) for order 4");
			}
		}}}}}}
	}


	////////////
	// Larger //
	////////////
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

		Petter::SymmetricPseudoBoolean<real> spbf;

		spbf.create_lp(pb);

		vector<label> x(nVars,0), y(nVars,0);
		for (int iter=0; iter<=100; ++iter) {
			for (int i=0;i<nVars;++i) {
				x[i] = rand()%2;
				y[i] = 1-x[i];
			}
			real f = pb.eval(x);
			real g = spbf.eval(x,y);
			if ( absolute(f-g) > 1e-6 ) {
				cout << "f = " << f << endl;
				cout << "g = " << g << endl;
				throw runtime_error("f(x) =/= g(x,bar(x)) for order 4");
			}
		}



		//int labeled;
		//double lowerbound = pb.minimize_reduction(x, labeled);
		//cout << " labeled = " << labeled  << "/" << nVars << " = " << double(labeled)/nVars*100 << "%" << endl;
		//cout << " bound = " << lowerbound << endl;

		//lowerbound = spbf.minimize(x, labeled);
		//cout << " labeled = " << labeled  << "/" << nVars << " = " << double(labeled)/nVars*100 << "%" << endl;
		//cout << " bound = " << lowerbound << endl;
	}

	

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

