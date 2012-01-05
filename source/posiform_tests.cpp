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
#include "Posiform.h"
using namespace Petter;




// Random number generator
namespace {
	mt19937 engine(unsigned(time(0)));
}

template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}


template<typename real>
void test_posiform()
{
	////////////////
	// Small test //
	////////////////
	{
		PseudoBoolean<real> pb;
		auto random_coef = bind(uniform_int_distribution<int>(-1000,1000), engine);

		vector<real> E(16,0);
		for (int i=0;i<16;++i) {
			E[i] = random_coef();
		}
		pb.add_clique(0,1,2,3, E);
		
		for (int i=0;i<16;++i) {
			E[i] = random_coef();
		}
		pb.add_clique(2,3,4,5, E);

		Posiform<real,4> phi(pb, true);

		vector<label> x(6,0);
		
		//
		// Verify that f(x) = phi(x)
		//
		for (x[0]=0;x[0]<=1;++x[0]) {
		for (x[1]=0;x[1]<=1;++x[1]) {
		for (x[2]=0;x[2]<=1;++x[2]) {
		for (x[3]=0;x[3]<=1;++x[3]) {
		for (x[4]=0;x[4]<=1;++x[4]) {
		for (x[5]=0;x[5]<=1;++x[5]) {
			real f = pb.eval(x);
			real p = -phi.eval(x);
			if ( absolute(f-p) > 1e-10 ) {
				cout << endl;
				for (int i=0;i<6;++i) {
					cout << int(x[i]);
				}
				cout << endl;
				cout << "f   = " << f << endl;
				cout << "phi = " << p << endl;
				throw runtime_error("f(x) =/= phi(x)");
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

		Posiform<real,4> phi(pb);

		vector<label> x(nVars,0), y(nVars,0);
		for (int iter=0; iter<=100; ++iter) {
			for (int i=0;i<nVars;++i) {
				x[i] = rand()%2;
			}
			real f = pb.eval(x);
			real p = -phi.eval(x);
			if ( absolute(f-p) > 1e-6 ) {
				cout << "f   = " << f << endl;
				cout << "phi = " << p << endl;
				throw runtime_error("f(x) =/= phi(x) for order 4");
			}
		}

	}

}

// Instantiate
template void test_posiform<double>();
template void test_posiform<int>();

