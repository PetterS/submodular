//
// Petter Strandmark 2011, 2014
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

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "Posiform.h"
using namespace Petter;

template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}

template<typename real>
void test_posiform_small()
{
		mt19937_64 engine(42ul);
		PseudoBoolean<real> pb;
		auto random_coef = bind(uniform_int_distribution<int>(-1000,1000), ref(engine));

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

			CAPTURE(f);
			CAPTURE(p);
			for (int i = 0; i<6; ++i) {
				CAPTURE(x[i]);
			}
			CHECK(absolute(f-p) < 1e-10);
		}}}}}}
}

TEST_CASE("Small")
{
	test_posiform_small<int>();
	test_posiform_small<double>();
}

template<typename real>
void test_posiform_large()
{
		mt19937_64 engine(42ul);
		Petter::PseudoBoolean<real> pb;

		int nTerms = 100;
		int nVars = 100;

		uniform_int_distribution<int> distribution(-100, 100);
		auto random_value = bind(distribution, ref(engine));

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

			CHECK( ! (i<0 || i>=j || j>=k || k>=nVars));

			pb.add_clique(i,j,k,l, random_value(),random_value(),random_value(),random_value(),
				random_value(),random_value(),random_value(),random_value(),
				random_value(),random_value(),random_value(),random_value(),
				random_value(),random_value(),random_value(),random_value());
		}

		Posiform<real,4> phi(pb);

		uniform_int_distribution<label> distribution01(0, 1);
		auto rand01 = bind(distribution01, ref(engine));

		vector<label> x(nVars,0), y(nVars,0);
		for (int iter=0; iter<=100; ++iter) {
			for (int i=0;i<nVars;++i) {
				x[i] = rand01();
			}
			real f = pb.eval(x);
			real p = -phi.eval(x);
			CAPTURE(f);
			CAPTURE(p);
			CHECK(absolute(f-p) <= 1e-6);
		}
}

TEST_CASE("Large")
{
	test_posiform_large<int>();
	test_posiform_large<double>();
}
