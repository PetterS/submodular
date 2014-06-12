// Petter Strandmark 2014.

#include <functional>
#include <random>
using namespace std;

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "PseudoBoolean.h"
using namespace Petter;

typedef double real;
template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}

TEST_CASE("Small-1")
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
		for (int i=0;i<4;++i) {
			CAPTURE(int(x[i]));
		}
		CHECK( absolute(f-g) < 1e-10 );
	}}}}
}

TEST_CASE("Small-2")
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
		for (int i = 0; i < 6; ++i) {
			CAPTURE(int(x[i]));
		}
		CHECK( absolute(f-g) < 1e-10 );
	}}}}}}
}

TEST_CASE("Large")
{
	Petter::PseudoBoolean<real> pb;

	int nTerms = 100;
	int nVars = 100;

	mt19937_64 engine(0ul);
	uniform_int_distribution<int> distribution(-100, 100);
	auto random_value = std::bind(distribution, engine);

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

		CAPTURE(i);
		CAPTURE(j);
		CAPTURE(k);
		CAPTURE(l);
		CHECK((i >= 0 && i < j && j < k && k < nVars));

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
		CAPTURE(f);
		CAPTURE(g);
		CHECK( absolute(f-g) < 1e-6 );
	}

	//int labeled;
	//double lowerbound = pb.minimize_reduction(x, labeled);
	//cout << " labeled = " << labeled  << "/" << nVars << " = " << double(labeled)/nVars*100 << "%" << endl;
	//cout << " bound = " << lowerbound << endl;

	//lowerbound = spbf.minimize(x, labeled);
	//cout << " labeled = " << labeled  << "/" << nVars << " = " << double(labeled)/nVars*100 << "%" << endl;
	//cout << " bound = " << lowerbound << endl;
}
