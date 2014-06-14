// Petter Strandmark 2014.

#include <functional>
#include <random>
using namespace std;

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "PseudoBoolean.h"
using namespace Petter;


template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}

template<typename real>
void reduction_test()
{
	Petter::PseudoBoolean<real> pb;

	mt19937_64 engine(0ul);
	uniform_int_distribution<int> distribution(-100, 100);
	auto random_value = std::bind(distribution, engine);

	int nTerms = 100;
	int nVars = 100;

	map< quad, bool> exists;
	for (int t = 0; t < nTerms; ++t) {
		int i, j, k, l;
		do {
			i = rand() % nVars;
			j = rand() % nVars;
			k = rand() % nVars;
			l = rand() % nVars;
		} while (i >= j || i >= k || j >= k || k >= l || exists[make_quad(i, j, k, l)]);

		exists[make_pair(make_pair(i, j), make_pair(k, l))] = true;

		if (i<0 || i >= j || j >= k || k >= nVars) {
			cout << "(" << i << "," << j << "," << k << ")" << endl;
			throw runtime_error("ijk failed");
		}

		pb.add_clique(i, j, k, l, random_value(), random_value(), random_value(), random_value(),
			random_value(), random_value(), random_value(), random_value(),
			random_value(), random_value(), random_value(), random_value(),
			random_value(), random_value(), random_value(), random_value());
	}


	//Random reduction
	int i1 = 1;
	int i2 = 10;
	int i3 = 30;
	int i4 = 40;
	int i5 = 50;
	ASSERT(nVars > 50);
	int x1 = rand() % 2;
	int x2 = rand() % 2;
	int x3 = rand() % 2;
	int x4 = rand() % 2;
	int x5 = rand() % 2;

	vector<label> x(nVars, -1);
	x[i1] = x1;
	x[i2] = x2;
	x[i3] = x3;
	x[i4] = x4;
	x[i5] = x5;
	Petter::PseudoBoolean<real> pb2 = pb;

	pb2.reduce(x);

	// Test that the two functions are equal
	// for the rest of the variables
	for (int iter = 0; iter <= 100; ++iter) {
		for (int i = 0; i<nVars; ++i) {
			x[i] = rand() % 2;
		}
		real f2 = pb2.eval(x);

		// Compute the original value
		x[i1] = x1;
		x[i2] = x2;
		x[i3] = x3;
		x[i4] = x4;
		x[i5] = x5;
		real f1 = pb.eval(x);
		CAPTURE(f1);
		CAPTURE(f2)
		CHECK(absolute(f1 - f2) < 1e-6);
	}
}

TEST_CASE("Reduction<double>")
{
	reduction_test<double>();
}

TEST_CASE("Reduction<int>")
{
	reduction_test<int>();
}
