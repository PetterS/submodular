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

TEST_CASE("PseudoBoolean-2")
{
	mt19937_64 engine(1ul);
	auto rand_value = bind(uniform_int_distribution<int>(-1000, 1000), ref(engine));

	real E00 = rand_value();
	real E01 = rand_value();
	real E10 = rand_value();
	real E11 = rand_value();

	Petter::PseudoBoolean<real> pb;
	pb.add_clique(0,1, E00, E01, E10, E11);
	vector<label> x(2);

	x[0]=0; x[1]=0;
	CHECK( absolute( pb.eval(x) - E00 ) < 1e-10 );
	x[0]=0; x[1]=1;
	CHECK( absolute( pb.eval(x) - E01 ) < 1e-10 );
	x[0]=1; x[1]=0;
	CHECK( absolute( pb.eval(x) - E10 ) < 1e-10 );
	x[0]=1; x[1]=1;
	CHECK( absolute( pb.eval(x) - E11 ) < 1e-10 );
}

TEST_CASE("PseudoBoolean-3")
{
	mt19937_64 engine(1ul);
	auto rand_value = bind(uniform_int_distribution<int>(-1000, 1000), ref(engine));

	vector<real> E(8,0);
	for (int i=0;i<8;++i) {
		E[i] = rand_value();
	}

	Petter::PseudoBoolean<real> pb;
	pb.add_clique(0,1,2, E);
	pb.add_clique(0,1,2, E);

	vector<label> x(3);
	for (x[0]=0;x[0]<=1;++x[0]) {
	for (x[1]=0;x[1]<=1;++x[1]) {
	for (x[2]=0;x[2]<=1;++x[2]) {
		real f = pb.eval(x);
		real e = 2*E[4*x[0] + 2*x[1] + x[2]];
		CHECK( absolute(f-e) <= 1e-10 );
	}}}
}

TEST_CASE("PseudoBoolean-4")
{
	mt19937_64 engine(1ul);
	auto rand_value = bind(uniform_int_distribution<int>(-1000, 1000), ref(engine));

	vector<real> E(16,0);
	for (int i=0;i<16;++i) {
		E[i] = rand_value();
	}

	Petter::PseudoBoolean<real> pb;
	pb.add_clique(0,1,2,3, E);
	pb.add_clique(0,1,2,3, E);

	vector<label> x(4);
	for (x[0]=0; x[0]<=1; ++x[0]) {
	for (x[1]=0; x[1]<=1; ++x[1]) {
	for (x[2]=0; x[2]<=1; ++x[2]) {
	for (x[3]=0; x[3]<=1; ++x[3]) {
		real f = pb.eval(x);
		real e = 2*E[8*x[0] + 4*x[1] + 2*x[2] + x[3]];
		CHECK( absolute(f-e) <= 1e-10 );
	}}}}
}
