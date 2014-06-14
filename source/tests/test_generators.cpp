// Petter Strandmark 2014.

#include <functional>
#include <random>
using namespace std;

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "PseudoBoolean.h"
using namespace Petter;

TEST_CASE("Generators")
{
	Generators<double> generators("generators/generators.txt");
	Petter::GeneratorPseudoBoolean<double> genpb(generators);

	cerr << "Number of generators : (" << generators.ngen2 << ", " << generators.ngen3 << ", " << generators.ngen4 << ")" << endl;
	PseudoBoolean<double> f("tests/quartic_paper.txt");
	genpb.create_lp(f);

	vector<label> x(f.nvars(), -1);
	int nlabelled=-1;
	genpb.minimize(x, nlabelled);
	CHECK(nlabelled == 4);
}
