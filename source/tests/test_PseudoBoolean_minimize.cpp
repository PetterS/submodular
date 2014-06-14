// Petter Strandmark 2014.

#include <functional>
#include <random>
using namespace std;

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "HOCR.h"
#include "graph.h"

#include "PseudoBoolean.h"
#include "Minimizer.h"
using namespace Petter;


template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}

template<typename MINIMIZER>
void create_function(MINIMIZER& minimizer)
{
	minimizer.AddUnaryTerm(0, 23, -64);
	minimizer.AddUnaryTerm(1, 80, -34);
	minimizer.AddUnaryTerm(2, -23, -64);

	minimizer.AddPairwiseTerm(0, 1, 23, 64, -20, -303);
	minimizer.AddPairwiseTerm(0, 2, -23, 60, -12, -302);
	minimizer.AddPairwiseTerm(1, 2, -20, -21, -22, -301);

	int ind[] = { 0, 1, 2 };
	int E[] = { 34, 21, -56, 23, 1, 2, -15, 1 };
	minimizer.AddHigherTerm(3, ind, E);
}

TEST_CASE("Minimization")
{
	PBF<int, 3> hocr;
	Petter::Minimizer<int> minimizer;

	create_function(hocr);
	create_function(minimizer);

	//
	// Minimize HOCR
	//
	int graph_constant = hocr.cnst();
	Graph<int, int, int> graph(hocr.maxID() + 1, hocr.size());
	graph.add_node(hocr.maxID() + 1);

	PBF<int, 2>::VID vars[2];
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
				graph.add_tweights(vars[0], 0, -c);
			}
		}
		else {
			// c*x(i)*x(j)

			// c has no be non-positive for this function to be 
			// submodular
			REQUIRE(c <= 0);

			graph_constant += c;
			graph.add_tweights(vars[1], 0, -c);

			graph.add_edge(vars[0], vars[1], -c, 0);
		}
	}

	int hocr_energy = graph_constant + graph.maxflow();
	CAPTURE(hocr_energy);

	//
	// Minimize using Minimizer
	//

	int minimizer_energy = minimizer.minimize();
	CAPTURE(minimizer_energy);

	CHECK(hocr_energy == minimizer_energy);
}
