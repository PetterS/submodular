// Petter Strandmark 2014.

#include <functional>
#include <random>
using namespace std;

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "Minimizer.h"
using namespace Petter;

typedef int real;

TEST_CASE("Single_term")
{
	real C = 0;
	Graph<real, real, real> graph1(100, 100);
	graph1.add_node(2);

	// -10*x0*x1
	add_monomial_2_to_graph(C, graph1, 0, 1, -10);
	// Force both variables to 1
	graph1.add_tweights(0, 0, 1000000);
	graph1.add_tweights(1, 0, 1000000);

	CHECK((C + graph1.maxflow()) == -10);
	CHECK(graph1.what_segment(0) == 1);
	CHECK(graph1.what_segment(1) == 1);
}

TEST_CASE("Two_nodes_11")
{
	real C = 0;
	Graph<real, real, real> graph2(100, 100);
	graph2.add_node(2);

	// -2*x0*x1 - 5*x0 + 7*x1 + 1;
	add_monomial_2_to_graph(C, graph2, 0, 1, -2);
	add_monomial_1_to_graph(C, graph2, 0, -5);
	add_monomial_1_to_graph(C, graph2, 1, 7);
	C += 1;

	// Force both variables to 1
	graph2.add_tweights(0, 0, 1000000);
	graph2.add_tweights(1, 0, 1000000);

	CHECK((C + graph2.maxflow()) == (-2 - 5 + 7 + 1));
	CHECK(graph2.what_segment(0) == 1);
	CHECK(graph2.what_segment(1) == 1);
}

TEST_CASE("One_node_1")
{
	real C = 0;
	Graph<real, real, real> graph2(100, 100);
	graph2.add_node(1);
	add_monomial_1_to_graph(C, graph2, 0, -5);
	graph2.add_tweights(0, 0, 1000000);
	CHECK((C + graph2.maxflow()) == -5);
	CHECK(graph2.what_segment(0) == 1);
}

TEST_CASE("One_node_0")
{
	real C = 0;
	Graph<real, real, real> graph2(100, 100);
	graph2.add_node(1);
	add_monomial_1_to_graph(C, graph2, 0, -5);
	graph2.add_tweights(0, 1000000, 0);
	CHECK((C + graph2.maxflow()) == 0);
	CHECK(graph2.what_segment(0) == 0);
}

TEST_CASE("One_node_1_positive_cost")
{
	real C = 0;
	Graph<real, real, real> graph2(100, 100);
	graph2.add_node(1);
	add_monomial_1_to_graph(C, graph2, 0, 5);
	graph2.add_tweights(0, 0, 1000000);
	CHECK((C + graph2.maxflow()) == 5);
	CHECK(graph2.what_segment(0) == 1);
}

TEST_CASE("One_node_0_positive_cost")
{
	real C = 0;
	Graph<real, real, real> graph2(100, 100);
	graph2.add_node(1);
	add_monomial_1_to_graph(C, graph2, 0, 5);
	graph2.add_tweights(0, 1000000, 0);
	CHECK((C + graph2.maxflow()) == 0);
	CHECK(graph2.what_segment(0) == 0);
}
