//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Minimizes a submodular function of degree <= 3. Each term does not have to 
// be submodular, but the complete function has to be.
//
// The code is suited for general purposes, except for resolve_different(), 
// which is used for minimizing symmetric g(x,y)
//

#ifndef PETTER_MINIMIZER_H
#define PETTER_MINIMIZER_H

#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>

// Vladimir Kolmogorov's maximum flow code
#include "graph.h"

#ifndef ASSERT
#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
                                    sout << "Error (line " << __LINE__ << " in " << __FILE__ << "): " << #cond; \
                                    throw std::runtime_error(sout.str()); }
#endif

namespace {
	void error_function(char* msg) 
	{
		throw std::runtime_error(msg);
	}
}


namespace Petter {

	template<typename T>
	std::pair<T,T> unordered_pair(T i, T j)
	{
		if (j<i) {
			return std::make_pair(j,i);
		}
		else {
			return std::make_pair(i,j);
		}
	}



	//
	// Adds a linear monomial to a graph
	//
	template<typename real>
	void add_monomial_1_to_graph( real& c, Graph<real,real,real>& graph, int i, real a) 
	{
		if (a<0) {
			// a*xi = -a*(1-xi) + a
			graph.add_tweights(i, 0,  -a);
			c += a;
		}
		else {
			graph.add_tweights(i, a,  0);
		}
	}

	//
	// Adds a quadratic monomial to a graph
	//
	template<typename real>
	void add_monomial_2_to_graph( real& c,  Graph<real,real,real>& graph, int i, int j, real a) 
	{
		ASSERT(a <= 0); // Submodularity

		// a*xi*xj = -a*(1-xi)*xj + a*xj
		graph.add_edge(i,j, -a, 0);
		add_monomial_1_to_graph(c, graph, j, a);
	}



	//
	// Takes a list of pairs (i,j) and tries to find a solution where
	// x[i] != x[j]. The function is assumed to be symmetric, i.e. 
	// (1,1) is never a unique optimal solution.
	//
	template<typename real>
	void resolve_different(Graph<real,real,real>& graph, std::vector<char>& x, const std::vector<std::pair<int,int> >& pairs)
	{
		real infinity = 100; //Does not have to be big

		// Fix the first element of the pair
		for (auto itr=pairs.begin(); itr != pairs.end(); ++itr) {
			int i = itr->first;
			int j = itr->second;

			// Extract solution
			x.at(i) = graph.what_segment(i, Graph<real,real,real>::SOURCE);
			x.at(j) = graph.what_segment(j, Graph<real,real,real>::SOURCE);

			//std::cout << i << ',' << j << " : " << int(x[i]) << ',' << int(x[j]) << "  ";

			//ASSERT(x[i]==0 || x[j]==0); // Should not happen, but might due to rounding errors

			// If x[i]==0 and y[i]==0, we don't know whether there is another 
			// solution where the two are different
			if (x[i]==0 && x[j]==0) {

				// Is there also an optimal solution for which x[i]==1?
				if (graph.what_segment(i, Graph<real,real,real>::SINK) == 1) {
					// Then fix x[i] to 1
					graph.add_tweights(i, 0, infinity);
					graph.mark_node(i);
					// Resolve the graph
					graph.maxflow(true);
					// Extract solution; different if possible
					x[i] = graph.what_segment(i, Graph<real,real,real>::SOURCE);
					x[j] = graph.what_segment(j, Graph<real,real,real>::SOURCE);
					ASSERT(x[i] == 1);

					if (x[j] == 0) {
						// We found a solution where x[i]==1 and x[j]==0
						// std::cout << "A\n";
						continue;
					}
				}

				// Now fix x[i] to 0
				graph.add_tweights(i, 2*infinity, 0);
				graph.mark_node(i);
				// Resolve the graph
				graph.maxflow(true);
				// Does there then exist a solution for which x[j]==1?
				if (graph.what_segment(j, Graph<real,real,real>::SINK) == 1) {
					// Then use that solution
					graph.add_tweights(j, 0, infinity);
					graph.mark_node(j);
					// Resolve the graph
					graph.maxflow(true);
					x[i] = graph.what_segment(i, Graph<real,real,real>::SOURCE);
					x[j] = graph.what_segment(j, Graph<real,real,real>::SOURCE);
					ASSERT(x[i]==0 && x[j]==1);
					// We found a solution where x[i]==0 and x[j]==1
					// std::cout << "B\n";
					continue;
				}

				// Otherwise, we only have the solution where x[i]==0 and x[j]==0
				// std::cout << "C\n";
				x[i] = 0;
				x[j] = 0;
			}
			else {
				// std::cout << '\n';
			}
		}

	}


	
	template<typename real>
	void test_graph_functions() 
	{
		{
			//
			// First test
			//
			real C = 0;
			Graph<real,real,real> graph1(100,100);
			graph1.add_node(2);
		
			// -10*x0*x1
			add_monomial_2_to_graph(C, graph1, 0,1, -10);
			// Force both variables to 1
			graph1.add_tweights(0,0,1000000);
			graph1.add_tweights(1,0,1000000);

			ASSERT( C+graph1.maxflow() == -10 );
			ASSERT( graph1.what_segment(0) == 1);
			ASSERT( graph1.what_segment(1) == 1);
		}

		{
			//
			// Second test
			//
			real C = 0;
			Graph<real,real,real> graph2(100,100);
			graph2.add_node(2);
		
			// -2*x0*x1 - 5*x0 + 7*x1 + 1;
			add_monomial_2_to_graph(C,graph2, 0,1, -2);
			add_monomial_1_to_graph(C,graph2, 0, -5);
			add_monomial_1_to_graph(C,graph2, 1,  7);
			C += 1;

			// Force both variables to 1
			graph2.add_tweights(0,0,1000000);
			graph2.add_tweights(1,0,1000000);
		
			ASSERT( C+graph2.maxflow() == -2 - 5 + 7 + 1);
			ASSERT( graph2.what_segment(0) == 1);
			ASSERT( graph2.what_segment(1) == 1);
		}

		{
			real C = 0;
			Graph<real,real,real> graph2(100,100);
			graph2.add_node(1);
			add_monomial_1_to_graph(C,graph2, 0, -5);
			graph2.add_tweights(0,0,1000000);		
			ASSERT( C+graph2.maxflow() == -5);
			ASSERT( graph2.what_segment(0) == 1);
		}

		{
			real C = 0;
			Graph<real,real,real> graph2(100,100);
			graph2.add_node(1);
			add_monomial_1_to_graph(C,graph2, 0, -5);
			graph2.add_tweights(0,1000000,0);		
			ASSERT( C+graph2.maxflow() == 0);
			ASSERT( graph2.what_segment(0) == 0);
		}

		{
			real C = 0;
			Graph<real,real,real> graph2(100,100);
			graph2.add_node(1);
			add_monomial_1_to_graph(C,graph2, 0, 5);
			graph2.add_tweights(0,0,1000000);		
			ASSERT( C+graph2.maxflow() == 5);
			ASSERT( graph2.what_segment(0) == 1);
		}

		{
			real C = 0;
			Graph<real,real,real> graph2(100,100);
			graph2.add_node(1);
			add_monomial_1_to_graph(C,graph2, 0, 5);
			graph2.add_tweights(0,1000000,0);		
			ASSERT( C+graph2.maxflow() == 0);
			ASSERT( graph2.what_segment(0) == 0);
		}
	}









	template<typename real>
	class Minimizer 
	{
	public:
		Minimizer(int nvars=0) 
		{
			this->graph    = 0;
			this->n = nvars;
			this->new_vars = 0;
			this->terms0   = 0;
		}

		~Minimizer() 
		{
			if (graph) {
				delete graph;
			}
		}

		void AddUnaryTerm(int i, real E0, real E1)
		{
			// Update number of variables
			if (i+1>n) n=i+1;

			terms0 += E0;
			terms1[i] += E1-E0;
		}

		void AddPairwiseTerm(int i, int j, real E00, real E01, real E10, real E11)
		{
			if (j<i) {
				AddPairwiseTerm(j,i, E00, E10, E01, E11);
				return;
			}

			// Update number of variables
			if (i+1>n) n=i+1;
			if (j+1>n) n=j+1;

			terms0 += E00;
			real ai = E10 - E00;
			real aj = E01 - E00;
			real aij = E11 - ai - aj - E00; 

			terms1[i] += ai;
			terms1[j] += aj;
			terms2[  unordered_pair(i,j)] += aij;
		}

		void AddHigherTerm(int nvars_add, int ind[], real E[])
		{
			ASSERT(nvars_add==3);

			real& E000 = E[0];
			real& E001 = E[1];
			real& E010 = E[2];
			real& E011 = E[3];
			real& E100 = E[4];
			real& E101 = E[5];
			real& E110 = E[6];
			real& E111 = E[7];

			int& i = ind[0];
			int& j = ind[1];
			int& k = ind[2];

			// Update number of variables
			if (i+1>n) n=i+1;
			if (j+1>n) n=j+1;
			if (k+1>n) n=k+1;

			terms0 += E000;

			real ai = E100 - E000;
			real aj = E010 - E000;
			real ak = E001 - E000;

			real aij = E110 - ai - aj - E000;
			real aik = E101 - ai - ak - E000;
			real ajk = E011 - aj - ak - E000;

			real aijk = E111 - aij - aik - ajk - ai - aj - ak - E000;

			terms1[i] += ai;
			terms1[j] += aj;
			terms1[k] += ak;
			terms2[unordered_pair(i,j)] += aij;
			terms2[unordered_pair(i,k)] += aik;
			terms2[unordered_pair(j,k)] += ajk;

			if (aijk != 0) {
				int l = new_vars++;

				if (aijk < 0) { 
					// Use reduction
					// -xi*xj*xk = min{z} z*(2-xi-xj-xk) 
					terms2_new[ std::make_pair(i,l)] += aijk;
					terms2_new[ std::make_pair(j,l)] += aijk;
					terms2_new[ std::make_pair(k,l)] += aijk;
					terms1_new[ l ] += -2*aijk;
				}
				else if (aijk > 0) { 
					// Use reduction
					// xi*xj*xk = min{z} z*(1-xi-xj-xk)  + xi*xj + xi*xk + xj*xk
					terms2_new[ std::make_pair(i,l)] += -aijk;
					terms2_new[ std::make_pair(j,l)] += -aijk;
					terms2_new[ std::make_pair(k,l)] += -aijk;
					terms1_new[l] += aijk;
					terms2[ unordered_pair(i,j)] += aijk;
					terms2[ unordered_pair(i,k)] += aijk;
					terms2[ unordered_pair(j,k)] += aijk;
				}
			}
		}


		real minimize()
		{
			using namespace std;

			ASSERT(!graph);
			int n_nodes = n + new_vars;
			graph = new Graph<real,real,real>(n_nodes,int(terms2.size()),error_function);
			graph->add_node(n_nodes);

			// Add degree-1 terms
			for (int i=0; i<n; ++i) {
				real c = terms1[i];

				// c*x(i)
				if (c >= 0) {
					graph->add_tweights(i, c, 0);
				}
				else {
					terms0 += c;
					graph->add_tweights(i, 0,  -c);
				}
			}
			// Add degree-1 terms for new variables
			for (int ind=0; ind<new_vars; ++ind) {
				real c = terms1_new[ind];
				int i = n + ind;
				if (c >= 0) {
					graph->add_tweights(i, c, 0);
				}
				else {
					terms0 += c;
					graph->add_tweights(i, 0,  -c);
				}
			}

			// Add degree-2 terms
			for (auto itr = terms2.begin(); itr != terms2.end(); ++itr) {
				int i,j;
				std::tie(i,j) = itr->first;
				real coef     = submodular_coef(itr->second);

				// Add to graph
				terms0 += coef;
				graph->add_tweights(j, 0,  -coef);
				graph->add_edge(i,j, -coef, 0);
			}
			// Add degree-2 terms for new variables
			for (auto itr = terms2_new.begin(); itr != terms2_new.end(); ++itr) {
				int i,ind;
				std::tie(i,ind) = itr->first;
				real coef       = submodular_coef(itr->second);
				int j = n+ind;

				// Add to graph
				terms0 += coef;
				graph->add_tweights(j, 0,  -coef);
				graph->add_edge(i,j, -coef, 0);
			}

			// Compute maximum flow
			real energy = graph->maxflow();

			// Extract solution
			x.resize(n);
			for (int i=0; i<n; ++i) {
				x[i] = graph->what_segment(i, Graph<real,real,real>::SOURCE);
			}

			return terms0 + energy;
		}


		char get_solution(int i)
		{
			return x.at(i);
		}



		//
		// Takes a list of pairs (i,j) and tries to find a solution where
		// x[i] != x[j]. The function is assumed to be symmetric, i.e. 
		// (1,1) is never a unique optimal solution.
		//
		void resolve_different(const std::vector<std::pair<int,int> >& pairs)
		{
			ASSERT(graph);
			Petter::resolve_different(*graph, x, pairs);
		}


	protected:

		inline real submodular_coef(real coef);

		int n;
		int new_vars;
		real terms0;
		std::map<int, real> terms1;
		std::map<int, real>  terms1_new;
		std::map< std::pair<int,int>, real> terms2;
		std::map< std::pair<int,int>, real> terms2_new;
		std::vector<char> x;
		Graph<real,real,real>* graph;
	};

	
	template<typename real>
	inline real Minimizer<real>::submodular_coef(real coef)
	{
		// Check submodularity
		ASSERT(coef <= 0);
		return coef;
	}

	template<>
	inline double Minimizer<double>::submodular_coef(double coef)
	{
		// Check submodularity
		ASSERT_STR(coef <= 1e-9, "Encountered non-submodular term when minimizing. This is probably due to floating-point arithmetic errors");
		return std::min(double(0),coef);
	}
}

#endif