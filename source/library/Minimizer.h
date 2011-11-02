//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Minimizes a submodular function of degree <= 3. Each term does not have to 
// be submodular, but the complete function has to be.
//

#ifndef PETTER_MINIMIZER_H
#define PETTER_MINIMIZER_H

#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>

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

	template<typename real>
	class Minimizer 
	{
	public:
		Minimizer(int nVars) : n(2*nVars), terms0(0), terms1(2*nVars, 0)
		{
			this->nVars = nVars;
			terms1.reserve(4*nVars);
		}

		void AddUnaryTerm(int i, real E0, real E1)
		{
			terms0 += E0;
			terms1.at(i) += E1-E0;
		}

		void AddPairwiseTerm(int i, int j, real E00, real E01, real E10, real E11)
		{
			if (j<i) {
				AddPairwiseTerm(j,i, E00, E10, E01, E11);
				return;
			}

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

			int l = nVars++;
			terms1.push_back(0);

			if (aijk < 0) { 
				// Use reduction
				// -xi*xj*xk = min{z} z*(2-xi-xj-xk) 
				terms2[ unordered_pair(i,l)] += aijk;
				terms2[ unordered_pair(j,l)] += aijk;
				terms2[ unordered_pair(k,l)] += aijk;
				terms1[ l ] += -2*aijk;
			}
			else if (aijk > 0) { 
				// Use reduction
				// xi*xj*xk = min{z} z*(1-xi-xj-xk)  + xi*xj + xi*xk + xj*xk
				terms2[ unordered_pair(i,l)] += -aijk;
				terms2[ unordered_pair(j,l)] += -aijk;
				terms2[ unordered_pair(k,l)] += -aijk;
				terms1.at(l) += aijk;
				terms2[ unordered_pair(i,j)] += aijk;
				terms2[ unordered_pair(i,k)] += aijk;
				terms2[ unordered_pair(j,k)] += aijk;
			}
		}


		real minimize()
		{
			using namespace std;

			Graph<real,real,real> graph(nVars,int(terms2.size()),error_function);
			graph.add_node(nVars);

			for (int i=0; i<nVars; ++i) {
				real c = terms1[i];
				//if (c!=0){ 
				//	cout << c << " * x[" << i << "]\n";
				//}
				// c*x(i)
				if (c >= 0) {
					graph.add_tweights(i, c, 0);
				}
				else {
					terms0 += c;
					graph.add_tweights(i, 0,  -c);
				}
			}

			for (auto itr = terms2.begin(); itr != terms2.end(); ++itr) {
				int i,j;
				std::tie(i,j) = itr->first;
				real coef     = itr->second;
				//if (coef!=0){ 
				//	cout << coef << " * x[" << i << "] * x[" << j << "]\n";
				//}
				// Check submodularity
				ASSERT(coef <= 1e-6);
				coef = min(real(0),coef);
				// Add to graph
				terms0 += coef;
				graph.add_tweights(j, 0,  -coef);
				graph.add_edge(i,j, -coef, 0);
			}

			// Compute maximum flow
			real energy = graph.maxflow();

			// Extract solution
			x.resize(n);
			for (int i=0; i<n; ++i) {
				x[i] = graph.what_segment(i);
			}

			//cout << "terms0 : " << terms0 << endl;
			//cout << "energy : " << energy << endl;

			return terms0 + energy;
		}

		char get_solution(int i)
		{
			return x.at(i);
		}

	protected:
		int n;
		int nVars;
		real terms0;
		std::vector<real> terms1;
		std::map< std::pair<int,int>, real> terms2;
		std::vector<char> x;
	};

}

#endif