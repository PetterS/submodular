//
// Petter Strandmark 2011
// petter@maths.lth.se
//
#ifndef PETTER_POSIFORM_H
#define PETTER_POSIFORM_H


#include <vector>
#include <map>
#include <stdexcept>
#include <iostream>
#include <tuple>

#include <sstream>

#include "PseudoBoolean.h" 
#include "VertexPacking.h"

namespace {
	mt19937 posi_engine(342); // Deterministic for testing purposes
}


namespace Petter
{
	using std::map;
	using std::vector;

	template<int n> 
	struct Ind
	{
		Ind() 
		{
			for (int i=0;i<n;++i){
				ind[i] = -1;
				bar[i] = false;
			}
		}

		bool operator<(const Ind<n>& rhs) const
		{
			// First compare indices
			for (int i=0;i<n;++i) {
				if (ind[i] < rhs.ind[i]) {
					return true;
				}
				if (ind[i] > rhs.ind[i]) {
					return false;
				}
			}

			// Then compare bars
			for (int i=0;i<n;++i) {
				if (bar[i] && !rhs.bar[i]) {
					return true;
				}
				if (!bar[i] && rhs.bar[i]) {
					return false;
				}
			}
			
			return false;
		}

		bool operator==(const Ind<n>& rhs) const
		{
			// First compare indices
			for (int i=0;i<n;++i) {
				if (ind[i] != rhs.ind[i]) {
					return false;
				}
			}

			// Then compare bars
			for (int i=0;i<n;++i) {
				if (bar[i] != rhs.bar[i]) {
					return false;
				}
			}
			
			return true;
		}

		int ind[n];
		bool bar[n];
	};


	template<typename real, int degree>
	class Posiform
	{
	public:
		typedef Ind<degree> I;

		Posiform();
		Posiform(const PseudoBoolean<real>& f, bool randomize=false);

		void print();

		// Adding monomials
		void add_monomial(int i, real a, bool randomize=false);
		void add_monomial(int i, int j, real a, bool randomize=false);
		void add_monomial(int i, int j, int k, real a, bool randomize=false);
		void add_monomial(int i, int j, int k, int l, real a, bool randomize=false);

		int nvars() const;

		// Evaluate the function
		real eval(const vector<label>& x) const;

		// Minimize the function
		real maximize(vector<label>& x) const;

		real constant;
		map<I,real> terms;
	protected:
	private:
	};




	template<typename real, int degree>
	Posiform<real,degree>::Posiform()
	{
		constant = 0;
	}

	template<typename real, int degree>
	Posiform<real,degree>::Posiform(const PseudoBoolean<real>& f, bool randomize)
	{
		ASSERT(degree >= 4);

		//TODO: this function now converts -f instead of f.
		//      should be documented

		constant = -f.constant;

		for (auto itr=f.ai.begin(); itr!=f.ai.end(); ++itr) {
			int i = itr->first;
			real a = -itr->second;
			add_monomial(i,a,randomize);
		}

		for (auto itr=f.aij.begin(); itr!=f.aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			real a = -itr->second;
			add_monomial(i,j,a,randomize);
		}

		for (auto itr=f.aijk.begin(); itr!=f.aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real a = -itr->second;
			add_monomial(i,j,k,a,randomize);
		}

		for (auto itr=f.aijkl.begin(); itr!=f.aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real a = -itr->second;
			add_monomial(i,j,k,l,a,randomize);
		}
	}


	template<typename real, int degree>
	void Posiform<real,degree>::print()
	{
		cout << constant << " + ";
		for (auto itr = terms.begin(); itr != terms.end(); ++itr) {
			const int* ind = itr->first.ind;
			const bool* bar = itr->first.bar;
			real alpha = itr->second;
			ASSERT(alpha >= 0);

			cout << alpha;
			for (int d=0;d<degree;++d) {
				int i = ind[d];
				if (i>=0) {
					if (!bar[d]) {
						cout << "x"<<i;
					}
					else {
						cout << "~x"<<i;
					}
				}
			}

			auto tmp = itr;
			tmp++;
			if (tmp != terms.end()) {
				cout << " + ";
			}
		}
		cout << endl;
	}


	template<typename real, int degree>
	void Posiform<real,degree>::add_monomial(int i, real a, bool randomize)
	{
		ASSERT(degree>=1);
		I ind;
		ind.ind[0] = i;

		if (a > 0) {
			terms[ind] += a;
		}
		else if (a < 0) {
			ind.bar[0] = true;

			terms[ind] += -a;
			constant += a;
		}
	}

	template<typename real, int degree>
	void Posiform<real,degree>::add_monomial(int i, int j, real a, bool randomize)
	{
		ASSERT(degree>=2);
		I ind;
		ind.ind[0] = i;
		ind.ind[1] = j;


		int flip = 0;
		static auto random_ind = bind(uniform_int_distribution<int>(0,1), posi_engine);
		if (randomize) flip = random_ind();

		if (a > 0) {
			terms[ind] += a;
		}
		else if (a < 0) {
			ind.bar[flip] = true;

			terms[ind] += -a;
			if (flip == 0) {
				add_monomial(j,a);
			}
			else {
				add_monomial(i,a);
			}
		}
	}

	template<typename real, int degree>
	void Posiform<real,degree>::add_monomial(int i, int j, int k, real a, bool randomize)
	{
		ASSERT(degree>=3);
		I ind;
		ind.ind[0] = i;
		ind.ind[1] = j;
		ind.ind[2] = k;

		int flip = 0;
		static auto random_ind = bind(uniform_int_distribution<int>(0,2), posi_engine);
		if (randomize) flip = random_ind();

		if (a > 0) {
			terms[ind] += a;
		}
		else if (a < 0) {
			ind.bar[flip] = true;

			terms[ind] += -a;

			if (flip == 0) {
				add_monomial(j,k,a);
			}
			else if (flip == 1) {
				add_monomial(i,k,a);
			}
			else {
				add_monomial(i,j,a);
			}
		}
	}

	template<typename real, int degree>
	void Posiform<real,degree>::add_monomial(int i, int j, int k, int l, real a, bool randomize)
	{
		ASSERT(degree>=4);
		I ind;
		ind.ind[0] = i;
		ind.ind[1] = j;
		ind.ind[2] = k;
		ind.ind[3] = l;

		int flip = 0;
		static auto random_ind = bind(uniform_int_distribution<int>(0,3), posi_engine);
		if (randomize) flip = random_ind();

		if (a > 0) {
			terms[ind] += a;
		}
		else if (a < 0) {
			ind.bar[flip] = true;

			terms[ind] += -a;

			if (flip == 0) {
				add_monomial(j,k,l,a);
			}
			else if (flip == 1) {
				add_monomial(i,k,l,a);
			}
			else if (flip == 2) {
				add_monomial(i,j,l,a);
			}
			else {
				add_monomial(i,j,k,a);
			}
		}
	}


	template<typename real, int degree>
	real Posiform<real,degree>::eval(const vector<label>& x) const
	{
		real value = constant;

		for (auto itr = terms.begin(); itr != terms.end(); ++itr) {
			const int* ind = itr->first.ind;
			const bool* bar = itr->first.bar;
			real alpha = itr->second;
			ASSERT(alpha >= 0);

			bool equal_to_one = true;
			for (int d=0;d<degree;++d) {
				int i = ind[d];
				if (i>=0) {
					if ( (x.at(i) == 0 && !bar[d])  ||
					     (x.at(i) == 1 && bar[d]) ) {
						equal_to_one = false;
						break;
					}
				}
			}

			if (equal_to_one) {
				value += alpha;
			}
		}



		return value;
	}


	template<typename real, int degree>
	int Posiform<real,degree>::nvars() const
	{
		int n = 0;
		for (auto itr = terms.begin(); itr != terms.end(); ++itr) {
			const int* ind = itr->first.ind;
			for (int d=0;d<degree;++d) {
				int i = ind[d];
				if (i>n) {
					n=i;
				}
			}
		}
		return n+1;
	}

	template<typename real, int degree>
	real Posiform<real,degree>::maximize(vector<label>& x) const
	{
		// * Find out the number of variables
		// * Find out the number of higher-order terms  
		// * Calculate a suitable M
		// * Find out which variables are actually used
		int n = -1;
		int m = 0;
		real M = 0;
		map<int,bool> var_used;
		for (auto itr = terms.begin(); itr != terms.end(); ++itr) {
			const int* ind = itr->first.ind;
			real alpha = itr->second;

			M += alpha<0 ? -alpha : alpha;

			if (ind[1] >= 0) {
				// Higher-degree term
				m++;
			}

			for (int d=0;d<degree;++d) {
				int i = ind[d];
				var_used[i] = true;
				if (i>n) {
					n=i;
				}
			}
		}
		n++;

		M *= 2;
		ASSERT(M>0);
		ASSERT(x.size()>=size_t(n));

		// Create vertex packing problem
		VertexPacking<real> packing(2*n + m);
		real packingconst = 0;

		// * Add weights
		// * Add edges
		int y = 0;
		for (auto itr = terms.begin(); itr != terms.end(); ++itr) {
			const int* ind = itr->first.ind;
			const bool* bar = itr->first.bar;
			real alpha = itr->second;

			if (ind[1] >= 0) {
				// Higher-degree term
				packing.add_weight(2*n+y, alpha);

				for (int d=0;d<degree;++d) {
					int i = ind[d];
					if (i<0) {
						break;
					}
					if (bar[d]) {
						packing.add_edge(2*n+y, i);
					}
					else {
						packing.add_edge(2*n+y, n+i);
					}
				}

				y++;
			}
			else {
				// Linear term
				int i = ind[0];

				if (!bar[0]) {
					packing.add_weight(i, alpha);
				}
				else {
					packing.add_weight(n+i, alpha);
				}			
			}
		}

		for (int i=0;i<n;++i) {
			packing.add_edge(i, n+i);
			packing.add_weight(i  , M);
			packing.add_weight(n+i, M);
			packingconst -= M;
		}

		// Assert that we have added the correct number of vars
		ASSERT(y == m);

		vector<signed char> sol;
		//real energy = packing.solve_lp(sol);
		real energy = packing.solve(sol);
		energy += constant + packingconst;

		for (int i=0;i<n;++i) {
			if (var_used[i] || x.at(i)<0) {
				x.at(i) = sol.at(i);
			}
		}

		//cout << " x = [";
		//for (size_t i=0;i<n;++i) {
		//	if (sol.at(i) >= 0) {
		//		cout << int(sol.at(i));
		//	}
		//	else {
		//		cout << "-";
		//	}
		//}
		//cout << "]\n";

		//cout << "~x = [";
		//for (size_t i=n;i<2*n;++i) {
		//	if (sol.at(i) >= 0) {
		//		cout << int(sol.at(i));
		//	}
		//	else {
		//		cout << "-";
		//	}
		//}
		//cout << "]\n";

		//cout << " y = [";
		//for (size_t i=2*n;i<2*n+m;++i) {
		//	if (sol.at(i) >= 0) {
		//		cout << int(sol.at(i));
		//	}
		//	else {
		//		cout << "-";
		//	}
		//}
		//cout << "]" << endl;

		return energy;
	}

}
#endif
