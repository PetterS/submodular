//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Reduces a function f(x) using known persistencies.
//


#include "PseudoBoolean.h"

namespace Petter
{
	template<typename real>
	void PseudoBoolean<real>::reduce(const vector<label>& x)
	{
		//
		// Reduce degree 4 monomials
		//
		for (auto itr = aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real& a = itr->second;

			if (x.at(i) == 0) {
				a = 0;
			}
			else if (x.at(i) == 1) {
				aijk[make_triple(j,k,l)] += a;
				a = 0;
			}

			if (x.at(j) == 0) {
				a = 0;
			}
			else if (x.at(j) == 1) {
				aijk[make_triple(i,k,l)] += a;
				a = 0;
			}

			if (x.at(k) == 0) {
				a = 0;
			}
			else if (x.at(k) == 1) {
				aijk[make_triple(i,j,l)] += a;
				a = 0;
			}

			if (x.at(l) == 0) {
				a = 0;
			}
			else if (x.at(l) == 1) {
				aijk[make_triple(i,j,k)] += a;
				a = 0;
			}
		}

		//
		// Reduce degree 3 monomials
		//
		for (auto itr = aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real& a = itr->second;

			if (x.at(i) == 0) {
				a = 0;
			}
			else if (x.at(i) == 1) {
				aij[make_pair(j,k)] += a;
				a = 0;
			}

			if (x.at(j) == 0) {
				a = 0;
			}
			else if (x.at(j) == 1) {
				aij[make_pair(i,k)] += a;
				a = 0;
			}

			if (x.at(k) == 0) {
				a = 0;
			}
			else if (x.at(k) == 1) {
				aij[make_pair(i,j)] += a;
				a = 0;
			}
		}

		//
		// Reduce degree 2 monomials
		//
		for (auto itr = aij.begin(); itr != aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			real& a = itr->second;

			if (x.at(i) == 0) {
				a = 0;
			}
			else if (x.at(i) == 1) {
				ai[j] += a;
				a = 0;
			}

			if (x.at(j) == 0) {
				a = 0;
			}
			else if (x.at(j) == 1) {
				ai[i] += a;
				a = 0;
			}
		}

		//
		// Reduce degree 1 monomials
		//
		for (auto itr = ai.begin(); itr != ai.end(); ++itr) {
			int i  = itr->first;
			real& a = itr->second;

			if (x.at(i) == 0) {
				a = 0;
			}
			else if (x.at(i) == 1) {
				constant += a;
				a = 0;
			}
		}


		// Remove unused monomials; care has to be taken 
		// so that lower degree ones needed by higher degrees 
		// are not removed.

		// Go through degree-1 monomials
		auto itr = ai.begin();
		while (itr != ai.end()) {
			bool should_remove = false;

			if (itr->second == 0) {
				should_remove = true;
			}

			if (should_remove) {
				itr = ai.erase(itr);
			}
			else {
				itr++;
			}
		}

		// Go through degree-4 monomials
		map<triple,bool> used3;
		{		
			auto itr = aijkl.begin();
			while (itr != aijkl.end()) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				int k = get_k(itr->first);
				int l = get_l(itr->first);
				bool should_remove = false;

				if (itr->second == 0) {
					should_remove = true;
				}

				auto prev = itr++;
				if (should_remove) {
					aijkl.erase(prev);
				}
				else {
					// We still need all lower-degre monomials
					used3[make_triple(i,j,k)] = true;
					used3[make_triple(i,j,l)] = true;
					used3[make_triple(i,k,l)] = true;
					used3[make_triple(j,k,l)] = true;
				}
			}
		}

		// Go through degree-3 monomials
		map<pair,bool> used2;
		{
			auto itr = aijk.begin();
			while (itr != aijk.end()) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				int k = get_k(itr->first);
				bool should_remove = false;

				if (itr->second == 0 && !used3[itr->first] ) {
					should_remove = true;
				}

				auto prev = itr++;
				if (should_remove) {
					aijk.erase(prev);
				}
				else {
					// We still need all lower-degre monomials
					used2[make_pair(i,j)] = true;
					used2[make_pair(i,k)] = true;
					used2[make_pair(j,k)] = true;
				}
			}
		}


		// Go through degree-2 monomials
		{
			auto itr = aij.begin();
			while (itr != aij.end()) {
				bool should_remove = false;

				if (itr->second == 0 && !used2[itr->first] ) {
					should_remove = true;
				}

				auto prev = itr++;
				if (should_remove) {
					aij.erase(prev);
				}
			}
		}



		///////////
		// Check //
		///////////
		for (auto itr = aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			ASSERT( aij.find(make_pair(i,j)) != aij.end() );
			ASSERT( aij.find(make_pair(i,k)) != aij.end() );
			ASSERT( aij.find(make_pair(i,l)) != aij.end() );
			ASSERT( aij.find(make_pair(j,k)) != aij.end() );
			ASSERT( aij.find(make_pair(j,l)) != aij.end() );
			ASSERT( aij.find(make_pair(k,l)) != aij.end() );
		}
		for (auto itr = aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			ASSERT( aij.find(make_pair(i,j)) != aij.end() );
			ASSERT( aij.find(make_pair(i,k)) != aij.end() );
			ASSERT( aij.find(make_pair(j,k)) != aij.end() );
		}

	} //reduce


} //namespace Petter

#include "pb_instances.inc"
