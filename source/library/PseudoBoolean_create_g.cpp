//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Creates a submodular relaxation g(x,y) using linear programming. Clp is required.
//

#include <typeinfo>

#include "PseudoBoolean.h"

#include <coin/ClpSimplex.hpp>

namespace Petter
{
	//
	// The default behaviour is to raise an error message
	//
	template<typename real> template<typename orgreal>
	void SymmetricPseudoBoolean<real>::create_lp(const PseudoBoolean<orgreal>& pbf) 
	{
		auto& id = typeid(real);
		std::string msg = "Linear programming not available for type ";
		msg += id.name();
		throw std::runtime_error(msg);
	}

	//
	// For doubles, LP is available
	//
	template<> template<typename orgreal>
	void SymmetricPseudoBoolean<double>::create_lp(const PseudoBoolean<orgreal>& pbf) 
	{
		// Convenient to have this name available
		typedef double real;

		clear();
		constant = pbf.constant;

		
		int nLPVars = int( 2*pbf.aij.size() + 4*pbf.aijk.size() + 8*pbf.aijkl.size() );

		int nConstraints = int( pbf.aij.size() + pbf.aijk.size() + pbf.aijkl.size() // equality
		                        + 2*pbf.aij.size()                               ); // submodularity 

		// Compute the number of entries in the 
		// constraint matrix
		size_t nEntries =  2 * pbf.aij.size()   +  // equality
		                   4 * pbf.aijk.size()  +  // equality
						   8 * pbf.aijkl.size() +  // equality
						   2 * pbf.aij.size()   +  // submodularity
						   6 * pbf.aijk.size();    // submodularity
		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			real aijkl = itr->second;
			if (aijkl >= 0) {
				nEntries += 36;
			}
			else {
				nEntries += 39;
			}
		}

		//Description of sparse matrix
		vector<int> rows;
		vector<int> cols;
		vector<double> values;
		rows.reserve(nEntries);
		cols.reserve(nEntries);
		values.reserve(nEntries);

		//Other LP parameters
		vector<double> rhs_lower(nConstraints, 0.0);
		vector<double> rhs_upper(nConstraints, 0.0);
		double var_limit = 100000000;
		vector<double> var_lb(nLPVars, -var_limit); 
		vector<double> var_ub(nLPVars, var_limit); 
		vector<double> cost(nLPVars, 0.0);

		//Keeps track of the current row
		int con = 0;
		//Adds a value to the sparse matrix
		auto add_element = [&rows,&cols,&values](size_t row, size_t col, double value) 
		{
			rows.push_back(int(row));
			cols.push_back(int(col));
			values.push_back(value);
		};
		//Changes the right-hand side of the constraints
		auto change_rhs = [&rhs_lower, &rhs_upper](size_t row, double lower, double upper) 
		{
			rhs_lower.at(row) = lower;
			rhs_upper.at(row) = upper;
		};


		//////////////////////////
		// Equality constraints //
		//////////////////////////

		// bi = ai
		// (there is no freedom in picking these)
		constant = pbf.constant;
		for (auto itr=pbf.ai.begin(); itr != pbf.ai.end(); ++itr) {
			bi[itr->first] = static_cast<real>( itr->second );
		}

		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			int i=get_i(itr->first);
			int j=get_j(itr->first);
			real aij = itr->second;

			// bij + cij = aij
			add_element(con, ib(i,j), 1);
			add_element(con, ic(i,j), 1);
			change_rhs(con, aij, aij);
			con++;

			//Objective function
			cost[ib(i,j)] = -1;
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			int i=get_i(itr->first);
			int j=get_j(itr->first);
			int k=get_k(itr->first);
			real aijk = itr->second;

			// bijk + ... + eijk = aijk
			add_element(con, ib(i,j,k), 1);
			add_element(con, ic(i,j,k), 1);
			add_element(con, id(i,j,k), 1);
			add_element(con, ie(i,j,k), 1);
			change_rhs(con, aijk, aijk);
			con++;

			//Objective function
			cost[ib(i,j,k)] = -1;

			// Positivity constraints 
			if (aijk >= 0) {
				var_lb[ ib(i,j,k) ] = 0;
				var_lb[ ic(i,j,k) ] = 0;
				var_lb[ id(i,j,k) ] = 0;
				var_lb[ ie(i,j,k) ] = 0;
			}
			else {
				var_ub[ ib(i,j,k) ] = 0;
				var_ub[ ic(i,j,k) ] = 0;
				var_ub[ id(i,j,k) ] = 0;
				var_ub[ ie(i,j,k) ] = 0;
			}
		}

		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			int i=get_i(itr->first);
			int j=get_j(itr->first);
			int k=get_k(itr->first);
			int l=get_l(itr->first);
			real aijkl = itr->second;

			// bijkl + ... + sijkl = aijkl
			add_element(con, ib(i,j,k,l), 1);
			add_element(con, ic(i,j,k,l), 1);
			add_element(con, id(i,j,k,l), 1);
			add_element(con, ie(i,j,k,l), 1);
			add_element(con, ip(i,j,k,l), 1);
			add_element(con, iq(i,j,k,l), 1);
			add_element(con, ir(i,j,k,l), 1);
			add_element(con, is(i,j,k,l), 1);
			change_rhs(con, aijkl, aijkl);
			con++;

			//Objective function
			cost[ib(i,j,k,l)] = -1;

			// Positivity constraints 
			if (aijkl >= 0) {
				var_lb[ ib(i,j,k,l) ] = 0;
				var_lb[ ic(i,j,k,l) ] = 0;
				var_lb[ id(i,j,k,l) ] = 0;
				var_lb[ ie(i,j,k,l) ] = 0;
				var_lb[ ip(i,j,k,l) ] = 0;
				var_lb[ iq(i,j,k,l) ] = 0;
				var_lb[ ir(i,j,k,l) ] = 0;
				var_lb[ is(i,j,k,l) ] = 0;
			}
			else {
				var_ub[ ib(i,j,k,l) ] = 0;
				var_ub[ ic(i,j,k,l) ] = 0;
				var_ub[ id(i,j,k,l) ] = 0;
				var_ub[ ie(i,j,k,l) ] = 0;
				var_ub[ ip(i,j,k,l) ] = 0;
				var_ub[ iq(i,j,k,l) ] = 0;
				var_ub[ ir(i,j,k,l) ] = 0;
				var_ub[ is(i,j,k,l) ] = 0;
			}
		}

		///////////////////////////////
		// Submodularity constraints //
		///////////////////////////////
		//const double zero = 1e-7; //Set to e.g. 1e-6 for margin
		const double zero = 0;
		map<pair,int> ind_to_con;

		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			ind_to_con[itr->first] = con;
			add_element(con, ib(i,j), 1);
			change_rhs(con, -var_limit, -zero); // LHS <= 0
			con++;
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real aijk = itr->second;
			int c;

			if (aijk >= 0) {
				// ij*
				c = ind_to_con[make_pair(i,j)];
				add_element(c, ib(i,j,k), 1);
				add_element(c, ic(i,j,k), 1);

				// i*j
				c = ind_to_con[make_pair(i,k)];
				add_element(c, ib(i,j,k), 1);
				add_element(c, id(i,j,k), 1);

				// *ij
				c = ind_to_con[make_pair(j,k)];
				add_element(c, ib(i,j,k), 1);
				add_element(c, ie(i,j,k), 1);
			}
		}
		
		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real aijkl = itr->second;
			int a = aijkl >= 0 ? 1 : -1;
			int c;
			
			// ij**
			c = ind_to_con[make_pair(i,j)];
			if (aijkl >= 0) {
				add_element(c, ib(i,j,k,l), 1);
				add_element(c, iq(i,j,k,l), 1+a);
			}
			else {
				add_element(c, iq(i,j,k,l), a);
			}
			add_element(c, ic(i,j,k,l), a);
			add_element(c, id(i,j,k,l), a);

			// i*j*
			c = ind_to_con[make_pair(i,k)];
			if (aijkl >= 0) {
				add_element(c, ib(i,j,k,l), 1);
				add_element(c, ir(i,j,k,l), 1+a);
			}
			else {
				add_element(c, ir(i,j,k,l), a);
			}
			add_element(c, ic(i,j,k,l), a);
			add_element(c, ie(i,j,k,l), a);
			

			// i**j
			c = ind_to_con[make_pair(i,l)];
			if (aijkl >= 0) {
				add_element(c, ib(i,j,k,l), 1);
				add_element(c, is(i,j,k,l), 1+a);
			}
			else {
				add_element(c, is(i,j,k,l), a);
			}
			add_element(c, id(i,j,k,l), a);
			add_element(c, ie(i,j,k,l), a);
			

			// *ij*
			c = ind_to_con[make_pair(j,k)];
			if (aijkl >= 0) {
				add_element(c, ib(i,j,k,l), 1);
				add_element(c, is(i,j,k,l), 1);
			}
			add_element(c, ic(i,j,k,l), a);
			add_element(c, ip(i,j,k,l), a);

			// *i*j
			c = ind_to_con[make_pair(j,l)];
			if (aijkl >= 0) {
				add_element(c, ib(i,j,k,l), 1);
				add_element(c, ir(i,j,k,l), 1);
			}
			add_element(c, id(i,j,k,l), a);
			add_element(c, ip(i,j,k,l), a);

			// **ij
			c = ind_to_con[make_pair(k,l)];
			if (aijkl >= 0) {
				add_element(c, ib(i,j,k,l), 1);
				add_element(c, iq(i,j,k,l), 1);
			}	
			add_element(c, ie(i,j,k,l), a);
			add_element(c, ip(i,j,k,l), a);
		}


		ind_to_con.clear();
		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			ind_to_con[itr->first] = con;
			add_element(con, ic(i,j), 1);
			change_rhs(con, zero, var_limit); // LHS >= 0
			con++;
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real aijk = itr->second;
			int c;

			if (aijk < 0) {
				// ij*
				c = ind_to_con[make_pair(i,j)];
				add_element(c, id(i,j,k), 1);
				add_element(c, ie(i,j,k), 1);

				// i*j
				c = ind_to_con[make_pair(i,k)];
				add_element(c, ic(i,j,k), 1);
				add_element(c, ie(i,j,k), 1);

				// *ij
				c = ind_to_con[make_pair(j,k)];
				add_element(c, ic(i,j,k), 1);
				add_element(c, id(i,j,k), 1);
			}
		}

		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real aijkl = itr->second;
			int a = aijkl >= 0 ? +1 : -1;
			int c;

			// ij**
			c = ind_to_con[make_pair(i,j)];
			if (aijkl < 0) {
				add_element(c, ie(i,j,k,l), 1);
				add_element(c, ip(i,j,k,l), 1);
			}
			add_element(c, ir(i,j,k,l), -a);
			add_element(c, is(i,j,k,l), -a);

			// i*j*
			c = ind_to_con[make_pair(i,k)];
			if (aijkl < 0) {
				add_element(c, id(i,j,k,l), 1);
				add_element(c, ip(i,j,k,l), 1);
			}
			add_element(c, iq(i,j,k,l), -a);
			add_element(c, is(i,j,k,l), -a);

			// i**j
			c = ind_to_con[make_pair(i,l)];
			if (aijkl < 0) {
				add_element(c, ic(i,j,k,l), 1);
				add_element(c, ip(i,j,k,l), 1);
			}
			add_element(c, iq(i,j,k,l), -a);
			add_element(c, ir(i,j,k,l), -a);

			// *ij*
			c = ind_to_con[make_pair(j,k)];
			if (aijkl < 0) {
				add_element(c, id(i,j,k,l), 1);
				add_element(c, ie(i,j,k,l), 1);
			}
			add_element(c, iq(i,j,k,l), -a);
			add_element(c, ir(i,j,k,l), -a);

			// *i*j
			c = ind_to_con[make_pair(j,l)];
			if (aijkl < 0) {
				add_element(c, ic(i,j,k,l), 1);
				add_element(c, ie(i,j,k,l), 1);
			}
			add_element(c, iq(i,j,k,l), -a);
			add_element(c, is(i,j,k,l), -a);

			// **ij
			c = ind_to_con[make_pair(k,l)];
			if (aijkl < 0) {
				add_element(c, ic(i,j,k,l), 1);
				add_element(c, id(i,j,k,l), 1);
			}	
			add_element(c, ir(i,j,k,l), -a);
			add_element(c, is(i,j,k,l), -a);
		}

		ASSERT(nLPVars == this->nlpvars);

		// Make sure that the precomputation of space is correct
		ASSERT( rows.size() == nEntries );
		// Make sure that the precomputation of rows is correct
		ASSERT( con == nConstraints );

		CoinPackedMatrix coinMatrix(false,&rows[0],&cols[0],&values[0], CoinBigIndex(values.size()) );
		ClpSimplex lpSolver;
		lpSolver.loadProblem (coinMatrix, &var_lb[0], &var_ub[0], &cost[0], &rhs_lower[0], &rhs_upper[0]);
		lpSolver.setLogLevel(0);
		
		int error = lpSolver.dual(); //primal, dual, barrier
		if (error != 0) {
			// I think this is what the error codes mean
			if (error == 1) {
				throw std::runtime_error("The linear program is infeasible");
			}
			else if (error == 2) {
				throw std::runtime_error("The linear program is unbounded");
			}
			else if (error == 3) {
				throw std::runtime_error("Maximum number of iterations reached when solving the linear programming problem");
			}
			else {
				throw std::runtime_error("Unknown problem with the linear programming solver (Clp)");
			}
		}

		double* lpvars = lpSolver.primalColumnSolution();

		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			const pair& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);

			bij[ind] = lpvars[ ib(i,j) ];
			cij[ind] = lpvars[ ic(i,j) ];
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			const triple& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);

			bijk[ind] = lpvars[ ib(i,j,k) ];
			cijk[ind] = lpvars[ ic(i,j,k) ];
			dijk[ind] = lpvars[ id(i,j,k) ];
			eijk[ind] = lpvars[ ie(i,j,k) ];
		}
	
		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			const quad& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);
			int l=get_l(ind);

			bijkl[ind] = lpvars[ ib(i,j,k,l) ];
			cijkl[ind] = lpvars[ ic(i,j,k,l) ];
			dijkl[ind] = lpvars[ id(i,j,k,l) ];
			eijkl[ind] = lpvars[ ie(i,j,k,l) ];
			pijkl[ind] = lpvars[ ip(i,j,k,l) ];
			qijkl[ind] = lpvars[ iq(i,j,k,l) ];
			rijkl[ind] = lpvars[ ir(i,j,k,l) ];
			sijkl[ind] = lpvars[ is(i,j,k,l) ];
		}

	}


	//TODO: find a better way of doing this
	template void SymmetricPseudoBoolean<double>::create_lp(const PseudoBoolean<double>& pbf);
	template void SymmetricPseudoBoolean<double>::create_lp(const PseudoBoolean<int>& pbf);
	template void SymmetricPseudoBoolean<int>::create_lp(const PseudoBoolean<double>& pbf);
	template void SymmetricPseudoBoolean<int>::create_lp(const PseudoBoolean<int>& pbf);

}

#include "pb_instances.inc"
