//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Functions to minimize a general f(x) with an LP relaxation
//
// This is work in progress. You can leave this file out while compiling.
//

#include <iomanip>

#include "PseudoBoolean.h"

#include "QPBO.h"
#include "HOCR.h"
#include "graph.h"

#include <coin/ClpSimplex.hpp>

namespace Petter
{

	template<typename real>
	real PseudoBoolean<real>::minimize_lp(bool verbose) const
	{
		using namespace std;

		// TODO
		int n = nvars();
		int nLPVars = n + aij.size() + aijk.size() + aijkl.size();
		size_t nConstraints = 3*(aij.size() + aijk.size() + aijkl.size());
		size_t nEntries = 100000;

		//Description of sparse matrix
		vector<int> rows;
		vector<int> cols;
		vector<double> values;
		rows.reserve(nEntries);
		cols.reserve(nEntries);
		values.reserve(nEntries);

		//Other LP parameters
		double var_limit = 100000000;
		vector<double> var_lb(nLPVars, -var_limit); 
		vector<double> var_ub(nLPVars, var_limit); 
		vector<double> rhs_lower(nConstraints,  -var_limit);
		vector<double> rhs_upper(nConstraints, var_limit);
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


		// Linear terms
		for (auto itr=ai.begin(); itr != ai.end(); ++itr) {
			int i = itr->first;
			real a = itr->second;
			
			var_lb[i] = 0;
			var_ub[i] = 1;
			cost[i] = a;
		}

		//Keeps track of the current var
		int var = n;

		// Quadratic terms
		for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			real a = itr->second;
			int u = var++;
			var_lb[u] = 0;
			var_ub[u] = 1;
			if (a < 0) {
				// u - xi <= 0
				add_element(con, u, 1);
				add_element(con, i, -1);
				rhs_upper.at(con) = 0;
				con++;
				// u - xj <= 0
				add_element(con, u, 1);
				add_element(con, j, -1);
				rhs_upper.at(con) = 0;
				con++;
			}
			else if (a > 0) {
				// -1 <= u - xi - xj 
				add_element(con, u, 1);
				add_element(con, i, -1);
				add_element(con, j, -1);
				rhs_lower.at(con) = -1;
				con++;
			}
			cost[u] = a;
		}

		// Cubic terms
		for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real a = itr->second;
			int u = var++;
			var_lb[u] = 0;
			var_ub[u] = 1;
			if (a < 0) {
				// u - xi <= 0
				add_element(con, u, 1);
				add_element(con, i, -1);
				rhs_upper.at(con) = 0;
				con++;
				// u - xj <= 0
				add_element(con, u, 1);
				add_element(con, j, -1);
				rhs_upper.at(con) = 0;
				con++;
				// u - xk <= 0
				add_element(con, u, 1);
				add_element(con, k, -1);
				rhs_upper.at(con) = 0;
				con++;
			}
			else if (a > 0) {
				// -2 <= u - xi - xj - xk
				add_element(con, u, 1);
				add_element(con, i, -1);
				add_element(con, j, -1);
				add_element(con, k, -1);
				rhs_lower.at(con) = -2;
				con++;
			}
			cost[u] = a;
		}

		// Quartic terms
		for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real a = itr->second;
			int u = var++;
			var_lb[u] = 0;
			var_ub[u] = 1;
			if (a < 0) {
				// u - xi <= 0
				add_element(con, u, 1);
				add_element(con, i, -1);
				rhs_upper.at(con) = 0;
				con++;
				// u - xj <= 0
				add_element(con, u, 1);
				add_element(con, j, -1);
				rhs_upper.at(con) = 0;
				con++;
				// u - xk <= 0
				add_element(con, u, 1);
				add_element(con, k, -1);
				rhs_upper.at(con) = 0;
				con++;
				// u - xl <= 0
				add_element(con, u, 1);
				add_element(con, l, -1);
				rhs_upper.at(con) = 0;
				con++;
			}
			else if (a > 0) {
				// -3 <= u - xi - xj - xk - xl
				add_element(con, u, 1);
				add_element(con, i, -1);
				add_element(con, j, -1);
				add_element(con, k, -1);
				add_element(con, l, -1);
				rhs_lower.at(con) = -3;
				con++;
			}
			cost[u] = a;
		}


		ASSERT(nLPVars == var);

		CoinPackedMatrix coinMatrix(false,&rows[0],&cols[0],&values[0], CoinBigIndex(values.size()) );
		ClpSimplex lpSolver;
		lpSolver.loadProblem (coinMatrix, &var_lb[0], &var_ub[0], &cost[0], &rhs_lower[0], &rhs_upper[0]);
		lpSolver.setLogLevel(0);
		
		int error = lpSolver.dual(); //primal, dual, barrier
		if (error != 0) {
			throw std::runtime_error("Clp failed");
		}

		double* lpvars = lpSolver.primalColumnSolution();

		if (verbose) {

			const double th = 0.999;

			// Linear terms
			cout << "xi    : ";
			for (auto itr=ai.begin(); itr != ai.end(); ++itr) {
				int i = itr->first;
				cout  << lpvars[i] << ' ';
			}
			cout << endl;
			var = n;

			// Quadratic terms
			cout << "yij   : ";
			for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				real a = itr->second;
				int u = var++;

				if (a>0) {
					if ( (lpvars[u] > th && (lpvars[i]<(1-th) || lpvars[j]<(1-th) )) ||
						 (lpvars[u] < (1-th) && (lpvars[i]>th && lpvars[j]>th )) ){
						stringstream sout;
						sout <<"Monomial x"<<i<<"*x"<<j<<" is inconsistent.";
						throw runtime_error(sout.str().c_str());
					}
				}

				cout << lpvars[u] << ' ';
			}
			cout << endl;

			// Cubic terms
			cout << "yijk  : ";
			for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				int k = get_k(itr->first);
				real a = itr->second;
				int u = var++;

				if (a>0) {
					if ( (lpvars[u] > th && (lpvars[i]<(1-th) || lpvars[j]<(1-th) || lpvars[k]<(1-th))) ||
						 (lpvars[u] < (1-th) && (lpvars[i]>th && lpvars[j]>th && lpvars[k]>th)) ){
						stringstream sout;
						sout <<"Monomial x"<<i<<"*x"<<j<<"*x" <<k << " is inconsistent.";
						throw runtime_error(sout.str().c_str());
					}
				}

				cout << lpvars[u] << ' ';
			}
			cout << endl;

			// Quartic terms
			cout << "yijkl : ";
			for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				int k = get_k(itr->first);
				int l = get_l(itr->first);
				real a = itr->second;
				int u = var++;

				if (a>0) {
					if ( (lpvars[u] > th && (lpvars[i]<(1-th) || lpvars[j]<(1-th) || lpvars[k]<(1-th) || lpvars[l]<(1-th))) ||
						 (lpvars[u] < (1-th) && (lpvars[i]>th && lpvars[j]>th && lpvars[k]>th && lpvars[l]>th)) ){
						stringstream sout;
						sout <<"Monomial x"<<i<<"*x"<<j<<"*x"<<k<<"*x"<<l<<" is inconsistent.";
						throw runtime_error(sout.str().c_str());
					}
				}

				cout << lpvars[u] << ' ';
			}
			cout << endl;
		}

		return constant + lpSolver.objectiveValue();
	}


}



#include "pb_instances.inc"
