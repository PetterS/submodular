//
// Petter Strandmark 2011
// petter@maths.lth.se
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include <coin/ClpSimplex.hpp>

#include "VertexPacking.h"
#include "Petter-Color.h"

#ifdef WIN32
    //Choose whether to use IBFS
    //#define USE_IBFS
#endif 

#ifdef USE_IBFS
	#include "ibfs.h"
#else
	#include "graph.h"
#endif

using namespace std;

void err_function(char* err) 
{
	throw runtime_error(err);
}

namespace Petter 
{


	template<typename real>
	VertexPacking<real>::VertexPacking(size_t n) :
		w(n,0)
	{
	}

	template<typename real>
	void VertexPacking<real>::set_weight(size_t i, real weight)
	{
		ASSERT(weight>=0);
		w.at(i)=weight;
	}

	template<typename real>
	void VertexPacking<real>::add_weight(size_t i, real weight)
	{
		ASSERT(weight>=0);
		w.at(i)+=weight;
	}

	template<typename real>
	void VertexPacking<real>::add_edge(size_t i, size_t j)
	{
		ASSERT(i<w.size());
		ASSERT(j<w.size());
		if (i==j) {
			return;
		}
		if (j<i) {
			add_edge(j,i);
			return;
		}
		E.push_back( std::make_pair(i,j) );
	}

	template<typename real>
	real VertexPacking<real>::solve_slower(std::vector<signed char>& x)
	{
		size_t n = w.size();

		//Encode as a bipartite vertex packing problem
		BipartiteVertexPacking<real> packing(n,n);

		for (size_t i=0;i<n;++i) {
			packing.set_left_weight(i, w[i]);
			packing.set_right_weight(i, w[i]);
		}

		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			packing.add_edge(itr->first, itr->second);
			packing.add_edge(itr->second, itr->first);
		}

		std::vector<signed char> xm(n),xn(n);
		//packing.solve_exhaustive(xm,xn);
		packing.solve(xm,xn);

		x.resize(n);
		real energy = 0;

		for (size_t i=0;i<n;++i) {
			if (xm[i]==xn[i]) {
				x[i] = xm[i];
				energy += x[i]*w[i];
			}
			else {
				x[i] = -1;
				energy += w[i]/real(2);
			}
		}

		return energy;
	}

	template<typename real>
	real VertexPacking<real>::solve(std::vector<signed char>& x)
	{
		//
		// Solves the LP relaxation of the vertex packing as a bipartite vertex 
		// packing. However, the bipartite problem is not create explicitly; instead,
		// the maximum flow graph is create directly.
		//
		// See solve_slow for details of what the bipartite problem looks like
		//
		size_t n = w.size();

		#ifdef USE_IBFS
			IBFSGraph<real,real,real> graph( int(m+n) , int(E.size()), err_function);
		#else
			Graph<real,real,real> graph( int(n+n) , int(2*E.size()), err_function);
		#endif

		graph.add_node( int(n+n) );

		real inf = 0;
		for (size_t i=0;i<n;++i) {
			graph.add_tweights( int(i), 0, w[i]);
			graph.add_tweights( int(n+i), w[i], 0);
			inf += w[i];
		}

		inf *= 4;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			graph.add_edge( int(n+itr->second), int(itr->first), inf, 0);
			graph.add_edge( int(n+itr->first), int(itr->second), inf, 0);
		}

		graph.maxflow();

		x.resize(n);
		real energy = 0;

		for (size_t i=0;i<n;++i) {
			signed char x1 = graph.what_segment(i);
			signed char x2 = graph.what_segment(i+n); 
			if (x1!=x2) {
				x[i] = x1;
				energy += x1*w[i];
			}
			else {
				x[i] = -1;
				energy += w[i]/real(2);
			}
		}

		return energy;
	}

	template<typename real>
	real VertexPacking<real>::solve_exhaustive(std::vector<signed char>& xsol)
	{
		size_t n = w.size();
		ASSERT(n<=30); // Otherwise too big

		vector<signed char> x(n,0);
		real best_energy = 0;
		while (true) {
			x[0]++;
			size_t i=0;
			while (x[i]>1) {
				x[i]=0;
				i++;
				if (i==n) {
					break;
				}
				x[i]+=1;
			}
			if (i==n) {
				break;
			}
			
			real energy = -1;
			if (is_feasible(x)) {
				energy = 0;
				for (size_t i=0; i<n; ++i) {
					energy += x[i]*w[i];
				}
			}

			if (energy > best_energy) {
				best_energy = energy;
				xsol = x;
			}
		}

		for (size_t i=0;i<n;++i) { x[i] = 0; }

		x[0]=-1;
		while (true) {
			x[0]++;
			size_t i=0;
			while (x[i]>1) {
				x[i]=0;
				i++;
				if (i==n) {
					break;
				}
				x[i]+=1;
			}
			if (i==n) {
				break;
			}
			
			real energy = -1;
			if (is_feasible(x)) {
				energy = 0;
				for (size_t i=0; i<n; ++i) {
					energy += x[i]*w[i];
				}
			}

			if (energy == best_energy) {

				cout << "Global maximum (";
				for (i=0;i<n;++i) {
					cout << int(x[i]);
				}
				cout << ") = ";
				cout << WHITE << energy << NORMAL << endl;
			}

		}
		return best_energy;
	}

	template<typename real>
	real VertexPacking<real>::solve_lp(std::vector<signed char>& x)
	{
		size_t n = w.size();
		size_t nConstraints = E.size() + n;
		size_t nEntries = 2*nConstraints;

		//Description of sparse matrix
		vector<int> rows;
		vector<int> cols;
		vector<double> values;
		rows.reserve(nEntries);
		cols.reserve(nEntries);
		values.reserve(nEntries);

		//Other LP parameters
		vector<double> var_lb(n, 0); 
		vector<double> var_ub(n, 1); 
		vector<double> rhs_lower(nConstraints,  0);
		vector<double> rhs_upper(nConstraints,  1);
		vector<double> cost(n, 0.0);

		//Keeps track of the current row
		int con = 0;
		//Adds a value to the sparse matrix
		auto add_element = [&rows,&cols,&values](size_t row, size_t col, double value) 
		{
			rows.push_back(int(row));
			cols.push_back(int(col));
			values.push_back(value);
		};
	

		for (size_t i=0;i<n;++i) {
			cost[i] = w[i];

			add_element(con, i, 1);
			con++;
		}

		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			add_element(con, itr->first, 1);
			add_element(con, itr->second, 1);
			con++;
		}

		CoinPackedMatrix coinMatrix(false,&rows[0],&cols[0],&values[0], CoinBigIndex(values.size()) );
		ClpSimplex lpSolver;
		lpSolver.loadProblem (coinMatrix, &var_lb[0], &var_ub[0], &cost[0], &rhs_lower[0], &rhs_upper[0]);
		lpSolver.setLogLevel(0);
		lpSolver.setOptimizationDirection(-1);

		int error = lpSolver.dual(); //primal, dual, barrier
		if (error != 0) {
			throw std::runtime_error("Clp failed");
		}

		double* lpvars = lpSolver.primalColumnSolution();

		const double th = 0.99999;

		x.resize(n);

		double energy = 0.0;

		for (size_t i=0;i<n;++i) {
			if (lpvars[i]<1-th) {
				x[i] = 0;
			}
			else if (lpvars[i]>th) {
				x[i] = 1;
			}
			else {
				x[i] = -1;
			}
			energy+=cost[i]*lpvars[i];
			//cout << lpvars[i] << ' ';
		}
		//cout << endl;

		return real(energy);
	}

	template<typename real>
	bool VertexPacking<real>::is_feasible(const std::vector<signed char>& x) const
	{
		bool feasible = true;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			if (x.at(itr->first)!=0 && x.at(itr->second)!=0 && itr->first != itr->second) {
				feasible = false;
				break;
			}
		}
		return feasible;
	}


	template<typename real>
	void VertexPacking<real>::load_from_file(std::string file)
	{
		ifstream fin(file.c_str());
		ASSERT(fin);

		size_t n;
		fin >> n;
		ASSERT(fin);
		w.resize(n);
		E.clear();

		for (size_t i=0;i<n;++i) {
			real val;
			fin >> val;
			ASSERT(fin);
			set_weight(i, val);
		}

		size_t m;
		fin >> m;
		ASSERT(fin);

		for (size_t i=0;i<m;++i) {
			int a,b;
			fin >> a >> b;
			ASSERT(fin);
			add_edge(a,b);
		}

	}

	template<typename real>
	void VertexPacking<real>::save_to_file(std::string file)
	{
		ofstream fout(file.c_str());

		fout << w.size() << endl;
		for (size_t i=0;i<w.size();++i) {
			fout << w[i] << ' ';
		}
		fout << endl;
		fout << E.size() << endl;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			fout << itr->first << ' ' << itr->second << endl;
		}
	}

	template<typename real>
	void VertexPacking<real>::print()
	{
		using namespace std;
		size_t n = w.size();
		for (size_t i=0;i<n;++i) {
			cout << setw(5) << w[i] << endl;
		}
		cout << endl;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			cout << itr->first << " <--> " << itr->second << endl;
		}
	}




	template<typename real>
	BipartiteVertexPacking<real>::BipartiteVertexPacking(size_t m, size_t n) :
		wn(n,0),
		wm(m,0)
	{

	}

	template<typename real>
	void BipartiteVertexPacking<real>::set_left_weight(size_t i, real w)
	{
		ASSERT(w>=0);
		wm.at(i)=w;
	}

	template<typename real>
	void BipartiteVertexPacking<real>::set_right_weight(size_t i, real w)
	{
		ASSERT(w>=0);
		wn.at(i)=w;
	}

	template<typename real>
	void BipartiteVertexPacking<real>::add_edge(size_t i, size_t j)
	{
		ASSERT(i<wm.size());
		ASSERT(j<wn.size());
		E.push_back( std::make_pair(i,j) );
	}

	template<typename real>
	real BipartiteVertexPacking<real>::solve(vector<signed char>& xm, vector<signed char>& xn)
	{
		size_t m = wm.size();
		size_t n = wn.size();

		#ifdef USE_IBFS
			IBFSGraph<real,real,real> graph( int(m+n) , int(E.size()), err_function);
		#else
			Graph<real,real,real> graph( int(m+n) , int(E.size()), err_function);
		#endif

		graph.add_node( int(m+n) );

		real inf = 0;
		for (size_t i=0;i<m;++i) {
			graph.add_tweights( int(i), 0, wm[i]);
			inf += wm[i];
		}

		for (size_t i=0;i<n;++i) {
			graph.add_tweights( int(m+i) , wn[i], 0);
			inf += wn[i];
		}

		inf *= 2;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			// x[i]=1 => x[j]=0
			graph.add_edge( int(itr->first), int(m+itr->second), 0, inf);
		}

		graph.maxflow();

		xm.resize(m);
		xn.resize(n);
		real w = 0;

		for (size_t i=0;i<m;++i) {
			xm[i] = graph.what_segment(i);
			w += xm.at(i)*wm.at(i);
		}

		for (size_t i=0;i<n;++i) {
			xn[i] = 1 - graph.what_segment( int(m+i) );
			w += xn.at(i)*wn.at(i);
		}

		return w;
	}

	template<typename real>
	real BipartiteVertexPacking<real>::solve_exhaustive(vector<signed char>& xmsol, vector<signed char>& xnsol)
	{
		size_t m = wm.size();
		size_t n = wn.size();
		ASSERT(n+m<=30); // Otherwise too big

		vector<signed char> x(n+m,0);
		vector<signed char> xm(m,0),xn(n,0);
		real best_energy = energy(xm,xn);
		while (true) {
			x[0]++;
			size_t i=0;
			while (x[i]>1) {
				x[i]=0;
				i++;
				if (i==n+m) {
					break;
				}
				x[i]+=1;
			}
			if (i==n+m) {
				break;
			}
			
			for (i=0;i<m;++i) {
				xm[i] = x[i];
			}
			for (i=0;i<n;++i) {
				xn[i] = x[m+i];
			}

			real en = energy(xm,xn);

			if (en > best_energy) {
				best_energy = en;
				xmsol = xm;
				xnsol = xn;
			}
		}

		//for (int i=0;i<n;++i) { x[i] = 0; }

		//x[0]=-1;
		//while (true) {
		//	x[0]++;
		//	int i=0;
		//	while (x[i]>1) {
		//		x[i]=0;
		//		i++;
		//		if (i==m+n) {
		//			break;
		//		}
		//		x[i]+=1;
		//	}
		//	if (i==n) {
		//		break;
		//	}
		//	
		//	for (i=0;i<m;++i) {
		//		xm[i] = x[i];
		//	}
		//	for (i=0;i<n;++i) {
		//		xn[i] = x[i];
		//	}

		//	real en = energy(xm,xn);

		//	if (en == best_energy) {

		//		cout << "Global maximum (";
		//		for (i=0;i<n;++i) {
		//			cout << int(x[i]);
		//		}
		//		cout << ") = ";
		//		cout << WHITE << en << NORMAL << endl;
		//	}

		//}
		return best_energy;
	}

	template<typename real>
	real BipartiteVertexPacking<real>::energy(const std::vector<signed char>& xm, const std::vector<signed char>& xn) const
	{
		size_t m = wm.size();
		size_t n = wn.size();
		if (!is_feasible(xm,xn)) {
			return -1;
		}
		real energy = 0;
		for (size_t i=0;i<m;++i) {
			energy += xm[i]*wm[i];
		}
		for (size_t i=0;i<n;++i) {
			energy += xn[i]*wn[i];
		}
		return energy;
	}

	template<typename real>
	bool BipartiteVertexPacking<real>::is_feasible(const std::vector<signed char>& xm, const std::vector<signed char>& xn) const
	{
		bool feasible = true;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			if (xm.at(itr->first)!=0 && xn.at(itr->second)!=0) {
				feasible = false;
				break;
			}
		}
		return feasible;
	}


	template<typename real>
	void BipartiteVertexPacking<real>::print()
	{
		using namespace std;
		size_t m = wm.size();
		size_t n = wn.size();
		size_t i=0,j=0;
		for (;i<m && j<n;++i,++j) {
			cout << setw(5) << wm[i] << setw(5) << wn[j] << endl;
		}
		for (;i<m;++i) {
			cout << setw(5) << wm[i] << endl;
		}
		for (;j<n;++j) {
			cout << setw(5) << ' ' << setw(5) << wn[j] << endl;
		}
		cout << endl;
		for (auto itr = E.begin(); itr != E.end(); ++itr) {
			// x[i]=1 => x[j]=0
			cout << itr->first << " <--> " << itr->second << endl;
		}
	}
}

#include "vertex_instances.inc"

#ifdef USE_IBFS
#include "ibfs.cpp"
#endif

