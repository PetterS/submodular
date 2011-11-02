//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Functions to minimize a general f(x) and a submodular g(x,y) using reductions
//

#include "PseudoBoolean.h"

#include "QPBO.h"
#include "HOCR.h"
#include "graph.h"

//#define USE_HOCR //Whether to use HOCR for representing quadratic functions
#ifndef USE_HOCR
	// Degree-3 minimizer
	#include "Minimizer.h"
#endif

namespace Petter
{

#ifdef USE_HOCR
	//
	// Converts a submodular quadratic pseudo-Boolean function into 
	// a graph for minimization via maxflow
	//
	real convert(Graph<real,real,real>& graph, PBF<real,2>& qpbf)
	{
		//The constant always seem to be zero
		real constant = qpbf.cnst();

		ASSERT_STR(graph.get_node_num() == 0, "Graph should be empty");
		graph.add_node(qpbf.maxID() + 1);

		PBF<real,2>::VID vars[2];
		real c;
		int size;
		qpbf.startEnum();
		while (qpbf.get(size, vars, c))
		{
			if (size == 1) {
				// c*x(i)
				if (c >= 0) {
					graph.add_tweights(vars[0], c, 0);
				}
				else {
					constant += c;
					graph.add_tweights(vars[0], 0,  -c);
				}
			}
			else {
				// c*x(i)*x(j)

				// c has no be non-positive for this function to be 
				// submodular
				ASSERT_STR(c <= 0, "Function not submodular");

				constant += c;
				graph.add_tweights(vars[1], 0,  -c);

				graph.add_edge(vars[0],vars[1], -c, 0);
			}
		}
		
		return constant;
	}
#endif

	// This is the error function used by the QPBO and Graph classes
	void err_function(char* msg)
	{
		throw std::runtime_error(msg);
	}

	real PseudoBoolean::minimize_reduction(vector<label>& x) const
	{
		int nlabelled;
		return minimize_reduction(x,nlabelled);
	}
	real PseudoBoolean::minimize_reduction(vector<label>& x, int& nlabelled) const
	{
		index nVars = index( x.size() );
		// Here it is important that all indices appear as pairs 
		// in aij or ai
		ASSERT_STR( nvars() <= nVars , "x too small");

		PBF<double, 4> hocr;

		for (auto itr=ai.begin(); itr != ai.end(); ++itr) {
			int i = itr->first;
			real a = itr->second;
			hocr.AddUnaryTerm(i, 0, a);
		}

		for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			real a = itr->second;
			hocr.AddPairwiseTerm(i,j, 0,0,0, a);
		}

		for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real a = itr->second;
			int ind[] = {i,j,k};
			double E[8] = {0,0,0,0, 0,0,0,0};
			E[7] = a;
			hocr.AddHigherTerm(3, ind, E);
		}

		for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real a = itr->second;
			int ind[] = {i,j,k,l};
			double E[16] = {0,0,0,0, 0,0,0,0,
			                0,0,0,0, 0,0,0,0};
			E[15] = a;
			hocr.AddHigherTerm(4, ind, E);
		}

		//Make sure all vars are used TODO
		hocr.AddUnaryTerm(nVars-1,1e-100,0);



		PBF<double,2> qpbf;
		hocr.toQuadratic(qpbf); 
		hocr.clear(); 
		QPBO<double> qpbo(nVars, qpbf.size(), err_function); 
		convert(qpbo, qpbf);
		qpbo.MergeParallelEdges();
		qpbo.Solve();
		qpbo.ComputeWeakPersistencies();

		nlabelled = 0;
		for (int i=0; i<nVars; ++i) {
			x[i] = qpbo.GetLabel(i);
			if (x[i] >= 0) {
				nlabelled++;
			}
		}

		double energy = constant + qpbo.ComputeTwiceLowerBound()/2;
		//double energy = eval(x); //Only when x is fully labelled
		return energy;
	}



	typedef double opttype;
#ifdef USE_HOCR
	typedef PBF<opttype,3> OPTIMIZER;
#else
	typedef Minimizer<real> OPTIMIZER;
#endif
	void reduce_bijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_cijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_dijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_eijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_pijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_qijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_rijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);
	void reduce_sijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr);

	real SymmetricPseudoBoolean::minimize(vector<label>& x) const
	{
		int l;
		return minimize(x,l);
	}
	real SymmetricPseudoBoolean::minimize(vector<label>& x, int& nlabelled) const
	{
		index nVars = index( x.size() );

		// Here the assumption that all indices appear as pairs 
		// in bij or bi is important
		map<int, bool> var_used;
		for (auto itr=bij.begin(); itr != bij.end(); ++itr) {
			var_used[get_i(itr->first)] = true;
			var_used[get_j(itr->first)] = true;
			ASSERT_STR(get_i(itr->first) < nVars, "x too small");
			ASSERT_STR(get_j(itr->first) < nVars, "x too small");
		}
		for (auto itr=bi.begin(); itr != bi.end(); ++itr) {
			var_used[itr->first] = true;
			ASSERT_STR(itr->first < nVars, "x too small");
		}


#ifdef USE_HOCR
		PBF<double, 3> hocr;
#else
		Minimizer<real> hocr(2*nVars);
#endif


		int var = 2*nVars; // Current variable

		// bi
		for (auto itr=bi.begin(); itr != bi.end(); ++itr) {
			int xi = itr->first;
			int yi = itr->first + nVars;
			real val = 0.5 * itr->second;
			hocr.AddUnaryTerm(xi, 0, val);
			hocr.AddUnaryTerm(yi, val, 0);
		}

		// bij
		for (auto itr=bij.begin(); itr != bij.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			real val = 0.5 * itr->second;

			hocr.AddPairwiseTerm(xi,xj, 0,  0,0, val);
			hocr.AddPairwiseTerm(yi,yj, val,0,0, 0);
		}

		// cij
		for (auto itr=cij.begin(); itr != cij.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			real val = 0.5 * itr->second;

			hocr.AddPairwiseTerm(xi,yj, 0, 0,   val, 0);
			hocr.AddPairwiseTerm(yi,xj, 0, val, 0,   0);
		}

		// bijk
		for (auto itr=bijk.begin(); itr != bijk.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			real w = 0.5 * itr->second;
			double vals1[] = {0,0,0,0,0,0,0,w}; //111 
			double vals2[] = {w,0,0,0,0,0,0,0}; //000
			int vars1[] = {xi,xj,xk};
			int vars2[] = {yi,yj,yk};

			hocr.AddHigherTerm(3,vars1,vals1);
			hocr.AddHigherTerm(3,vars2,vals2);
		}

		// cijk
		for (auto itr=cijk.begin(); itr != cijk.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			real w = 0.5 * itr->second;
			double val1[] = {0, //000
			                 0, //001
			                 0, //010
			                 0, //011
			                 0, //100
			                 0, //101
			                 w, //110
			                 0  //111
			                };
			double val2[] = {0, //000
			                 w, //001
			                 0, //010
			                 0, //011
			                 0, //100
			                 0, //101
			                 0, //110
			                 0  //111
			                };
			int vars1[] = {xi,xj,yk};
			int vars2[] = {yi,yj,xk};

			hocr.AddHigherTerm(3,vars1,val1);
			hocr.AddHigherTerm(3,vars2,val2);
		}

		// dijk
		for (auto itr=dijk.begin(); itr != dijk.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			real w = 0.5 * itr->second;
			double val1[] = {0, //000
			                 0, //001
			                 0, //010
			                 0, //011
			                 0, //100
			                 w, //101
			                 0, //110
			                 0  //111
			                };
			double val2[] = {0, //000
			                 0, //001
			                 w, //010
			                 0, //011
			                 0, //100
			                 0, //101
			                 0, //110
			                 0  //111
			                };
			int vars1[] = {xi,yj,xk};
			int vars2[] = {yi,xj,yk};

			hocr.AddHigherTerm(3,vars1,val1);
			hocr.AddHigherTerm(3,vars2,val2);
		}

		// eijk
		for (auto itr=eijk.begin(); itr != eijk.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			real w = 0.5 * itr->second;
			double val1[] = {0, //000
			                 0, //001
			                 0, //010
			                 w, //011
			                 0, //100
			                 0, //101
			                 0, //110
			                 0  //111
			                };
			double val2[] = {0, //000
			                 0, //001
			                 0, //010
			                 0, //011
			                 w, //100
			                 0, //101
			                 0, //110
			                 0  //111
			                };
			int vars1[] = {yi,xj,xk};
			int vars2[] = {xi,yj,yk};

			hocr.AddHigherTerm(3,vars1,val1);
			hocr.AddHigherTerm(3,vars2,val2);
		}


		// bijkl
		for (auto itr=bijkl.begin(); itr != bijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_bijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// cijkl
		for (auto itr=cijkl.begin(); itr != cijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_cijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// dijkl
		for (auto itr=dijkl.begin(); itr != dijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_dijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// eijkl
		for (auto itr=eijkl.begin(); itr != eijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_eijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// pijkl
		for (auto itr=pijkl.begin(); itr != pijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_pijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// qijkl
		for (auto itr=qijkl.begin(); itr != qijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_qijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// rijkl
		for (auto itr=rijkl.begin(); itr != rijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_rijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		// sijkl
		for (auto itr=sijkl.begin(); itr != sijkl.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int xk = get_k(itr->first);
			int xl = get_l(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			int yk = xk + nVars;
			int yl = xl + nVars;
			real w = 0.5 * itr->second;

			reduce_sijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		//Make sure all vars are used TODO
		hocr.AddUnaryTerm(2*nVars-1,1e-100,0);

		//Hack to make the solver stay away from (0,0) solutions if others
		//exist TODO
		for (int i=0;i<nVars; ++i) {
			hocr.AddUnaryTerm(i, 0, 1e-6);
			hocr.AddUnaryTerm(i+nVars, 1e-6, 0);
		}


		double maxflowconstant;
		double energy2;
		vector<label> y(nVars,0);

#ifdef USE_HOCR

		// Convert the submodular, reducible polynomial
		// to a quadratic one
		PBF<double,2> qpbf;
		hocr.toQuadratic(qpbf); 
		hocr.clear();
		
		// Convert the quadratic problem into a maximum flow problem
		Graph<real,real,real> graph(qpbf.maxID() + 1, qpbf.size(), err_function); 
		maxflowconstant = convert(graph, qpbf);

		//Solve maximum flow problem
		energy2 = constant + maxflowconstant + graph.maxflow();

		//Extract solution
		nlabelled = 0;
		for (int i=0; i<nVars; ++i) {
			if (var_used[i]) {
				x[i] = graph.what_segment(i);
				y[i] = graph.what_segment(i+nVars);
				if (x[i] != y[i]) {
					nlabelled++;
				}
			}
			else {
				// This variable is not part of the polynomial,
				// therefore labelled
				nlabelled++;
			}
		}


#else
		energy2 = constant + hocr.minimize();

		// Extract solution
		nlabelled = 0;
		for (int i=0; i<nVars; ++i) {
			if (var_used[i]) {
				x[i] = hocr.get_solution(i);
				y[i] = hocr.get_solution(i+nVars);
				if (x[i] != y[i]) {
					nlabelled++;
				}
			}
			else {
				// This variable is not part of the polynomial,
				// therefore labelled
				nlabelled++;
			}
		}
#endif

		//Calculate energy from solution
		double energy = eval(x,y);
		//Make sure the two computed energies agree
		//std::cout << "\n\n";
		//std::cout << energy << "\n";
		//std::cout << energy2 << "\n";
		//If this assertion fails, check that the version of HOCR is >= 1.02
		ASSERT( abs(energy-energy2)/abs(energy) < 1e-5 || abs(energy) < 1e-5);

		// Mark unlabeled nodes as -1
		for (int i=0; i<nVars; ++i) {
			if (x[i] == y[i] && var_used[i]) {
				x[i] = -1;
			}
		}

		return energy;
	}





	void reduce_bijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1,z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x2,z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x3,z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x4,z, 0,0,0, -2*abs(a));

			hocr.AddPairwiseTerm(x1,x2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x3,x4, 0,0,0, abs(a));

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -2*abs(a), 0,0,0);

			hocr.AddPairwiseTerm(y1, y2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1, y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1, y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, y4, abs(a), 0,0,0);
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z,0, abs(a)*3);
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -abs(a));

			hocr.AddUnaryTerm(w, abs(a)*3, 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);
		}
	}

	void reduce_cijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y4, z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = a; //E111
			int ind[3] = {x1,x2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = a; //E000
			int ind[3] = {y1,y2,y3};
			hocr.AddHigherTerm(3, ind, E);}
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y4, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x1,x2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,y4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,y4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x3,y4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x2,x3};
			hocr.AddHigherTerm(3, ind, E);}


			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,x4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3,x4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,y3};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	void reduce_dijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = a; //E111
			int ind[3] = {x1,x2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = a; //E000
			int ind[3] = {y1,y2,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x1,x2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,y3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,y3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y3,x4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x2,x4};
			hocr.AddHigherTerm(3, ind, E);}


			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,x3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3,y4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	void reduce_eijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = a; //E111
			int ind[3] = {x1,x3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = a; //E000
			int ind[3] = {y1,y3,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x1,y2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y2,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y2,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x3,x4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x3,x4};
			hocr.AddHigherTerm(3, ind, E);}


			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3,y4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y3,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	void reduce_pijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(y1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x2,x3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(x1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y2,y3,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(y1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y1,x2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y1,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y1,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x3,x4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x2,x3,x4};
			hocr.AddHigherTerm(3, ind, E);}


			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(x1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x1,y2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x1,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x1,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3,y4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y2,y3,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	void reduce_qijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y4, z, 0,0,0, -2*abs(a));

			hocr.AddPairwiseTerm(x1,x2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,y3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,y4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,y3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x2,y4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y3,y4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x2, 0,0,0, abs(a));



			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -2*abs(a), 0,0,0);

			hocr.AddPairwiseTerm(y1,y2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,x3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,x4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3,x4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y2, abs(a), 0,0,0);
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1,z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x2,z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y3,z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y4,z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x2, 0,0,0, -abs(a));


			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1,w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3,w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4,w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y2, -abs(a), 0,0,0);
		}
	}

	void reduce_rijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y4, z, 0,0,0, -2*abs(a));

			hocr.AddPairwiseTerm(x1,y2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,y4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y2,x3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y2,y4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x3,y4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x3,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x3, 0,0,0, abs(a));



			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -2*abs(a), 0,0,0);

			hocr.AddPairwiseTerm(y1,x2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2,y3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2,x4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3,x4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y3, abs(a), 0,0,0);
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z,0, abs(a)*3);
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y4, z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,x3,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x3, 0,0,0, -abs(a));



			hocr.AddUnaryTerm(w, abs(a)*3, 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,y3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y3, -abs(a), 0,0,0);
		}
	}

	void reduce_sijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER& hocr)
	{
		if (a > 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z, 0, 3*abs(a));
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y2, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(y3, z, 0,0,0, -2*abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -2*abs(a));

			hocr.AddPairwiseTerm(x1,y2, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,y3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(x1,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y2,y3, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y2,x4, 0,0,0, abs(a));
			hocr.AddPairwiseTerm(y3,x4, 0,0,0, abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,y3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x4, 0,0,0, abs(a));



			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -2*abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -2*abs(a), 0,0,0);

			hocr.AddPairwiseTerm(y1,x2, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,x3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y1,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2,x3, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2,y4, abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3,y4, abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,x3,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y4, abs(a), 0,0,0);
		}
		else if (a < 0) {
			//Extra variables needed
			int z = var++;
			int w = var++;

			hocr.AddUnaryTerm(z,0, abs(a)*3);
			hocr.AddPairwiseTerm(x1, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y2, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(y3, z, 0,0,0, -abs(a));
			hocr.AddPairwiseTerm(x4, z, 0,0,0, -abs(a));

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,y3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x4, 0,0,0, -abs(a));

			hocr.AddUnaryTerm(w, abs(a)*3, 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			{opttype E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,x3,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y4, -abs(a), 0,0,0);
		}
	}
}

