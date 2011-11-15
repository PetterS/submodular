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
	template<typename real>
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
				ASSERT_STR(c <= 1e-9, "Function not submodular");
				c = std::min(real(0), c);

				constant += c;
				graph.add_tweights(vars[1], 0,  -c);

				graph.add_edge(vars[0],vars[1], -c, 0);
			}
		}
		
		return constant;
	}
#endif

	// This is the error function used by the graph class
	void err_function(char* msg)
	{
		throw std::runtime_error(msg);
	}

#ifdef USE_HOCR
	//typedef PBF<opttype,3> OPTIMIZER;
	#define OPTIMIZER(real) PBF<real,3>
#else
	//typedef Minimizer<real> OPTIMIZER;
	#define OPTIMIZER(real) Minimizer<real>
#endif
	template<typename real> void reduce_bijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_cijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_dijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_eijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_pijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_qijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_rijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);
	template<typename real> void reduce_sijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);

	template<typename real>
	real SymmetricPseudoBoolean<real>::minimize(vector<label>& x) const
	{
		int l;
		return minimize(x,l);
	}

	template<typename real>
	real SymmetricPseudoBoolean<real>::minimize(vector<label>& x, int& nlabelled) const
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
			real val = itr->second;
			hocr.AddUnaryTerm(xi, 0, val);
			hocr.AddUnaryTerm(yi, val, 0);
		}

		// bij
		for (auto itr=bij.begin(); itr != bij.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			real val = itr->second;

			hocr.AddPairwiseTerm(xi,xj, 0,  0,0, val);
			hocr.AddPairwiseTerm(yi,yj, val,0,0, 0);
		}

		// cij
		for (auto itr=cij.begin(); itr != cij.end(); ++itr) {
			int xi = get_i(itr->first);
			int xj = get_j(itr->first);
			int yi = xi + nVars;
			int yj = xj + nVars;
			real val = itr->second;

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
			real w = itr->second;
			real vals1[] = {0,0,0,0,0,0,0,w}; //111 
			real vals2[] = {w,0,0,0,0,0,0,0}; //000
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
			real w = itr->second;
			real   val1[] = {0, //000
			                 0, //001
			                 0, //010
			                 0, //011
			                 0, //100
			                 0, //101
			                 w, //110
			                 0  //111
			                };
			real   val2[] = {0, //000
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
			real w = itr->second;
			real   val1[] = {0, //000
			                 0, //001
			                 0, //010
			                 0, //011
			                 0, //100
			                 w, //101
			                 0, //110
			                 0  //111
			                };
			real   val2[] = {0, //000
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
			real w = itr->second;
			real   val1[] = {0, //000
			                 0, //001
			                 0, //010
			                 w, //011
			                 0, //100
			                 0, //101
			                 0, //110
			                 0  //111
			                };
			real   val2[] = {0, //000
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
			real w = itr->second;

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
			real w = itr->second;

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
			real w = itr->second;

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
			real w = itr->second;

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
			real w = itr->second;

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
			real w = itr->second;

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
			real w = itr->second;

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
			real w = itr->second;

			reduce_sijkl(xi,xj,xk,xl,yi,yj,yk,yl, w, var, hocr);
		}

		
		real energy2;
		vector<label> y(nVars,0);

#ifdef USE_HOCR

		real maxflowconstant;

		//Make sure all vars are used
		hocr.AddUnaryTerm(2*nVars-1,1e-100,0);

		// Hack to make the solver stay away from (0,0) solutions if 
		// other solutions exist 
		for (int i=0;i<nVars; ++i) {
			hocr.AddUnaryTerm(i, 0, 1e-6);
			hocr.AddUnaryTerm(i+nVars, 1e-6, 0);
		}

		// Convert the submodular, reducible polynomial
		// to a quadratic one
		PBF<double,2> qpbf;
		hocr.toQuadratic(qpbf); 
		hocr.clear();
		
		// Convert the quadratic problem into a maximum flow problem
		Graph<real,real,real> graph(qpbf.maxID() + 1, qpbf.size(), err_function); 
		maxflowconstant = convert(graph, qpbf);

		//Solve maximum flow problem
		// Divide by two because we didn't include all all the 1/2's 
		// in front of the summations
		energy2 = constant + (maxflowconstant + graph.maxflow()) / real(2);

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

		// Divide by two because we didn't include all all the 1/2's 
		// in front of the summations
		energy2 = constant + hocr.minimize() / real(2);


		vector<std::pair<int,int> > pairs(nVars);
		for (int i=0;i<nVars;++i) {
			pairs.at(i).first = i;
			pairs.at(i).second = i+nVars;
		}

		hocr.resolve_different(pairs);

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
				if (x[i]<0) {
					x[i]=0;
				}
				nlabelled++;
			}
		}
#endif


		//Calculate energy from solution
		real energy = eval(x,y);

		// Make sure the two computed energies agree
		// If this assertion fails, it might be due to floating-point rounding
		// errors, or that a bug has been found.
		if ( abs(energy) > 1e-5 && abs(energy-energy2)/abs(energy) > 1e-5 ) {
			std::cout << "\n\n";
			std::cout << "eval(x,y) : " << energy << "\n";
			std::cout << "maxflow   : " << energy2 << "\n";
			throw std::runtime_error("Energy mismatch");
		}

		// Mark unlabeled nodes as -1
		for (int i=0; i<nVars; ++i) {
			if (x[i] == y[i] && var_used[i]) {
				x[i] = -1;
			}
		}

		return energy;
	}




	template<typename real> 
	void reduce_bijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

	template<typename real> 
	void reduce_cijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = a; //E111
			int ind[3] = {x1,x2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,y3};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	template<typename real> 
	void reduce_dijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = a; //E111
			int ind[3] = {x1,x2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	template<typename real> 
	void reduce_eijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = a; //E111
			int ind[3] = {x1,x3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y3,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	template<typename real> 
	void reduce_pijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x2,x3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(x1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y2,y3,y4};
			hocr.AddHigherTerm(3, ind, E);}
		}
	}

	template<typename real> 
	void reduce_qijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x2, 0,0,0, -abs(a));


			hocr.AddUnaryTerm(w, 3*abs(a), 0);
			hocr.AddPairwiseTerm(y1,w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y2,w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3,w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4,w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y2, -abs(a), 0,0,0);
		}
	}

	template<typename real> 
	void reduce_rijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,y2,x3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,x3,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x3, 0,0,0, -abs(a));



			hocr.AddUnaryTerm(w, abs(a)*3, 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x4, w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,x2,y3};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,y3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y3, -abs(a), 0,0,0);
		}
	}

	template<typename real> 
	void reduce_sijkl(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr)
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = -abs(a); //E111
			int ind[3] = {x1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = -abs(a); //E000
			int ind[3] = {y1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
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

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,y2,x4};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[7] = abs(a); //E111
			int ind[3] = {x1,y3,x4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(x1,x4, 0,0,0, -abs(a));

			hocr.AddUnaryTerm(w, abs(a)*3, 0);
			hocr.AddPairwiseTerm(y1, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x2, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(x3, w, -abs(a), 0,0,0);
			hocr.AddPairwiseTerm(y4, w, -abs(a), 0,0,0);

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,x2,y4};
			hocr.AddHigherTerm(3, ind, E);}

			{real E[8]={0,0,0,0, 0,0,0,0};
			E[0] = abs(a); //E000
			int ind[3] = {y1,x3,y4};
			hocr.AddHigherTerm(3, ind, E);}

			hocr.AddPairwiseTerm(y1,y4, -abs(a), 0,0,0);
		}
	}





//TODO: change this into nicer code
#define INSTANTIATE_REDUCTIONS(real) \
	template void reduce_bijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_cijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_dijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_eijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_pijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_qijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_rijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr); \
	template void reduce_sijkl<real>(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4, real a, int& var, OPTIMIZER(real)& hocr);

INSTANTIATE_REDUCTIONS(int);
INSTANTIATE_REDUCTIONS(double);

}



#include "pb_instances.inc"