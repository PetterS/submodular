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


namespace Petter
{

	// This is the error function used by the QPBO class
	void err_function_qpbo(char* msg)
	{
		throw std::runtime_error(msg);
	}

	template<typename real>
	real PseudoBoolean<real>::minimize_reduction(vector<label>& x) const
	{
		int nlabelled;
		return minimize_reduction(x,nlabelled);
	}

	template<typename real>
	real PseudoBoolean<real>::minimize_reduction(vector<label>& x, int& nlabelled) const
	{
		index nVars = index( x.size() );
		// Here it is important that all indices appear as pairs 
		// in aij or ai
		ASSERT_STR( nvars() <= nVars , "x too small");

		PBF<real, 4> hocr;

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
			real E[8] = {0,0,0,0, 0,0,0,0};
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
			real   E[16] = {0,0,0,0, 0,0,0,0,
			                0,0,0,0, 0,0,0,0};
			E[15] = a;
			hocr.AddHigherTerm(4, ind, E);
		}


		PBF<real,2> qpbf;
		hocr.toQuadratic(qpbf); 
		hocr.clear(); 
		QPBO<real> qpbo(nVars, qpbf.size(), err_function_qpbo); 
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

		real energy = constant + qpbo.ComputeTwiceLowerBound()/2;
		//double energy = eval(x); //Only when x is fully labelled
		return energy;
	}


	
	template<typename real>
	real PseudoBoolean<real>::minimize_reduction_fixetal(vector<label>& x) const
	{
		int nlabelled;
		return minimize_reduction_fixetal(x,nlabelled);
	}

	template<typename real>
	real PseudoBoolean<real>::minimize_reduction_fixetal(vector<label>& x, int& nlabelled) const
	{
		int n = index( x.size() );

		// Work with a copy of the coefficients
		map<triple, real> aijk  = this->aijk;
		map<pair, real> aij  = this->aij;
		map<int, real> ai  = this->ai;
		real constant = this->constant;

		// The quadratic function we are reducing to
		int nedges = aij.size() + aijk.size() + aijkl.size() + 1000;
		QPBO<real> qpbo(n + aijkl.size() + aijk.size() + 1000, nedges, err_function_qpbo);
		qpbo.AddNode(n);


		//
		// Step 1: reduce all positive higher-degree terms
		//

		// Go through the variables one by one and perform the 
		// reductions of the quartic terms

		//for (int ind=0; ind<n; ++ind) {

		//	// Collect all terms with positive coefficients
		//	// containing variable i
		//	real alpha_sum = 0;

		//	// Holds new var
		//	int y = -1;

		//	for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
		//		int i = get_i(itr->first);
		//		int j = get_j(itr->first);
		//		int k = get_k(itr->first);
		//		int l = get_l(itr->first);
		//		real a = itr->second;

		//		// We only have to test for ind==i because of the 
		//		// order we process the indices
		//		if (ind==i && a > 0) {
		//			alpha_sum += a;

		//			// Add term of degree 3
		//			aijk[ make_triple(j,k,l) ]  += a;

		//			// Add negative term of degree 4
		//			// -a*y*xj*xk*xl
		//			if (y<0) y = qpbo.AddNode();
		//			int z = qpbo.AddNode();
		//			qpbo.AddUnaryTerm(z, 0, 3*a);
		//			qpbo.AddPairwiseTerm(z,y, 0,0,0, -a);
		//			qpbo.AddPairwiseTerm(z,j, 0,0,0, -a);
		//			qpbo.AddPairwiseTerm(z,k, 0,0,0, -a);
		//			qpbo.AddPairwiseTerm(z,l, 0,0,0, -a);
		//		}
		//	}

		//	for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
		//		int i = get_i(itr->first);
		//		int j = get_j(itr->first);
		//		int k = get_k(itr->first);
		//		real a = itr->second;

		//		// We only have to test for ind==i because of the 
		//		// order we process the indices
		//		if (ind==i && a > 0) {
		//			alpha_sum += a;

		//			// Add term of degree 2
		//			qpbo.AddPairwiseTerm(j,k, 0,0,0, a);

		//			// Add negative term of degree 3
		//			// -a*y*xj*xk
		//			if (y<0) y = qpbo.AddNode();
		//			int z = qpbo.AddNode();
		//			qpbo.AddUnaryTerm(z, 0, 2*a);
		//			qpbo.AddPairwiseTerm(z,y, 0,0,0, -a);
		//			qpbo.AddPairwiseTerm(z,j, 0,0,0, -a);
		//			qpbo.AddPairwiseTerm(z,k, 0,0,0, -a);
		//		}
		//	}

		//	if (alpha_sum > 0) {
		//		// Add the new quadratic term
		//		qpbo.AddPairwiseTerm(y,ind, 0,0,0, alpha_sum);
		//	}
		//}

		//
		// This code should be equivalent to the commented
		// block above, but faster
		//

		vector<real> alpha_sum(n, 0);
		vector<int>  y(n, -1);

		for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real a = itr->second;

			if (a > 0) {
				alpha_sum[i] += a;

				// Add term of degree 3
				aijk[ make_triple(j,k,l) ]  += a;

				// Add negative term of degree 4
				// -a*y*xj*xk*xl
				if (y[i]<0) y[i] = qpbo.AddNode();
				int z = qpbo.AddNode();
				qpbo.AddUnaryTerm(z, 0, 3*a);
				qpbo.AddPairwiseTerm(z,y[i], 0,0,0, -a);
				qpbo.AddPairwiseTerm(z,j, 0,0,0, -a);
				qpbo.AddPairwiseTerm(z,k, 0,0,0, -a);
				qpbo.AddPairwiseTerm(z,l, 0,0,0, -a);
			}
		}

		for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real a = itr->second;

			// We only have to test for ind==i because of the 
			// order we process the indices
			if (a > 0) {
				alpha_sum[i] += a;

				// Add term of degree 2
				//qpbo.AddPairwiseTerm(j,k, 0,0,0, a);
				aij[ make_pair(j,k) ] += a;

				// Add negative term of degree 3
				// -a*y*xj*xk
				if (y[i]<0) y[i] = qpbo.AddNode();
				/*int z = qpbo.AddNode();
				qpbo.AddUnaryTerm(z, 0, 2*a);
				qpbo.AddPairwiseTerm(z,y[i], 0,0,0, -a);
				qpbo.AddPairwiseTerm(z,j, 0,0,0, -a);
				qpbo.AddPairwiseTerm(z,k, 0,0,0, -a);*/
				aijk[ make_triple(y[i],j,k) ] += -a;
			}
		}

		for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			real& a = itr->second;

			// We only have to test for ind==i because of the 
			// order we process the indices
			if (a > 0) {
				alpha_sum[i] += a;

				// Add term of degree 1
				ai[ j ] += a;

				// Add negative term of degree 2
				// -a*y*xj
				if (y[i]<0) y[i] = qpbo.AddNode();
				aij[ make_pair(y[i],j) ] += -a;

				// Now remove this term
				a = 0;
			}
		}

		for (auto itr=ai.begin(); itr != ai.end(); ++itr) {
			int i =itr->first;
			real& a = itr->second;

			// We only have to test for ind==i because of the 
			// order we process the indices
			if (a > 0) {
				alpha_sum[i] += a;

				// Add term of degree 0
				constant += a;

				// Add negative term of degree 1
				// -a*y*xj
				if (y[i]<0) y[i] = qpbo.AddNode();
				ai[ y[i] ] += -a;

				// Now remove this term
				a = 0;
			}
		}


		for (int i=0;i<n;++i) {
			if (alpha_sum[i] > 0) {
				// Add the new quadratic term
				qpbo.AddPairwiseTerm(y[i],i, 0,0,0, alpha_sum[i]);
			}
		}


		//
		// Done with reducing all positive higher-degree terms
		//

		// Add all negative quartic terms
		for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
			real a = itr->second;
			if (a < 0) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				int k = get_k(itr->first);
				int l = get_l(itr->first);

				int z = qpbo.AddNode();
				qpbo.AddUnaryTerm(z, 0, -3*a);
				qpbo.AddPairwiseTerm(z,i, 0,0,0, a);
				qpbo.AddPairwiseTerm(z,j, 0,0,0, a);
				qpbo.AddPairwiseTerm(z,k, 0,0,0, a);
				qpbo.AddPairwiseTerm(z,l, 0,0,0, a);
			}
		}

		// Add all negative cubic terms
		for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
			real a = itr->second;
			if (a < 0) {
				int i = get_i(itr->first);
				int j = get_j(itr->first);
				int k = get_k(itr->first);

				int z = qpbo.AddNode();
				qpbo.AddUnaryTerm(z, 0, -2*a);
				qpbo.AddPairwiseTerm(z,i, 0,0,0, a);
				qpbo.AddPairwiseTerm(z,j, 0,0,0, a);
				qpbo.AddPairwiseTerm(z,k, 0,0,0, a);
			}
		}


		// Add all quadratic terms
		for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
			real a = itr->second;
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			qpbo.AddPairwiseTerm(i,j, 0,0,0, a);
		}

		// Add all linear terms
		for (auto itr=ai.begin(); itr != ai.end(); ++itr) {
			real a = itr->second;
			int i = itr->first;
			qpbo.AddUnaryTerm(i, 0, a);
		}

		qpbo.MergeParallelEdges();
		qpbo.Solve();
		qpbo.ComputeWeakPersistencies();

		nlabelled = 0;
		for (int i=0; i<n; ++i) {
			x[i] = qpbo.GetLabel(i);
			if (x[i] >= 0) {
				nlabelled++;
			}
		}

		real energy = constant + qpbo.ComputeTwiceLowerBound()/2;

		return energy;
	}

}



#include "pb_instances.inc"
