//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Heuristics for creating a submodular relaxation g(x,y)
// Note that this entire file contains unpublished work. 
//


#include "PseudoBoolean.h"

#include <tuple>
#include <algorithm>


namespace Petter
{

	// Helper function which solves
	// [0 1 1] [x1]   [b1]
	// [1 0 1] [x2] = [b2]
	// [1 1 0] [x3]   [b3]
	//
	template<typename real>
	void solve_system(real& x1, real& x2, real& x3, real b1, real b2, real b3) 
	{
		// TODO: change representation of g so that 
		// this can be done in integer arithmetic.

		x1 = (-b1 + b2 + b3) / 2;
		x2 = ( b1 - b2 + b3) / 2;
		x3 = ( b1 + b2 - b3) / 2;
	}

	template<typename real> template<typename orgreal>
	void SymmetricPseudoBoolean<real>::create_heuristic(PseudoBoolean<orgreal>& pbf)
	{
		using namespace std;

		// There is no freedom in picking these
		constant = real(pbf.constant);
		for (auto itr=pbf.ai.begin(); itr != pbf.ai.end(); ++itr) {
			bi[itr->first] = static_cast<real>( itr->second );
		}	

		// Right-hand sides in the submodular inequalities
		map<pair,real> bRHS;
		map<pair,real> cRHS;
		for (auto itr=pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			bRHS[itr->first]  = 0;
			cRHS[itr->first] = real(itr->second);
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			triple ind = itr->first;
			int i = get_i(ind);
			int j = get_j(ind);
			int k = get_k(ind);
			real a123 = real(itr->second);
		
			if (a123 > 0) {

				real m01 = max(bRHS[make_pair(i,j)] - cRHS[make_pair(i,j)] , real(0));
				real m02 = max(bRHS[make_pair(i,k)] - cRHS[make_pair(i,k)] , real(0));
				real m12 = max(bRHS[make_pair(j,k)] - cRHS[make_pair(j,k)] , real(0));

				real bcanuse = max(m01+m02+m12 - a123, real(0)) / 2;
				real bmargin = min(m01, min(m02,m12));
				if (bmargin > 0) {
					real m = min(bmargin, min(a123, bcanuse));
					bijk[ind] +=  m;
					a123 -= m;
					m01 -= m;
					m02 -= m;
					m12 -= m;
					bRHS[make_pair(i,j)] -= m;
					bRHS[make_pair(i,k)] -= m;
					bRHS[make_pair(i,k)] -= m;
				}
		

				while (a123 > 0 && (m01>0 || m02>0 || m12>0) ) {  
					if (m01 > 0) {
						real m = min(m01, a123);
						cijk[ind] += m;
						a123 -= m;
						m01  -= m;
						bRHS[make_pair(i,j)] -= m;
					}
					else if (m02 > 0) {
						real m = min(m02, a123);
						dijk[ind] +=  m;
						a123 -= m;
						m02  -= m;
						bRHS[make_pair(i,k)] -= m;
					}
					else if (m12 > 0) {
						real m = min(m12, a123);
						eijk[ind] +=  m;
						a123 -= m;
						m12  -= m;
						bRHS[make_pair(j,k)] -= m;
					}
				}
				
				dijk[ind] += a123;
				bRHS[make_pair(i,k)] -= a123;
			}
			else {
				
				real m01 = max(cRHS[make_pair(i,j)] - bRHS[make_pair(i,j)] , real(0));
				real m02 = max(cRHS[make_pair(i,k)] - bRHS[make_pair(i,k)] , real(0));
				real m12 = max(cRHS[make_pair(j,k)] - bRHS[make_pair(j,k)] , real(0));


				real x0,x1,x2;
				solve_system(x0,x1,x2, -m01, -m02, -m12);


				if (x0<=0 && x1<=0 && x2<=0) {
					real m = max(real(0),min(-a123, -x0));
					cijk[ind] -= m;
					a123 += m;
					cRHS[make_pair(i,k)] -= m;
					cRHS[make_pair(j,k)] -= m;

					m = max(real(0),min(-a123, -x1));
					dijk[ind] -= m;
					a123 += m;
					cRHS[make_pair(i,j)] -= m;
					cRHS[make_pair(j,k)] -= m;

					m = max(real(0),min(-a123, -x2));
					eijk[ind] -= m;
					a123 += m;
					cRHS[make_pair(i,j)] -= m;
					cRHS[make_pair(i,k)] -= m;
				}
				else {
					real mct,mdt,met;
					do {
						mct = min(m02,m12);
						mdt = min(m01,m12);
						met = min(m01,m02);
						if (mct>=mdt && mct>=met) {
							real m = min(-a123,mct);
							cijk[ind] -= m;
							a123 += m;
							m02 -= m;
							m12 -= m;
							cRHS[make_pair(i,k)] -= m;
							cRHS[make_pair(j,k)] -= m;
						}
						else if (mdt>=mct && mdt>=met) {
							real m = min(-a123,mdt);
							dijk[ind] -= m;
							a123 += m;
							m01 -= m;
							m12 -= m;
							cRHS[make_pair(i,j)] -= m;
							cRHS[make_pair(j,k)] -= m;
						}
						else if (met>=mct && met>=mdt) {
							real m = min(-a123,met);
							eijk[ind] -= m;
							a123 += m;
							m01 -= m;
							m02 -= m;
							cRHS[make_pair(i,j)] -= m;
							cRHS[make_pair(i,k)] -= m;
						}

					} while (max(mct,max(mdt,met)) > 0 && a123<0);
				}
				

				bijk[ind] += a123;

			} // a123>=0
		}// for aijk


		// Quartic case
		// TODO: make this a little more sophisticated :-)
		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			quad ijkl = itr->first;
			real a0123 = real(itr->second);

			if (a0123 > 0) {
				dijkl[ijkl] = a0123;
			}
			else {
				bijkl[ijkl] = a0123;
			}
		}


		bRHS.clear();
		cRHS.clear();
		for (auto itr=pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			bRHS[itr->first]  = 0;
			cRHS[itr->first] = real(itr->second);
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			triple ind = itr->first;
			int i = get_i(ind);
			int j = get_j(ind);
			int k = get_k(ind);

			bRHS[make_pair(i,j)] -= max(bijk[ind],real(0));
			bRHS[make_pair(i,k)] -= max(bijk[ind],real(0));
			bRHS[make_pair(j,k)] -= max(bijk[ind],real(0));
			bRHS[make_pair(i,j)] -= max(cijk[ind],real(0));
			bRHS[make_pair(i,k)] -= max(dijk[ind],real(0));
			bRHS[make_pair(j,k)] -= max(eijk[ind],real(0));

			cRHS[make_pair(i,j)] += min(dijk[ind],real(0));
			cRHS[make_pair(i,j)] += min(eijk[ind],real(0));
			cRHS[make_pair(i,k)] += min(cijk[ind],real(0));
			cRHS[make_pair(i,k)] += min(eijk[ind],real(0));
			cRHS[make_pair(j,k)] += min(cijk[ind],real(0));
			cRHS[make_pair(j,k)] += min(dijk[ind],real(0));
		}

		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			quad ijkl = itr->first;
			int i = get_i(ijkl);
			int j = get_j(ijkl);
			int k = get_k(ijkl);
			int l = get_l(ijkl);

			bRHS[make_pair(i,j)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(i,k)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(i,l)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(j,k)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(j,l)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(k,l)] -= max(bijkl[ijkl],real(0));

			bRHS[make_pair(i,j)] -= abs(cijkl[ijkl]);
			bRHS[make_pair(i,k)] -= abs(cijkl[ijkl]);
			bRHS[make_pair(j,k)] -= abs(cijkl[ijkl]);

			bRHS[make_pair(i,j)] -= abs(dijkl[ijkl]);
			bRHS[make_pair(i,l)] -= abs(dijkl[ijkl]);
			bRHS[make_pair(j,l)] -= abs(dijkl[ijkl]);

			bRHS[make_pair(i,k)] -= abs(eijkl[ijkl]);
			bRHS[make_pair(i,l)] -= abs(eijkl[ijkl]);
			bRHS[make_pair(k,l)] -= abs(eijkl[ijkl]);

			bRHS[make_pair(j,k)] -= abs(pijkl[ijkl]);
			bRHS[make_pair(j,l)] -= abs(pijkl[ijkl]);
			bRHS[make_pair(k,l)] -= abs(pijkl[ijkl]);
		
			bRHS[make_pair(i,j)] -= abs(qijkl[ijkl]) + max(qijkl[ijkl], real(0));
			bRHS[make_pair(k,l)] -= max(qijkl[ijkl], real(0));

			bRHS[make_pair(i,k)] -= abs(rijkl[ijkl]) + max(rijkl[ijkl], real(0));
			bRHS[make_pair(j,l)] -= max(rijkl[ijkl], real(0));

			bRHS[make_pair(i,l)] -= abs(sijkl[ijkl]) + max(sijkl[ijkl], real(0));
			bRHS[make_pair(j,k)] -= max(sijkl[ijkl], real(0));



			cRHS[make_pair(i,l)] += min(cijkl[ijkl], real(0));
			cRHS[make_pair(j,l)] += min(cijkl[ijkl], real(0));
			cRHS[make_pair(k,l)] += min(cijkl[ijkl], real(0));

			cRHS[make_pair(i,k)] += min(dijkl[ijkl], real(0));
			cRHS[make_pair(j,k)] += min(dijkl[ijkl], real(0));
			cRHS[make_pair(k,l)] += min(dijkl[ijkl], real(0));

			cRHS[make_pair(i,j)] += min(eijkl[ijkl], real(0));
			cRHS[make_pair(j,k)] += min(eijkl[ijkl], real(0));
			cRHS[make_pair(j,l)] += min(eijkl[ijkl], real(0));

			cRHS[make_pair(i,j)] += min(pijkl[ijkl], real(0));
			cRHS[make_pair(i,k)] += min(pijkl[ijkl], real(0));
			cRHS[make_pair(i,l)] += min(pijkl[ijkl], real(0));

			cRHS[make_pair(i,k)] -= abs(qijkl[ijkl]);
			cRHS[make_pair(i,l)] -= abs(qijkl[ijkl]);
			cRHS[make_pair(j,k)] -= abs(qijkl[ijkl]);
			cRHS[make_pair(j,l)] -= abs(qijkl[ijkl]);

			cRHS[make_pair(i,j)] -= abs(rijkl[ijkl]);
			cRHS[make_pair(i,l)] -= abs(rijkl[ijkl]);
			cRHS[make_pair(j,k)] -= abs(rijkl[ijkl]);
			cRHS[make_pair(k,l)] -= abs(rijkl[ijkl]);

			cRHS[make_pair(i,j)] -= abs(sijkl[ijkl]);
			cRHS[make_pair(i,k)] -= abs(sijkl[ijkl]);
			cRHS[make_pair(j,l)] -= abs(sijkl[ijkl]);
			cRHS[make_pair(k,l)] -= abs(sijkl[ijkl]);
		}

		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			pair ind = itr->first;
			real aij = real(itr->second);

			bij[ind] = min( bRHS[ind], cRHS[ind] );
			cij[ind] = aij - bij[ind];
		}

	} 


	template<typename real>
	void SymmetricPseudoBoolean<real>::make_submodular(const PseudoBoolean<real>& pbf)
	{
		using namespace std;

		map<pair,real> bRHS;
		map<pair,real> cRHS;

		for (auto itr=pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			bRHS[itr->first]  = 0;
			cRHS[itr->first] = real(itr->second);
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			triple ind = itr->first;
			int i = get_i(ind);
			int j = get_j(ind);
			int k = get_k(ind);

			bRHS[make_pair(i,j)] -= max(bijk[ind],real(0));
			bRHS[make_pair(i,k)] -= max(bijk[ind],real(0));
			bRHS[make_pair(j,k)] -= max(bijk[ind],real(0));
			bRHS[make_pair(i,j)] -= max(cijk[ind],real(0));
			bRHS[make_pair(i,k)] -= max(dijk[ind],real(0));
			bRHS[make_pair(j,k)] -= max(eijk[ind],real(0));

			cRHS[make_pair(i,j)] += min(dijk[ind],real(0));
			cRHS[make_pair(i,j)] += min(eijk[ind],real(0));
			cRHS[make_pair(i,k)] += min(cijk[ind],real(0));
			cRHS[make_pair(i,k)] += min(eijk[ind],real(0));
			cRHS[make_pair(j,k)] += min(cijk[ind],real(0));
			cRHS[make_pair(j,k)] += min(dijk[ind],real(0));
		}

		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			quad ijkl = itr->first;
			int i = get_i(ijkl);
			int j = get_j(ijkl);
			int k = get_k(ijkl);
			int l = get_l(ijkl);

			bRHS[make_pair(i,j)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(i,k)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(i,l)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(j,k)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(j,l)] -= max(bijkl[ijkl],real(0));
			bRHS[make_pair(k,l)] -= max(bijkl[ijkl],real(0));

			bRHS[make_pair(i,j)] -= abs(cijkl[ijkl]);
			bRHS[make_pair(i,k)] -= abs(cijkl[ijkl]);
			bRHS[make_pair(j,k)] -= abs(cijkl[ijkl]);

			bRHS[make_pair(i,j)] -= abs(dijkl[ijkl]);
			bRHS[make_pair(i,l)] -= abs(dijkl[ijkl]);
			bRHS[make_pair(j,l)] -= abs(dijkl[ijkl]);

			bRHS[make_pair(i,k)] -= abs(eijkl[ijkl]);
			bRHS[make_pair(i,l)] -= abs(eijkl[ijkl]);
			bRHS[make_pair(k,l)] -= abs(eijkl[ijkl]);

			bRHS[make_pair(j,k)] -= abs(pijkl[ijkl]);
			bRHS[make_pair(j,l)] -= abs(pijkl[ijkl]);
			bRHS[make_pair(k,l)] -= abs(pijkl[ijkl]);
		
			bRHS[make_pair(i,j)] -= abs(qijkl[ijkl]) + max(qijkl[ijkl], real(0));
			bRHS[make_pair(k,l)] -= max(qijkl[ijkl], real(0));

			bRHS[make_pair(i,k)] -= abs(rijkl[ijkl]) + max(rijkl[ijkl], real(0));
			bRHS[make_pair(j,l)] -= max(rijkl[ijkl], real(0));

			bRHS[make_pair(i,l)] -= abs(sijkl[ijkl]) + max(sijkl[ijkl], real(0));
			bRHS[make_pair(j,k)] -= max(sijkl[ijkl], real(0));



			cRHS[make_pair(i,l)] += min(cijkl[ijkl], real(0));
			cRHS[make_pair(j,l)] += min(cijkl[ijkl], real(0));
			cRHS[make_pair(k,l)] += min(cijkl[ijkl], real(0));

			cRHS[make_pair(i,k)] += min(dijkl[ijkl], real(0));
			cRHS[make_pair(j,k)] += min(dijkl[ijkl], real(0));
			cRHS[make_pair(k,l)] += min(dijkl[ijkl], real(0));

			cRHS[make_pair(i,j)] += min(eijkl[ijkl], real(0));
			cRHS[make_pair(j,k)] += min(eijkl[ijkl], real(0));
			cRHS[make_pair(j,l)] += min(eijkl[ijkl], real(0));

			cRHS[make_pair(i,j)] += min(pijkl[ijkl], real(0));
			cRHS[make_pair(i,k)] += min(pijkl[ijkl], real(0));
			cRHS[make_pair(i,l)] += min(pijkl[ijkl], real(0));

			cRHS[make_pair(i,k)] -= abs(qijkl[ijkl]);
			cRHS[make_pair(i,l)] -= abs(qijkl[ijkl]);
			cRHS[make_pair(j,k)] -= abs(qijkl[ijkl]);
			cRHS[make_pair(j,l)] -= abs(qijkl[ijkl]);

			cRHS[make_pair(i,j)] -= abs(rijkl[ijkl]);
			cRHS[make_pair(i,l)] -= abs(rijkl[ijkl]);
			cRHS[make_pair(j,k)] -= abs(rijkl[ijkl]);
			cRHS[make_pair(k,l)] -= abs(rijkl[ijkl]);

			cRHS[make_pair(i,j)] -= abs(sijkl[ijkl]);
			cRHS[make_pair(i,k)] -= abs(sijkl[ijkl]);
			cRHS[make_pair(j,l)] -= abs(sijkl[ijkl]);
			cRHS[make_pair(k,l)] -= abs(sijkl[ijkl]);
		}

		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			pair ind = itr->first;
			real aij = real(itr->second);

			bij[ind] = min( bRHS[ind], cRHS[ind] );
			cij[ind] = aij - bij[ind];
		}
	}


	// TODO: better solution to this?
	template void SymmetricPseudoBoolean<double>::create_heuristic(PseudoBoolean<double>& pbf);
	template void SymmetricPseudoBoolean<double>::create_heuristic(PseudoBoolean<int>& pbf);
	template void SymmetricPseudoBoolean<int>::create_heuristic(PseudoBoolean<double>& pbf);
	template void SymmetricPseudoBoolean<int>::create_heuristic(PseudoBoolean<int>& pbf);


} //namespace Petter


#include "pb_instances.inc"
