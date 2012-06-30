//
// Contains the definition of basic functions, such as creating and evaluating
// the polynomials.
//

#include <iomanip>
#include <fstream>
#include <cmath>

#include "PseudoBoolean.h"

using std::size_t;

namespace Petter
{

	const char* str(Method method)
	{
		#define switch_helper(m) case m: return #m;
		switch (method) {
			switch_helper(HOCR)
			switch_helper(Fix)
			switch_helper(GRD)
			switch_helper(GRD_heur)
			switch_helper(GRD_gen)
			switch_helper(M_reduction)
			switch_helper(LP)
		}
		return "Unknown";
	}

	void check_Boolean(const vector<label>& x)
	{
		for (std::size_t i=0; i<x.size(); ++i) {
			ASSERT( x[i] == 0 || x[i] == 1 );
		}
	}


	template<typename real>
	PseudoBoolean<real>::PseudoBoolean()
	{
		constant = 0;
	}

	template<typename real>
	PseudoBoolean<real>::PseudoBoolean(std::string filename)
	{
		constant = 0;
		std::ifstream fin(filename.c_str());
		ASSERT_STR(fin,"Could not open file");

		int n,i,j,k,l;
		real a;
		fin >> n;
		for (i=0;i<n;++i) {
			fin >> a;
			add_monomial(i, a);
		}
		ASSERT_STR(fin,"Could not read degree 1 monomials");
		fin >> n;
		if (!fin) {
			return;
		}
		for (int c=1;c<=n;++c) {
			fin >> i >> j;
			i--; j--;
			fin >> a;
			add_monomial(i,j,a);
		}
		ASSERT_STR(fin,"Could not read degree 2 monomials");
		fin >> n;
		if (!fin) {
			return;
		}
		for (int c=1;c<=n;++c) {
			fin >> i >> j >> k;
			i--; j--; k--;
			fin >> a;
			add_monomial(i,j,k,a);
		}
		ASSERT_STR(fin,"Could not read degree 3 monomials");
		fin >> n;
		if (!fin) {
			return;
		}
		for (int c=1;c<=n;++c) {
			fin >> i >> j >> k >> l; 
			i--; j--; k--; l--;
			fin >> a;
			add_monomial(i,j,k,l,a);
		}
		ASSERT_STR(fin,"Could not read degree 4 monomials");
	}
	
	template<typename real>
	void PseudoBoolean<real>::save_to_file(std::string filename)
	{
		std::ofstream fout(filename.c_str());
		ASSERT(fout);
		int n = nvars();

		fout << n << ' ';
		for (int i=0;i<n;++i) {
			fout << ai[i] << ' ';
		}
		fout << std::endl;

		fout << aij.size() << ' ';
		for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			real aij = itr->second;
			fout << i+1 << ' ' << j+1 << ' ' << aij << ' ';
		}
		fout << std::endl;

		fout << aijk.size() << ' ';
		for (auto itr=aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			real aijk = itr->second;
			fout << i+1 << ' ' << j+1 << ' ' << k+1 << ' ' << aijk << ' ';
		}
		fout << std::endl;

		fout << aijkl.size() << ' ';
		for (auto itr=aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			real aijkl = itr->second;
			fout << i+1 << ' ' << j+1 << ' ' << k+1 << ' ' << l+1 << ' ' << aijkl << ' ';
		}
		fout << std::endl;
	}

	template<typename real>
	void PseudoBoolean<real>::clear()
	{
		constant = 0;
		ai.clear();
		aij.clear();
		aijk.clear();
		aijkl.clear();
	}

	template<typename real>
	int  PseudoBoolean<real>::nvars() const
	{
		int nVars = -1;
		// Here it is important that all indices appear as pairs 
		// in aij or ai
		for (auto itr=aij.begin(); itr != aij.end(); ++itr) {
			nVars = std::max(nVars, get_i(itr->first));
			nVars = std::max(nVars, get_j(itr->first));
		}
		for (auto itr=ai.begin(); itr != ai.end(); ++itr) {
			nVars = std::max(nVars, itr->first);
		}
		nVars++;
		return nVars;
	}

	template<typename real>
	void PseudoBoolean<real>::print_helper(std::ostream& out) const
	{
		const PseudoBoolean<real>& pbf = *this;

		size_t limit = 7;
		size_t size = pbf.ai.size();

		if (pbf.constant != 0) {
			out << pbf.constant << " ";
		}

		size_t c=0;
		for (auto itr = pbf.ai.begin(); itr != pbf.ai.end(); ++itr) {
			if (c<limit || c==size-1) {
				if (itr->second != 0) {
					if (c==0 && itr->second<0) {
						out << "- ";
					}
					else if (c==0 && pbf.constant != 0) {
						out << "+ ";
					}
					else if (c>0) {
						if (itr->second > 0) out << "+ ";
						else out << "- ";
					}

					if (std::abs(std::abs(itr->second) - 1) > 1e-7) {
						out << std::abs(itr->second) << '*';
					}
					out << "x"<<itr->first<<" "; 
				}
			}
			else if (c == limit) {
				out << "+ ... ";
			}

			c++;
		}

		out << std::endl;

		limit = 5;
		size = pbf.aij.size();
		c=0;
		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			if (c<limit || c==size-1) {
				if (itr->second != 0) {
					if (itr->second > 0) out << "+ ";
					else out << "- ";

					if (std::abs(std::abs(itr->second) - 1) > 1e-7) {
						out << std::abs(itr->second) << '*';
					}
					out << "x"<<i<<"*x"<<j<<" "; 
				}
			}
			else if (c == limit) {
				out << "+ ... ";
			}

			c++;
		}

		out << std::endl;

		limit = 4;
		size = pbf.aijk.size();
		c=0;
		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			if (c<limit || c==size-1) {
				if (itr->second != 0) {
					if (itr->second > 0) out << "+ ";
					else out << "- ";

					if (std::abs(std::abs(itr->second) - 1) > 1e-7) {
						out << std::abs(itr->second) << '*';
					}
					out << "x"<<i<<"*x"<<j<<"*x"<<k<<" "; 
				}
			}
			else if (c == limit) {
				out << "+ ... ";
			}

			c++;
		}

		out << std::endl;

		limit = 2;
		size = pbf.aijkl.size();
		c=0;
		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			if (c<limit || c==size-1) {
				if (itr->second != 0) {
					if (itr->second > 0) out << "+ ";
					else out << "- ";
					if (std::abs(std::abs(itr->second) - 1) > 1e-7) {
						out << std::abs(itr->second) << '*';
					}
					out << "x"<<i<<"*x"<<j<<"*x"<<k<<"*x"<<l<<" "; 
				}
			}
			else if (c == limit) {
				out << "+ ... ";
			}

			c++;
		}

	}

	template<typename real>
	void SymmetricPseudoBoolean<real>::print_helper(std::ostream& out) 
	{
		size_t limit = 7;
		size_t size = bi.size();

		if (constant != 0) {
			out << "C = " << constant << " ";
		}

		size_t c=0;
		for (auto itr = bi.begin(); itr != bi.end(); ++itr) {
			if (c<limit || c==size-1) {
				out << "b(" << itr->first << ") = " << itr->second << "  ";
			}
			else if (c==limit) {
				out << " ... ";
			}
			c++;
		}
		out << std::endl;
		c =0;
		for (auto itr = bij.begin(); itr != bij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			if (c<limit || c==size-1) {
				out << "{b,c}(" << i << "," << j << ") = {" 
					<< bij[itr->first]  << "," << cij[itr->first]
					<< "}  ";
			}
			else if (c==limit) {
				out << " ... ";
			}
			c++;
		}
		out << std::endl;
		c=0;
		for (auto itr = bijk.begin(); itr != bijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			if (c<limit || c==size-1) {
				out << "{b,c,d,e}(" << i << "," << j << "," << k << ") = {" 
					<< bijk[itr->first] << "," << cijk[itr->first] << ","
					<< dijk[itr->first] << "," << eijk[itr->first] 
					<< "}  ";
			}
			else if (c==limit) {
				out << " ... ";
			}
			c++;
		}
		out << std::endl;
		c=0;
		for (auto itr = bijkl.begin(); itr != bijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			if (c<limit || c==size-1) {
				out << "{b,c,d,e,p,q,r,s}(" << i << "," << j << "," << k << "," << l <<") = {" 
					<< bijkl[itr->first] << "," << cijkl[itr->first] << ","
					<< dijkl[itr->first] << "," << eijkl[itr->first] << ","
					<< pijkl[itr->first] << "," << qijkl[itr->first] << ","
					<< rijkl[itr->first] << "," << sijkl[itr->first] 
					<< "}  ";
			}
			else if (c==limit) {
				out << " ... ";
			}
			c++;
		}
	}


	pair make_pair(int i, int j)
	{
		//
		// We only want pairs where i<j, so we
		// permute accordingly
		//
		if (j<i) {
			return make_pair(j,i);
		}
		ASSERT(i<j);
		return std::make_pair(i,j);
	}

	triple make_triple(int i, int j, int k)
	{
		//
		// We only want triplets where i<j<k, so we
		// permute accordingly
		//
		if (j<i) {
			return make_triple(j,i,k);
		}
		if (k<i) {
			return make_triple(k,j,i);
		}
		if (k<j) {
			return make_triple(i,k,j);
		}
		ASSERT(i<j && j<k);


		return std::make_pair(i, std::make_pair(j,k));
	}

	quad make_quad(int i, int j, int k, int l)
	{
		//
		// We only want quads where i<j<k<l, so we
		// permute accordingly
		//
		if (j<i) {
			return make_quad(j,i,k,l);
		}
		if (k<i) {
			return make_quad(k,j,i,l);
		}
		if (l<i) {
			return make_quad(l,j,k,i);
		}
		if (k<j) {
			return make_quad(i,k,j,l);
		}
		if (l<j) {
			return make_quad(i,l,k,j);
		}
		if (l<k) {
			return make_quad(i,j,l,k);
		}
		ASSERT(i<j && j<k && k<l);
		return std::make_pair(std::make_pair(i,j), std::make_pair(k,l));
		//return std::make_tuple<int,int,int,int>(i,j,k,l);
	}

	int get_i(const pair& p)
	{
		return p.first;
	}
	int get_j(const pair& p)
	{
		return p.second;
	}

	int get_i(const triple& t)
	{
		return t.first;
	}
	int get_j(const triple& t)
	{
		return t.second.first;
	}
	int get_k(const triple& t)
	{
		return t.second.second;
	}

	int get_i(const quad& q)
	{
		return q.first.first;
		//return std::get<0>(q);
	}
	int get_j(const quad& q)
	{
		return q.first.second;
		//return std::get<1>(q);
	}
	int get_k(const quad& q)
	{
		return q.second.first;
		//return std::get<2>(q);
	}
	int get_l(const quad& q)
	{
		return q.second.second;
		//return std::get<3>(q);
	}


	//
	// Monomials
	//
	template<typename real>
	void PseudoBoolean<real>::add_monomial(int i, real a)
	{
		ASSERT(i>=0);
		ai[i] += a;
	}
	template<typename real>
	void PseudoBoolean<real>::add_monomial(int i, int j, real a)
	{
		ASSERT(0<=i && 0<=j);
		aij[ make_pair(i,j) ] += a;
	}
	template<typename real>
	void PseudoBoolean<real>::add_monomial(int i, int j, int k, real a)
	{
		ASSERT(0<=i && 0<=j && 0<=k);
		aijk[ make_triple(i,j,k) ] += a;
		//It is important to also add the lower order monomials
		add_monomial(i,j, 0);
		add_monomial(i,k, 0);
		add_monomial(j,k, 0);
	}
	template<typename real>
	void PseudoBoolean<real>::add_monomial(int i, int j, int k, int l, real a)
	{
		ASSERT(0<=i && 0<=j && 0<=k && 0<=l);
		aijkl[ make_quad(i,j,k,l) ] += a;
		//It is important to also add the lower order monomials
		add_monomial(i,j,k, 0);
		add_monomial(i,j,l, 0);
		add_monomial(i,k,l, 0);
		add_monomial(j,k,l, 0);
	}

	//
	// Order-1 clique
	//
	template<typename real>
	void PseudoBoolean<real>::add_clique(int i, real E0, real E1)
	{
		ASSERT(0<=i);

		this->constant += E0;
		this->ai[i] += E1 - E0;
	}

	//
	// Order-2 clique
	//
	template<typename real>
	void PseudoBoolean<real>::add_clique(int i, int j, real E00, real E01, real E10, real E11)
	{
		ASSERT(0<=i && 0<=j);

		this->constant += E00;
		real ai = E10 - E00;
		real aj = E01 - E00;
		real aij = E11 - ai - aj - E00; 
		this->ai[i] += ai;
		this->ai[j] += aj;
		this->aij[make_pair(i,j)] += aij;
	}

	//
	// Order-3 clique
	//
	template<typename real>
	void PseudoBoolean<real>::add_clique(int i, int j, int k, const vector<real>& E)
	{
		add_clique(i,j,k, E.at(0), E.at(1), E.at(2), E.at(3), E.at(4),
		                  E.at(5), E.at(6), E.at(7));
	}
	template<typename real>
	void PseudoBoolean<real>::add_clique(int i, int j, int k, real E000, real E001, real E010, real E011,
	                                                          real E100, real E101, real E110, real E111)
	{
		ASSERT(0<=i && 0<=j && 0<=k);

		this->constant += E000;

		real ai = E100 - E000;
		real aj = E010 - E000;
		real ak = E001 - E000;

		real aij = E110 - ai - aj - E000;
		real aik = E101 - ai - ak - E000;
		real ajk = E011 - aj - ak - E000;

		real aijk = E111 - aij - aik - ajk - ai - aj - ak - E000;

		this->ai[i] += ai;
		this->ai[j] += aj;
		this->ai[k] += ak;
		this->aij[make_pair(i,j)] += aij;
		this->aij[make_pair(i,k)] += aik;
		this->aij[make_pair(j,k)] += ajk;
		this->aijk[make_triple(i,j,k)] += aijk;
	}

	//
	// Order 4 clique
	//
	template<typename real>
	void PseudoBoolean<real>::add_clique(int i, int j, int k, int l, const vector<real>& E)
	{
		add_clique(i,j,k,l, E.at(0), E.at(1), E.at(2), E.at(3), E.at(4),
		                    E.at(5), E.at(6), E.at(7), E.at(8), E.at(9),
		                    E.at(10), E.at(11), E.at(12), E.at(13), E.at(14), E.at(15));
	}

	template<typename real>
	void PseudoBoolean<real>::add_clique(int i, int j, int k, int l, real E0000, real E0001, real E0010, real E0011,
	                                                                 real E0100, real E0101, real E0110, real E0111,
	                                                                 real E1000, real E1001, real E1010, real E1011,
	                                                                 real E1100, real E1101, real E1110, real E1111)
	{
		ASSERT(0<=i && 0<=j && 0<=k && 0<=l);
		ASSERT(0<=i);

		this->constant += E0000;

		real ai = E1000 - E0000;
		real aj = E0100 - E0000;
		real ak = E0010 - E0000;
		real al = E0001 - E0000;

		real aij = E1100 - ai - aj - E0000;
		real aik = E1010 - ai - ak - E0000;
		real ail = E1001 - ai - al - E0000;
		real ajk = E0110 - aj - ak - E0000;
		real ajl = E0101 - aj - al - E0000;
		real akl = E0011 - ak - al - E0000;

		real aijk = E1110 - aij - aik - ajk - ai - aj - ak - E0000;
		real aijl = E1101 - aij - ail - ajl - ai - aj - al - E0000;
		real aikl = E1011 - aik - ail - akl - ai - ak - al - E0000;
		real ajkl = E0111 - ajk - ajl - akl - aj - ak - al - E0000;

		real aijkl = E1111 - aijk - aijl - aikl - ajkl - aij - aik - ail - ajk - ajl - akl - ai - aj - ak - al - E0000;

		this->ai[i] += ai;
		this->ai[j] += aj;
		this->ai[k] += ak;
		this->ai[l] += al;

		this->aij[make_pair(i,j)] += aij;
		this->aij[make_pair(i,k)] += aik;
		this->aij[make_pair(i,l)] += ail;
		this->aij[make_pair(j,k)] += ajk;
		this->aij[make_pair(j,l)] += ajl;
		this->aij[make_pair(k,l)] += akl;

		this->aijk[make_triple(i,j,k)] += aijk;
		this->aijk[make_triple(i,j,l)] += aijl;
		this->aijk[make_triple(i,k,l)] += aikl;
		this->aijk[make_triple(j,k,l)] += ajkl;

		this->aijkl[make_quad(i,j,k,l)] += aijkl;
	}

	

	//
	// Evaluation
	//
	template<typename real>
	real PseudoBoolean<real>::eval(const vector<label>& x) const
	{
		real val = constant;
		for (auto itr = ai.begin(); itr != ai.end(); ++itr) {
			val += itr->second * x.at( itr->first );
		}
		for (auto itr = aij.begin(); itr != aij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			val += itr->second * x.at(i) * x.at(j);
		}
		for (auto itr = aijk.begin(); itr != aijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			val += itr->second * x.at(i) * x.at(j) * x.at(k);
		}
		for (auto itr = aijkl.begin(); itr != aijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);
			val += itr->second * x.at(i) * x.at(j) * x.at(k) * x.at(l);
		}

		return val;
	}



	template<typename real>
	SymmetricPseudoBoolean<real>::SymmetricPseudoBoolean()
	{
		nlpvars = 0;
		constant = 0;
	}

	template<typename real>
	void SymmetricPseudoBoolean<real>::clear()
	{
		nlpvars = 0;
		constant = 0;

		bi.clear();
		bij.clear();
		cij.clear();
		bijk.clear();
		cijk.clear();
		dijk.clear();
		eijk.clear();
		bijkl.clear();
		cijkl.clear();
		dijkl.clear();
		eijkl.clear();
		pijkl.clear();
		qijkl.clear();
		rijkl.clear();
		sijkl.clear();

		indbi.clear();
		indbij.clear();
		indcij.clear();
		indbijk.clear();
		indcijk.clear();
		inddijk.clear();
		indeijk.clear();
		indbijkl.clear();
		indcijkl.clear();
		inddijkl.clear();
		indeijkl.clear();
		indpijkl.clear();
		indqijkl.clear();
		indrijkl.clear();
		indsijkl.clear();			
	}

	template<typename real>
	real SymmetricPseudoBoolean<real>::eval(const vector<label>& x, const vector<label>& y) const
	{
		real val = 0;
		// bi
		for (auto itr = bi.begin(); itr != bi.end(); ++itr) {
			int i = itr->first;

			val += itr->second * (x.at(i) +  (1-y.at(i)));
		}
		// bij
		for (auto itr = bij.begin(); itr != bij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);

			val += itr->second * (x.at(i) * x.at(j)  +  (1-y.at(i)) * (1-y.at(j)));
		}
		// cij
		for (auto itr = cij.begin(); itr != cij.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);

			val += itr->second * (x.at(i) * (1-y.at(j))  +  (1-y.at(i)) * x.at(j));
		}
		

		// bijk
		for (auto itr = bijk.begin(); itr != bijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);

			val += itr->second * (  x.at(i)     * x.at(j)     * x.at(k)      +  
			                        (1-y.at(i)) * (1-y.at(j)) * (1-y.at(k))  );
		}
		// cijk
		for (auto itr = cijk.begin(); itr != cijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);

			val += itr->second * (  x.at(i)     * x.at(j)     * (1-y.at(k))      +  
			                        (1-y.at(i)) * (1-y.at(j)) * x.at(k)  );
		}
		// dijk
		for (auto itr = dijk.begin(); itr != dijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);

			val += itr->second * (  x.at(i)     * (1-y.at(j))  * x.at(k)      +  
			                        (1-y.at(i)) * x.at(j)      * (1-y.at(k))  );
		}
		// eijk
		for (auto itr = eijk.begin(); itr != eijk.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);

			val += itr->second * (  (1-y.at(i)) * x.at(j)     * x.at(k)      +  
			                        x.at(i)     * (1-y.at(j)) * (1-y.at(k))  );
		}


		// bijkl
		for (auto itr = bijkl.begin(); itr != bijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * x.at(j)     * x.at(k)     * x.at(l)     +  
			                        (1-y.at(i)) * (1-y.at(j)) * (1-y.at(k)) * (1-y.at(l)) );
		}
		// cijkl
		for (auto itr = cijkl.begin(); itr != cijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * x.at(j)     * x.at(k)     * (1-y.at(l))    +  
			                        (1-y.at(i)) * (1-y.at(j)) * (1-y.at(k)) * x.at(l) );
		}
		// dijkl
		for (auto itr = dijkl.begin(); itr != dijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * x.at(j)     * (1-y.at(k)) * x.at(l)     +  
			                        (1-y.at(i)) * (1-y.at(j)) * x.at(k)     * (1-y.at(l)) );
		}
		// eijkl
		for (auto itr = eijkl.begin(); itr != eijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * (1-y.at(j)) * x.at(k)     * x.at(l)     +  
			                        (1-y.at(i)) * x.at(j)     * (1-y.at(k)) * (1-y.at(l)) );
		}
		// pijkl
		for (auto itr = pijkl.begin(); itr != pijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  (1-y.at(i)) * x.at(j)     * x.at(k)     * x.at(l)     +  
			                        x.at(i)     * (1-y.at(j)) * (1-y.at(k)) * (1-y.at(l)) );
		}
		// qijkl
		for (auto itr = qijkl.begin(); itr != qijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * x.at(j)     * (1-y.at(k)) * (1-y.at(l)) +  
			                        (1-y.at(i)) * (1-y.at(j)) * x.at(k)     *  x.at(l)    );
		}
		// rijkl
		for (auto itr = rijkl.begin(); itr != rijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * (1-y.at(j)) * x.at(k)     * (1-y.at(l))  +  
			                        (1-y.at(i)) * x.at(j)     * (1-y.at(k)) * x.at(l)      );
		}
		// sijkl
		for (auto itr = sijkl.begin(); itr != sijkl.end(); ++itr) {
			int i = get_i(itr->first);
			int j = get_j(itr->first);
			int k = get_k(itr->first);
			int l = get_l(itr->first);

			val += itr->second * (  x.at(i)     * (1-y.at(j)) * (1-y.at(k)) * x.at(l)     +  
			                        (1-y.at(i)) * x.at(j)     * x.at(k)     * (1-y.at(l)) );
		}

		return constant + val/real(2);
	}

	//////////////////////////////
	// Functions to get indices //
	//////////////////////////////
	//int SymmetricPseudoBoolean::ib(int i) { return getindex(indbi, i); }
	template<typename real> int SymmetricPseudoBoolean<real>::ib(int i, int j){ return getindex(indbij, make_pair(i,j)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ic(int i, int j){ return getindex(indcij, make_pair(i,j)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ib(int i, int j, int k){ return getindex(indbijk, make_triple(i,j,k)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ic(int i, int j, int k){ return getindex(indcijk, make_triple(i,j,k)); }
	template<typename real> int SymmetricPseudoBoolean<real>::id(int i, int j, int k){ return getindex(inddijk, make_triple(i,j,k)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ie(int i, int j, int k){ return getindex(indeijk, make_triple(i,j,k)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ib(int i, int j, int k, int l){ return getindex(indbijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ic(int i, int j, int k, int l){ return getindex(indcijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::id(int i, int j, int k, int l){ return getindex(inddijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ie(int i, int j, int k, int l){ return getindex(indeijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ip(int i, int j, int k, int l){ return getindex(indpijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::iq(int i, int j, int k, int l){ return getindex(indqijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::ir(int i, int j, int k, int l){ return getindex(indrijkl, make_quad(i,j,k,l)); }
	template<typename real> int SymmetricPseudoBoolean<real>::is(int i, int j, int k, int l){ return getindex(indsijkl, make_quad(i,j,k,l)); }



	template<typename Map>
	int maxindextemp(const Map& m)
	{
		int maxval = -1;
		for (auto itr = m.begin(); itr!= m.end(); ++itr) {
			if (itr->second > maxval) {
				maxval = itr->second;
			}
		}
		return maxval;
	}

	template<typename real>
	void SymmetricPseudoBoolean<real>::test()
	{

	}

	template<typename real>
	void SymmetricPseudoBoolean<real>::test_clear()
	{
		clear();

		int b123 = ib(1,2,3);
		int b2 = ib(1,2);
		int e3415 = ie(2,3,4,5);
		int p3415 = ip(2,3,4,5);
		int q3415 = iq(2,3,4,5);
		int b0 = ib(0,1);
		b2 = ib(1,2);
		int b45 = ib(4,5);
		int c45 = ic(4,5);
		int c345 = ic(3,4,5);
		int d345 = id(3,4,5);
		int e345 = ie(3,4,5);
		int b3415 = ib(1,3,4,5);
		int c3415 = ic(1,3,4,5);
		int d3415 = id(1,3,4,5);
		int r3415 = ir(1,3,4,5);
		int s3415 = is(1,3,4,5);

		ASSERT(b123 == ib(1,2,3));
		ASSERT(b2 == ib(1,2));
		ASSERT(e3415 == ie(2,3,4,5));
		ASSERT( p3415 == ip(2,3,4,5));
		ASSERT( q3415 == iq(2,3,4,5));
		ASSERT(  b0 == ib(0,1));
		ASSERT(  b45 == ib(4,5));
		ASSERT(  c45 == ic(4,5));
		ASSERT(  c345 == ic(3,4,5));
		ASSERT(  d345 == id(3,4,5));
		ASSERT(  e345 == ie(3,4,5));
		ASSERT(  b3415 == ib(1,3,4,5));
		ASSERT(  c3415 == ic(1,3,4,5));
		ASSERT(  d3415 == id(1,3,4,5));
		ASSERT(  r3415 == ir(1,3,4,5));
		ASSERT(  s3415 == is(1,3,4,5));

		ASSERT( b123 == 0);
		ASSERT( b2 != b123 );
		ASSERT( b45 != c45 );
		ASSERT( r3415 != s3415);
		ASSERT( r3415 != d3415 );


		int maxval = -1;
		maxval = std::max(maxval, maxindextemp(indbi));
		maxval = std::max(maxval, maxindextemp(indbij));
		maxval = std::max(maxval, maxindextemp(indcij));
		maxval = std::max(maxval, maxindextemp(indbijk));
		maxval = std::max(maxval, maxindextemp(indcijk));
		maxval = std::max(maxval, maxindextemp(inddijk));
		maxval = std::max(maxval, maxindextemp(indeijk));
		maxval = std::max(maxval, maxindextemp(indbijkl));
		maxval = std::max(maxval, maxindextemp(indcijkl));
		maxval = std::max(maxval, maxindextemp(inddijkl));
		maxval = std::max(maxval, maxindextemp(indeijkl));
		maxval = std::max(maxval, maxindextemp(indpijkl));
		maxval = std::max(maxval, maxindextemp(indqijkl));
		maxval = std::max(maxval, maxindextemp(indrijkl));
		maxval = std::max(maxval, maxindextemp(indsijkl));
		maxval++;
		ASSERT( maxval == 16 );
		ASSERT( maxval == nlpvars );
	}
}


#include "pb_instances.inc"
