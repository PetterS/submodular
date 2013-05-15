
#include <fstream>
#include <memory>
#include <typeinfo>

// BK maxflow
#include "graph.h"
// LP solver
#include <coin/ClpSimplex.hpp>

#include "PseudoBoolean.h"
#include "Minimizer.h"

#include "StdAfx.h"
#include "apgc/APGC.h"
#include "prgc/PRGC.h"
#include "../../kod/insert_clique2.h"
#include "Petter-Color.h"

using namespace std;

namespace Petter
{

	size_t read_vec(vector<int>& vec, ifstream& fin, int n)
	{
		size_t cnt=0;
		vec.resize(n);
		for (int i=0;i<n;++i) {
			fin >> vec[i];
			if (vec[i]!=0){
				cnt++;
			}
		}
		return cnt;
	}

	template<typename real, typename vectype>
	void read_reduced_polynomial(vector<vectype>& polynomial, ifstream& fin)
	{
		real c;
		int n,i,j;
		typedef typename Generators<real>::Monomial Monomial;

		// Start with an empty polynomial
		polynomial.clear();

		// Read constant term
		fin >> c;
		if (c!=0) {
			polynomial.push_back( Monomial(c) );
		}
		ASSERT(fin);

		// Read linear terms
		fin >> n;
		ASSERT(fin);
		for (i=0; i<n; ++i) {
			fin >> c;
			ASSERT(fin);
			if (c!=0) {
				polynomial.push_back( Monomial(i,c) );
			}
		}

		// Read quadratic terms
		fin >> n;
		for (int ind=0;ind<n;++ind) {
			fin >> i >> j >> c;
			i--; // Convert between Maple indexing and C++ indexing
			j--; // Convert between Maple indexing and C++ indexing
			ASSERT(fin);
			ASSERT(i>=0 && j>=0);
			if (c!=0) {
				polynomial.push_back( Monomial(i,j,c) );
			}
		}
	}

	// Creates a generator of a specified degree and
	// reads it from the input stream.
	template<typename real>
	SymmetricGenerator<real>::SymmetricGenerator(int degree, std::ifstream& fin)
	{
		ASSERT(fin);
		ASSERT(degree == 2 || degree == 3 || degree == 4);
		for (int i = 0; i < degree; ++i) {
			int tmp;
			fin >> tmp;
			tmp--; // Convert between Maple indexing and C++ indexing
			this->indices1.push_back(tmp);
		}

		int nvalues = 1 << degree;
		for (int i = 0; i < nvalues; ++i) {
			real tmp;
			fin >> tmp;
			this->values1.push_back(tmp);
		}

		for (int i = 0; i < degree; ++i) {
			int tmp;
			fin >> tmp;
			tmp--; // Convert between Maple indexing and C++ indexing
			this->indices2.push_back(tmp);
		}

		for (int i = 0; i < nvalues; ++i) {
			real tmp;
			fin >> tmp;
			this->values2.push_back(tmp);
		}
		ASSERT_STR(fin, "Could not read generator.");
	}

	template<typename real>
	Generators<real>::Generators(const string& filename)
	{
		ifstream fin(filename);

		// First try to read the version of the file.
		string V;
		version = 0;
		fin >> V >> version;
		if (V != "Version") {
			// There was no version string. Reset the file stream.
			fin.close();
			fin.open(filename);
			version = 1;
		}

		ASSERT_STR(version == 1 || version == 2 || version == 3, "Invalid version number.");

		fin >> ngen2 >> ngen3 >> ngen4pos;
		ASSERT_STR(fin, "Could not read file");

		ngen4=2*ngen4pos; //Half of generators required to be positive!
		ngen4neg = ngen4pos;

		std::cout << "\nGenerator version : " << version << "\n";

		if (version == 1) {
			//
			// Read reduced forms
			//

			// Read quadratic generators
			gen2red.resize(ngen2);
			for (int i=0;i<ngen2;++i) {
				read_reduced_polynomial<real>( gen2red.at(i), fin);
			}

			// Read cubic generators
			gen3red.resize(ngen3);
			for (int i=0;i<ngen3;++i) {
				read_reduced_polynomial<real>( gen3red.at(i), fin);
			}

			// Read quadratic generators
			gen4redpos.resize(ngen4pos);
			for (int i=0;i<ngen4pos;++i) {
				read_reduced_polynomial<real>( gen4redpos.at(i), fin);
			}
			gen4redneg.resize(ngen4neg);
			for (int i=0;i<ngen4neg;++i) {
				read_reduced_polynomial<real>( gen4redneg.at(i), fin);
			}
		}
		if (version == 2) {
			//
			// Read actual generator.
			//
			for (int i=0;i<ngen2;++i) {
				SymmetricGenerator<real> generator(2, fin);
				this->gen2.push_back(generator);
			}

			for (int i=0;i<ngen3;++i) {
				SymmetricGenerator<real> generator(3, fin);
				this->gen3.push_back(generator);
			}

			for (int i=0;i<ngen4pos;++i) {
				SymmetricGenerator<real> generator(4, fin);
				this->gen4pos.push_back(generator);
			}

			for (int i=0;i<ngen4neg;++i) {
				SymmetricGenerator<real> generator(4, fin);
				this->gen4neg.push_back(generator);
			}
		}
		if (version == 3) {
			//
			// Read actual generator.
			//
			for (int i=0;i<ngen2;++i) {
				SymmetricGenerator<real> generator(2, fin);
				this->gen2.push_back(generator);
			}

			for (int i=0;i<ngen3;++i) {
				SymmetricGenerator<real> generator(3, fin);
				this->gen3.push_back(generator);
			}

			for (int i=0;i<ngen4pos;++i) {
				SymmetricGenerator<real> generator(4, fin);
				this->gen4pos.push_back(generator);
			}

			for (int i=0;i<ngen4neg;++i) {
				SymmetricGenerator<real> generator(4, fin);
				this->gen4neg.push_back(generator);
			}
		}




		//
		// Read aa,bb and cc
		//

		nentries2 = read_vec(cc,fin,ngen2); ASSERT(nentries2 == 2);
		read_vec(obj2,fin,ngen2);

		ASSERT_STR(fin, "Could not read cc");

		nentries3  = read_vec(bb12,fin,ngen3); ASSERT(nentries3 == 4);
		nentries3 += read_vec(bb13,fin,ngen3); ASSERT(nentries3 == 8);
		nentries3 += read_vec(bb23,fin,ngen3); ASSERT(nentries3 == 12);
		nentries3 += read_vec(bb123,fin,ngen3); ASSERT(nentries3 == 20);

		read_vec(obj3,fin,ngen3);


		ASSERT_STR(fin, "Could not read bb");

		nentries4  = read_vec(aa12pos,fin,ngen4pos);
		nentries4 += read_vec(aa13pos,fin,ngen4pos);
		nentries4 += read_vec(aa14pos,fin,ngen4pos);
		nentries4 += read_vec(aa23pos,fin,ngen4pos);
		nentries4 += read_vec(aa24pos,fin,ngen4pos);
		nentries4 += read_vec(aa34pos,fin,ngen4pos);
		nentries4 += read_vec(aa123pos,fin,ngen4pos);
		nentries4 += read_vec(aa124pos,fin,ngen4pos);
		nentries4 += read_vec(aa134pos,fin,ngen4pos);
		nentries4 += read_vec(aa234pos,fin,ngen4pos);
		nentries4 += read_vec(aa1234pos,fin,ngen4pos);

		read_vec(obj4pos,fin,ngen4pos);

		size_t ndummy = read_vec(aa12neg,fin,ngen4neg);
		ndummy += read_vec(aa13neg,fin,ngen4neg);
		ndummy += read_vec(aa14neg,fin,ngen4neg);
		ndummy += read_vec(aa23neg,fin,ngen4neg);
		ndummy += read_vec(aa24neg,fin,ngen4neg);
		ndummy += read_vec(aa34neg,fin,ngen4neg);
		ndummy += read_vec(aa123neg,fin,ngen4neg);
		ndummy += read_vec(aa124neg,fin,ngen4neg);
		ndummy += read_vec(aa134neg,fin,ngen4neg);
		ndummy += read_vec(aa234neg,fin,ngen4neg);
		ndummy += read_vec(aa1234neg,fin,ngen4neg);

		read_vec(obj4neg,fin,ngen4neg);

		ASSERT_STR(ndummy==nentries4,"Not same number of entries in pos and neg");
		ASSERT_STR(fin, "Could not read aa");

		// The file should now be empty.
		int tmp;
		fin >> tmp;
		ASSERT_STR( !fin, "There is more information in the file.");
	}


	template<typename real>
	GeneratorPseudoBoolean<real>::GeneratorPseudoBoolean(const Generators<real>& generators) :
		gen(generators)
	{
		constant = 0;
		nlpvars = 0;
	}

	template<typename real>
	void GeneratorPseudoBoolean<real>::clear()
	{
		nlpvars = 0;

		constant = 0;
		alphai.clear();
		alphaij.clear();
		alphaijk.clear();
		alphaijkl.clear();

		indaa.clear();
		indbb.clear();
		indcc.clear();

		var_used.clear();
	}









	//////////////////////////////
	// Functions to get indices //
	//////////////////////////////
	template<typename real> int GeneratorPseudoBoolean<real>::icc(int i, int j){ return getindex(indcc, make_pair(i,j)); }
	template<typename real> int GeneratorPseudoBoolean<real>::ibb(int i, int j, int k){ return getindex(indbb, make_triple(i,j,k)); }
	template<typename real> int GeneratorPseudoBoolean<real>::iaa(int i, int j, int k, int l){ return getindex(indaa, make_quad(i,j,k,l)); }



	//
	// The default behaviour is to raise an error message
	//
	template<typename real>
	void GeneratorPseudoBoolean<real>::create_lp(const PseudoBoolean<real>& pbf)
	{
		auto& id = typeid(real);
		std::string msg = "Linear programming not available for type ";
		msg += id.name();
		throw std::runtime_error(msg);
	}

	//
	// For doubles, LP is available
	//
	template<>
	void GeneratorPseudoBoolean<double>::create_lp(const PseudoBoolean<double>& pbf)
	{
		// Convenient to have this name available
		typedef double real;
		clear();

		// Save the constant
		constant = pbf.constant;

		// Save the linear coefficients
		real linear_coef_sum = 0;
		for (auto itr = pbf.ai.begin(); itr != pbf.ai.end(); ++itr) {
			alphai[itr->first] =  itr->second / 2;
			linear_coef_sum    += itr->second / 2;
			var_used[itr->first] = true;
		}


		int nLPVars = int( gen.ngen2*pbf.aij.size() + gen.ngen3*pbf.aijk.size() + gen.ngen4pos*pbf.aijkl.size() );

		int nConstraints = int( pbf.aij.size() + pbf.aijk.size() + pbf.aijkl.size() ); // equality constraints


		// Compute the number of entries in the
		// constraint matrix
		size_t nEntries = gen.nentries2*pbf.aij.size() + gen.nentries3*pbf.aijk.size() + gen.nentries4*pbf.aijkl.size();

		//Description of sparse matrix
		vector<int> rows;
		vector<int> cols;
		vector<double> values;
		rows.reserve(nEntries);
		cols.reserve(nEntries);
		values.reserve(nEntries);

		//Other LP parameters
		vector<double> rhs_eq(nConstraints, 0.0); // only equality constraints
		double var_limit = 100000000;
		vector<double> var_lb(nLPVars, 0.0);      // only non-negative variables
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
		auto change_rhs = [&rhs_eq](size_t row, double eq)
		{
			rhs_eq.at(row) = eq;
		};

		//////////////////////////
		// Equality constraints //
		//////////////////////////
		map<pair,int> indpair_to_con;
		map<triple,int> indtriple_to_con;


		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			int i=get_i(itr->first);
			int j=get_j(itr->first);
			real aij = itr->second;

			indpair_to_con[itr->first] = con;
			int ind=icc(i,j); //column index of cc
			for (int ii=0;ii<gen.ngen2;++ii){
				int colind=ind+ii;
				add_element(con, colind, gen.cc.at(ii)); // cc[ii] should always be non-zero, since degree 2 generator
				//Objective function
				cost[colind] = -gen.obj2[ii];
			}
			change_rhs(con, aij/2);
			con++;
			var_used[i] = true;
			var_used[j] = true;
		}

		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			int i=get_i(itr->first);
			int j=get_j(itr->first);
			int k=get_k(itr->first);
			real aijk = itr->second;

			indtriple_to_con[itr->first] = con;

			int ind=ibb(i,j,k);
			// ij,ik,jk
			int con_12 = indpair_to_con[make_pair(i,j)];
			int con_13 = indpair_to_con[make_pair(i,k)];
			int con_23 = indpair_to_con[make_pair(j,k)];

			for (int ii=0;ii<gen.ngen3;++ii){
				int colind=ind+ii;
				add_element(con, colind, gen.bb123[ii]); // bb123[ii] should always be non-zero, since degree 3 generator
				if (gen.bb12[ii]!=0){
					add_element(con_12, colind, gen.bb12[ii]);
				}
				if (gen.bb13[ii]!=0){
					add_element(con_13, colind, gen.bb13[ii]);
				}
				if (gen.bb23[ii]!=0){
					add_element(con_23, colind, gen.bb23[ii]);
				}
				//Objective function
				cost[colind] = -gen.obj3[ii];
			}
			change_rhs(con, aijk/2);
			con++;

			var_used[i] = true;
			var_used[j] = true;
			var_used[k] = true;
		}

		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			int i=get_i(itr->first);
			int j=get_j(itr->first);
			int k=get_k(itr->first);
			int l=get_l(itr->first);
			real aijkl = itr->second;

			int ind=iaa(i,j,k,l);

			// 12,13,14,23,24,34
			int con_12 = indpair_to_con[make_pair(i,j)];
			int con_13 = indpair_to_con[make_pair(i,k)];
			int con_14 = indpair_to_con[make_pair(i,l)];
			int con_23 = indpair_to_con[make_pair(j,k)];
			int con_24 = indpair_to_con[make_pair(j,l)];
			int con_34 = indpair_to_con[make_pair(k,l)];

			int con_123 = indtriple_to_con[make_triple(i,j,k)];
			int con_124 = indtriple_to_con[make_triple(i,j,l)];
			int con_134 = indtriple_to_con[make_triple(i,k,l)];
			int con_234 = indtriple_to_con[make_triple(j,k,l)];

			if (aijkl>=0){

				// Save that we are using positive generators
				posgen4[itr->first] = true;

				for (int ii=0;ii<gen.ngen4pos;++ii){
					int colind=ind+ii;
					add_element(con, colind, gen.aa1234pos[ii]); // aa1234[ii] should always be non-zero, since degree 4 generator
					if (gen.aa12pos[ii]!=0){
						add_element(con_12, colind, gen.aa12pos[ii]);
					}
					if (gen.aa13pos[ii]!=0){
						add_element(con_13, colind, gen.aa13pos[ii]);
					}
					if (gen.aa14pos[ii]!=0){
						add_element(con_14, colind, gen.aa14pos[ii]);
					}
					if (gen.aa23pos[ii]!=0){
						add_element(con_23, colind, gen.aa23pos[ii]);
					}
					if (gen.aa24pos[ii]!=0){
						add_element(con_24, colind, gen.aa24pos[ii]);
					}
					if (gen.aa34pos[ii]!=0){
						add_element(con_34, colind, gen.aa34pos[ii]);
					}
					if (gen.aa123pos[ii]!=0){
						add_element(con_123, colind, gen.aa123pos[ii]);
					}
					if (gen.aa124pos[ii]!=0){
						add_element(con_124, colind, gen.aa124pos[ii]);
					}
					if (gen.aa134pos[ii]!=0){
						add_element(con_134, colind, gen.aa134pos[ii]);
					}
					if (gen.aa234pos[ii]!=0){
						add_element(con_234, colind, gen.aa234pos[ii]);
					}
					//Objective function
					cost[colind] = -gen.obj4pos[ii];
				}
			}
			else{

				// Save that we are using negative generators
				posgen4[itr->first] = false;

				for (int ii=0;ii<gen.ngen4neg;++ii){
					int colind=ind+ii;
					add_element(con, colind, gen.aa1234neg[ii]); // aa1234[ii] should always be non-zero, since degree 4 generator
					if (gen.aa12neg[ii]!=0){
						add_element(con_12, colind, gen.aa12neg[ii]);
					}
					if (gen.aa13neg[ii]!=0){
						add_element(con_13, colind, gen.aa13neg[ii]);
					}
					if (gen.aa14neg[ii]!=0){
						add_element(con_14, colind, gen.aa14neg[ii]);
					}
					if (gen.aa23neg[ii]!=0){
						add_element(con_23, colind, gen.aa23neg[ii]);
					}
					if (gen.aa24neg[ii]!=0){
						add_element(con_24, colind, gen.aa24neg[ii]);
					}
					if (gen.aa34neg[ii]!=0){
						add_element(con_34, colind, gen.aa34neg[ii]);
					}
					if (gen.aa123neg[ii]!=0){
						add_element(con_123, colind, gen.aa123neg[ii]);
					}
					if (gen.aa124neg[ii]!=0){
						add_element(con_124, colind, gen.aa124neg[ii]);
					}
					if (gen.aa134neg[ii]!=0){
						add_element(con_134, colind, gen.aa134neg[ii]);
					}
					if (gen.aa234neg[ii]!=0){
						add_element(con_234, colind, gen.aa234neg[ii]);
					}
					//Objective function
					cost[colind] = -gen.obj4neg[ii];
				}
			}
			change_rhs(con, aijkl/2);
			con++;

			var_used[i] = true;
			var_used[j] = true;
			var_used[k] = true;
			var_used[l] = true;
		}

		ASSERT(nLPVars == this->nlpvars);

		// Make sure that the precomputation of space is correct
		ASSERT( rows.size() == nEntries );
		// Make sure that the precomputation of rows is correct
		ASSERT( con == nConstraints );


		CoinPackedMatrix coinMatrix(false,&rows[0],&cols[0],&values[0], CoinBigIndex(values.size()) );
		ClpSimplex lpSolver;
		lpSolver.loadProblem (coinMatrix, &var_lb[0], &var_ub[0], &cost[0], &rhs_eq[0], &rhs_eq[0]);
		lpSolver.setLogLevel(0);

		int error = lpSolver.dual(); //primal, dual, barrier
		if (error != 0) {
			throw std::runtime_error("Clp failed");
		}

		double* lpvars = lpSolver.primalColumnSolution();



		double obj = 0;
		for (int ii=0;ii<nLPVars;++ii){
			obj+= (lpvars[ii]*cost[ii]);
		}

		//cout << endl << "Objective                     : " << -obj << endl;
		//cout <<         "Objective (with linear terms) : " << -obj+linear_coef_sum << endl;

		for (auto itr = pbf.aij.begin(); itr != pbf.aij.end(); ++itr) {
			const pair& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);

			alphaij[ind].resize(gen.ngen2);
			int icctmp = icc(i,j);
			for (int ii=0;ii<gen.ngen2;++ii){
				alphaij[ind][ii] = lpvars[ icctmp+ii ];

				//if (alphaij[ind][ii] != 0) {
				//	// +1 to reflect Maple indexing
				//	cout << "alphaij (" << i+1 << "," << j+1 << "; " << ii+1 << "): " << alphaij[ind][ii] << endl;
				//}
			}
		}


		for (auto itr = pbf.aijk.begin(); itr != pbf.aijk.end(); ++itr) {
			const triple& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);

			alphaijk[ind].resize(gen.ngen3);
			int ibbtmp = ibb(i,j,k);
			for (int ii=0;ii<gen.ngen3;++ii){
				alphaijk[ind][ii] = lpvars[ ibbtmp+ii ];

				//if (alphaijk[ind][ii] != 0) {
				//	// +1 to reflect Maple indexing
				//	cout << "alphaijk (" << i+1 << "," << j+1 << "," << k+1 << "; " << ii+1 << "): " << alphaijk[ind][ii] << endl;
				//}
			}
		}

		for (auto itr = pbf.aijkl.begin(); itr != pbf.aijkl.end(); ++itr) {
			const quad& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);
			int l=get_l(ind);

			alphaijkl[ind].resize(gen.ngen4pos);
			int iaatmp = iaa(i,j,k,l);
			for (int ii=0;ii<gen.ngen4pos;++ii){
				alphaijkl[ind][ii] = lpvars[ iaatmp+ii ];

				//if (alphaijkl[ind][ii] != 0) {
				//	// +1 to reflect Maple indexing
				//	cout << "alphaijkl (" << i+1 << "," << j+1 << "," << k+1 << "," << l+1 << "; " << ii+1 << "): " << alphaijkl[ind][ii] << endl;
				//}
			}
		}




	} //end of create_lp



	template<typename Monomial>
	void print_monomial(std::ostream& out, const Monomial& monomial)
	{
		out << "Monomial ";
		out << monomial.c;
		if (monomial.i>=0) {
			out << "*x" << monomial.i;
		}
		if (monomial.j>=0) {
			out << "*x" << monomial.j;
		}
	}


	//
	// Adds alpha times a monomial to a graph
	//
	template<typename real>
	void add_monomial_to_graph( real& C, Graph<real,real,real>& graph, real alpha, const typename Generators<real>::Monomial& monomial, const vector<int>& idx)
	{
		monomial.check(); //Debug

		//cout << "Adding ";
		//print_monomial(cout,monomial);
		//cout << " to the graph." << endl;

		real beta = monomial.c; // Constant in front of monomial

		if (monomial.j >= 0) {
			// Quadraic term, add to graph
			add_monomial_2_to_graph(C,graph, idx.at(monomial.i) ,idx.at(monomial.j), alpha*beta);
		}
		else if (monomial.i >= 0) {
			// Linear term, add to graph
			add_monomial_1_to_graph(C,graph, idx.at(monomial.i), alpha*beta);
		}
		else {
			// Constant
			C += alpha*beta;
		}
	}

	//
	// Adds alpha times a quadratic generator to a graph
	//
	template<typename real, typename vectype>
	void add_generator_to_graph( int nVars, real& C, Graph<real,real,real>& graph, real alpha, const vector<vectype>& polynomial, int i, int j)
	{
		// Go though the number of indices used
		// Not needed because there will be no extra variables

		vector<int> idx(4); // Translates from "local" indices to "global"
		idx[0] = i; // x variables
		idx[1] = j;
		idx[2] = i + nVars; // y variables
		idx[3] = j + nVars;

		for (auto itr = polynomial.begin(); itr != polynomial.end(); ++itr) {
			add_monomial_to_graph(C,graph,alpha,*itr,idx);
		}
	}


	//
	// Adds alpha times a cubic generator to a graph
	//
	template<typename real, typename vectype>
	void add_generator_to_graph( int nVars, real& C, Graph<real,real,real>& graph, real alpha, const vector<vectype>& polynomial, int i, int j, int k)
	{
		// Go though the number of indices used
		int maxi = -1;
		for (auto itr = polynomial.begin(); itr != polynomial.end(); ++itr) {
			maxi = max( maxi, itr->i);
			maxi = max( maxi, itr->j);
		}

		vector<int> idx(maxi+1); // Translates from "local" indices to "global"
		idx.at(0) = i; // x variables
		idx.at(1) = j;
		idx.at(2) = k;
		idx.at(3) = i + nVars; // y variables
		idx.at(4) = j + nVars;
		idx.at(5) = k + nVars;

		for (int n=6; n<=maxi; ++n) {
			idx.at(n) = graph.add_node(); // Extra variable
		}

		for (auto itr = polynomial.begin(); itr != polynomial.end(); ++itr) {
			add_monomial_to_graph(C,graph,alpha,*itr,idx);
		}
	}

	//
	// Adds alpha times a quartic generator to a graph
	//
	template<typename real, typename vectype>
	void add_generator_to_graph( int nVars, real& C, Graph<real,real,real>& graph, real alpha, const vector<vectype>& polynomial, int i, int j, int k, int l)
	{
		// Go though the number of indices used
		int maxi = -1;
		for (auto itr = polynomial.begin(); itr != polynomial.end(); ++itr) {
			maxi = max( maxi, itr->i);
			maxi = max( maxi, itr->j);
		}

		vector<int> idx(maxi+1); // Translates from "local" indices to "global"
		idx.at(0) = i;  // x variables
		idx.at(1) = j;
		idx.at(2) = k;
		idx.at(3) = l;
		idx.at(4) = i + nVars; // y variables
		idx.at(5) = j + nVars;
		idx.at(6) = k + nVars;
		idx.at(7) = l + nVars;

		for (int n=8; n<=maxi; ++n) {
			idx.at(n) = graph.add_node(); // Extra variable
		}

		for (auto itr = polynomial.begin(); itr != polynomial.end(); ++itr) {
			add_monomial_to_graph(C,graph,alpha,*itr,idx);
		}
	}


	template<typename real>
	real GeneratorPseudoBoolean<real>::minimize(vector<label>& x, int& nlabelled) const
	{
		ASSERT(this->gen.version == 1 || this->gen.version == 2  || this->gen.version == 3);
		if (this->gen.version == 1) {
			return minimize_version_1(x, nlabelled);
		}
		if (this->gen.version == 2) {
			return minimize_version_2(x, nlabelled);
		}
		if (this->gen.version == 3) {
			return minimize_version_3(x, nlabelled);
		}



		return 0;  // Not used.
	}

	float make_clique_positive(int clique_size, float* E)
	{
		float min_value = 0.0f;
		int num_values = 1 << clique_size;
		for (int i = 0; i < num_values; ++i) {
			min_value = std::min(min_value, E[i]);
		}

		for (int i = 0; i < num_values; ++i) {
			E[i] -= min_value;
		}

		return min_value;
	}

	template<typename GC, typename real>
	void add_generator_to_graph(GC* graph, real* C, int i, int j, int k, int l,
		const std::vector<real>& values,
		real alpha, std::unique_ptr<PseudoBoolean<real>>& f_debug)
	{
		ASSERT(values.size() == 16);
		float E[] = {alpha * values.at(0), // E0000
			alpha * values.at(1), // E0001
			alpha * values.at(2), // E0010
			alpha * values.at(3), // E0011
			alpha * values.at(4), // E0100
			alpha * values.at(5), // E0101
			alpha * values.at(6), // E0110
			alpha * values.at(7), // E0111
			alpha * values.at(8), // E1000
			alpha * values.at(9), // E1001
			alpha * values.at(10), // E1010
			alpha * values.at(11), // E1011
			alpha * values.at(12), // E1100
			alpha * values.at(13), // E1101
			alpha * values.at(14), // E1110
			alpha * values.at(15)};// E1111
		int indices[] = {i, j, k, l};
		(*C) += make_clique_positive(4, E);
		graph->AddHigherTerm(indices, E);

		if (f_debug) {
			std::vector<real> Ev(E, E+16);
			f_debug->add_clique(indices[0], indices[1], indices[2], indices[3], Ev);
		}
	}


	//instead of just one, add all combinations of cliques.
	// valuesabc is clique of size 3 with already inserted alpha_abc
	// valuesab is clique of size 2 with already inserted alpha_abc

	template< typename real>
	void add_generators_to_clique(real alpha, float* E,  const std::vector<real>& values1234)
	{
		if (values1234.size() == 16){
			E[0] =  alpha * values1234.at(0); // E0000
			E[1] =  alpha * values1234.at(1); // E0001
			E[2] =  alpha * values1234.at(2); // E0010
			E[3] =  alpha * values1234.at(3); // E0011
			E[4] =  alpha * values1234.at(4); // E0100
			E[5] =  alpha * values1234.at(5); // E0101
			E[6] =  alpha * values1234.at(6); // E0110
			E[7] = alpha * values1234.at(7); // E0111
			E[8] =  alpha * values1234.at(8); // E0000
			E[9] =  alpha * values1234.at(9); // E0001
			E[10] =  alpha * values1234.at(10); // E0010
			E[11] =  alpha * values1234.at(11); // E0011
			E[12] =  alpha * values1234.at(12); // E0100
			E[13] =  alpha * values1234.at(13); // E0101
			E[14] =  alpha * values1234.at(14); // E0110
			E[15] = alpha * values1234.at(15); // E0111
		}
		if (values1234.size() == 8){
			E[0] =  alpha * values1234.at(0); // E0000
			E[1] =  alpha * values1234.at(1); // E0001
			E[2] =  alpha * values1234.at(2); // E0010
			E[3] =  alpha * values1234.at(3); // E0011
			E[4] =  alpha * values1234.at(4); // E0100
			E[5] =  alpha * values1234.at(5); // E0101
			E[6] =  alpha * values1234.at(6); // E0110
			E[7] = alpha * values1234.at(7); // E0111

		}
		if (values1234.size() == 4){
			E[0] =  alpha * values1234.at(0); // E0000
			E[1] =  alpha * values1234.at(1); // E0001
			E[2] =  alpha * values1234.at(2); // E0010
			E[3] =  alpha * values1234.at(3); // E0011

		}



	}


	template<typename real>
	real GeneratorPseudoBoolean<real>::minimize_version_2(vector<label>& x, int& nlabelled) const
	{
		//ASSERT_STR(this->gen.ngen4 == 0, "Degree-4 generators are not yet supported.");

		index nVars = index( x.size() ); // Number of variables.
		int n = 2 * nVars;  // Because we have x and y.
		int num_cliques = 0;

		for (auto itr = alphaij.begin(); itr != alphaij.end(); ++itr) {
			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen2;++ii) {
				real alpha = vec.at(ii);
				if (alpha > 0) {
					num_cliques++;
				}
			}
		}

		for (auto itr = alphaijk.begin(); itr != alphaijk.end(); ++itr) {
			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen3;++ii) {
				real alpha = vec.at(ii);
				if (alpha > 0) {
					// Add monomials for this generator to the graph
					num_cliques++;
				}
			}
		}

		for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
			const auto& vec = itr->second;
			for (int ii = 0; ii < std::max(gen.ngen4pos, gen.ngen4neg); ++ii) {
				real alpha = vec.at(ii);
				if (alpha > 0) {
					// Add monomials for this generator to the graph
					num_cliques++;
				}
			}
		}



		real C = 0; // Constant in objective function.
		int clique_size = 4;
		int num_cliques_per_node = 2 * num_cliques; // TODO: Fix this. (Is this parameter used by GC?)

		// We add two extra variables in order to be able to add degree-2 cliques
		// as degree-4 cliques.

		typedef PRGC GCType;
		//typedef APGC GCType;

		GCType graph(n + 2,
			2 * num_cliques, // Each generator gives two cliques.
			clique_size,
			num_cliques_per_node);
		int extra1 = n;
		int extra2 = n + 1;

		std::unique_ptr<PseudoBoolean<real>> f_debug;

		// Uncomment this line to also minimize g exhaustively.
		//f_debug.reset(new PseudoBoolean<real>);

		//
		// Degree-1 terms.
		//
		for (auto itr = alphai.begin(); itr != alphai.end(); ++itr) {
			int i = itr->first;
			real alpha = itr->second;

			graph.AddUnaryTerm(i,         0,     alpha);
			graph.AddUnaryTerm(i + nVars, alpha,     0);

			if (f_debug) {
				f_debug->add_clique(i, 0, alpha);
				f_debug->add_clique(i + nVars, alpha, 0);
			}
		}

		//
		// Go through all alphas which correspond to quadratic generators
		//
		for (auto itr = alphaij.begin(); itr != alphaij.end(); ++itr) {
			const pair& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);

			vector<int> idx(4); // Translates from "local" indices to "global"
			idx[0] = i; // x variables
			idx[1] = j;
			idx[2] = i + nVars; // y variables
			idx[3] = j + nVars;

			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen2;++ii) {
				real alpha = vec.at(ii);
				auto& generator = gen.gen2.at(ii);
				if (alpha > 0) {
					// Add cliques for this generator to the graph
					{
						float E1[]= {alpha * generator.values1.at(0), // E0000
							alpha * generator.values1.at(1), // E0001
							alpha * generator.values1.at(2), // E0010
							alpha * generator.values1.at(3), // E0011
							alpha * generator.values1.at(0), // E0100
							alpha * generator.values1.at(1), // E0101
							alpha * generator.values1.at(2), // E0110
							alpha * generator.values1.at(3), // E0111
							alpha * generator.values1.at(0), // E1000
							alpha * generator.values1.at(1), // E1001
							alpha * generator.values1.at(2), // E1010
							alpha * generator.values1.at(3), // E1011
							alpha * generator.values1.at(0), // E1100
							alpha * generator.values1.at(1), // E1101
							alpha * generator.values1.at(2), // E1110
							alpha * generator.values1.at(3)};// E1111
						int indices1[] = {extra1,
							extra2,
							idx.at(generator.indices1.at(0)),
							idx.at(generator.indices1.at(1))};
						C += make_clique_positive(clique_size, E1);
						graph.AddHigherTerm(indices1, E1);
						cout <<"minimize_2 2g_1:         " << "ii2: " << generator.indices1.at(0) <<" " << "jj2: " << generator.indices1.at(1) << endl;


						if (f_debug) {
							std::vector<real> Ev(E1, E1+16);
							f_debug->add_clique(indices1[0], indices1[1], indices1[2], indices1[3], Ev);
						}
					}

					{
						float E2[]= {alpha * generator.values2.at(0), // E0000
							alpha * generator.values2.at(1), // E0001
							alpha * generator.values2.at(2), // E0010
							alpha * generator.values2.at(3), // E0011
							alpha * generator.values2.at(0), // E0100
							alpha * generator.values2.at(1), // E0101
							alpha * generator.values2.at(2), // E0110
							alpha * generator.values2.at(3), // E0111
							alpha * generator.values2.at(0), // E1000
							alpha * generator.values2.at(1), // E1001
							alpha * generator.values2.at(2), // E1010
							alpha * generator.values2.at(3), // E1011
							alpha * generator.values2.at(0), // E1100
							alpha * generator.values2.at(1), // E1101
							alpha * generator.values2.at(2), // E1110
							alpha * generator.values2.at(3)};// E1111
						int indices2[] = {extra1,
							extra2,
							idx.at(generator.indices2.at(0)),
							idx.at(generator.indices2.at(1))};
						C += make_clique_positive(clique_size, E2);
						graph.AddHigherTerm(indices2, E2);

						cout <<"minimize_2 2g_2:         " << "ii2: " << generator.indices2.at(0) <<" " << "jj2: " << generator.indices2.at(1) << endl;

						if (f_debug) {
							std::vector<real> Ev(E2, E2+16);
							f_debug->add_clique(indices2[0], indices2[1], indices2[2], indices2[3], Ev);
						}
					}
				}
			}
		}

		//
		// Go through all alphas which correspond to cubic generators
		//
		for (auto itr = alphaijk.begin(); itr != alphaijk.end(); ++itr) {


			const triple& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);

			//	cout<<"#################################:       " << "i: " << i << " j: " << j << " k: " << k << endl;


			vector<int> idx(6); // Translates from "local" indices to "global"
			idx.at(0) = i; // x variables
			idx.at(1) = j;
			idx.at(2) = k;
			idx.at(3) = i + nVars; // y variables
			idx.at(4) = j + nVars;
			idx.at(5) = k + nVars;

			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen3;++ii) {
				real alpha = vec.at(ii);
				auto& generator = gen.gen3.at(ii);
				if (alpha > 0) {
					// Add cliques for this generator to the graph

					{
						float E1[]= {alpha * generator.values1.at(0), // E0000
							alpha * generator.values1.at(1), // E0001
							alpha * generator.values1.at(2), // E0010
							alpha * generator.values1.at(3), // E0011
							alpha * generator.values1.at(4), // E0100
							alpha * generator.values1.at(5), // E0101
							alpha * generator.values1.at(6), // E0110
							alpha * generator.values1.at(7), // E0111
							alpha * generator.values1.at(0), // E1000
							alpha * generator.values1.at(1), // E1001
							alpha * generator.values1.at(2), // E1010
							alpha * generator.values1.at(3), // E1011
							alpha * generator.values1.at(4), // E1100
							alpha * generator.values1.at(5), // E1101
							alpha * generator.values1.at(6), // E1110
							alpha * generator.values1.at(7)};// E1111
						int indices1[] = {extra1,
							idx.at(generator.indices1.at(0)),
							idx.at(generator.indices1.at(1)),
							idx.at(generator.indices1.at(2))};
						C += make_clique_positive(clique_size, E1);
						graph.AddHigherTerm(indices1, E1);

						cout <<"minimize_2 3g_1:         " << "extra: " << extra1 <<" ii: " << idx.at(generator.indices1.at(0)) <<" " << "jj: " << idx.at(generator.indices1.at(1)) <<" " << "kk: " << idx.at(generator.indices1.at(2)) << endl;
						cout << "E1" << endl;
						for(int i = 0; i< 16; i++){
							cout << "i: "  << i <<" : " << E1[i] << endl;
						}
						cout << endl;


						if (f_debug) {
							std::vector<real> Ev(E1, E1+16);
							f_debug->add_clique(indices1[0], indices1[1], indices1[2], indices1[3], Ev);
						}
					}

					{
						float E2[]= {alpha * generator.values2.at(0), // E0000
							alpha * generator.values2.at(1), // E0001
							alpha * generator.values2.at(2), // E0010
							alpha * generator.values2.at(3), // E0011
							alpha * generator.values2.at(4), // E0100
							alpha * generator.values2.at(5), // E0101
							alpha * generator.values2.at(6), // E0110
							alpha * generator.values2.at(7), // E0111
							alpha * generator.values2.at(0), // E1000
							alpha * generator.values2.at(1), // E1001
							alpha * generator.values2.at(2), // E1010
							alpha * generator.values2.at(3), // E1011
							alpha * generator.values2.at(4), // E1100
							alpha * generator.values2.at(5), // E1101
							alpha * generator.values2.at(6), // E1110
							alpha * generator.values2.at(7)};// E1111
						int indices2[] = {extra1,
							idx.at(generator.indices2.at(0)),
							idx.at(generator.indices2.at(1)),
							idx.at(generator.indices2.at(2))};
						C += make_clique_positive(clique_size, E2);
						graph.AddHigherTerm(indices2, E2);

						cout <<"minimize_2 3g_2:         " << "extra " << extra1 <<   " ii: " << idx.at(generator.indices2.at(0)) <<" " << "jj: " << idx.at(generator.indices2.at(1)) <<" " << "kk: " << idx.at(generator.indices2.at(2)) << endl;
						cout << "E2" << endl;
						for(int i = 0; i< 16; i++){
							cout << "i: "  << i << " : " << E2[i] << endl;
						}
						cout << endl;

						if (f_debug) {
							std::vector<real> Ev(E2, E2+16);
							f_debug->add_clique(indices2[0], indices2[1], indices2[2], indices2[3], Ev);
						}
					}
				}
			}
		}

		//
		// Go through all alphas which correspond to quartic generators
		//

		for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
			const quad& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);
			int l=get_l(ind);
			const auto& vec = itr->second;

			vector<int> idx(8); // Translates from "local" indices to "global"
			idx.at(0) = i; // x variables
			idx.at(1) = j;
			idx.at(2) = k;
			idx.at(3) = l;
			idx.at(4) = i + nVars; // y variables
			idx.at(5) = j + nVars;
			idx.at(6) = k + nVars;
			idx.at(7) = l + nVars;



			// Was a positive or negative generator used?
			bool pos = false;
			auto mitr = posgen4.find(ind);
			if (mitr != posgen4.end()) {
				pos = mitr->second;
			}

			ASSERT(gen.ngen4pos == gen.ngen4neg);
			// The code below assumes a clique size of 4.
			ASSERT(clique_size == 4);

			for (int ii=0;ii<gen.ngen4pos;++ii) {
				real alpha = vec.at(ii);

				if (alpha > 0) {

					// Add monomials for this generator to the graph
					if (pos) {
						// Positive generator was used
						auto& generator = gen.gen4pos.at(ii);
						int ii, jj, kk, ll;

						ii = idx.at(generator.indices1.at(0));
						jj = idx.at(generator.indices1.at(1));
						kk = idx.at(generator.indices1.at(2));
						ll = idx.at(generator.indices1.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values1, alpha, f_debug);
						cout <<"minimize_2 4g_1:         " << "ii: " <<ii <<" " << "jj: " << jj <<" " << "kk: " << kk <<" " << "ll: " << ll << endl;


						ii = idx.at(generator.indices2.at(0));
						jj = idx.at(generator.indices2.at(1));
						kk = idx.at(generator.indices2.at(2));
						ll = idx.at(generator.indices2.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values2, alpha, f_debug);
						cout <<"minimize_2 4g_2:         " << "ii: " << ii <<" " << "jj: " << jj <<" " << "kk: " << kk <<" " << "ll: " << ll << endl;

					}
					else {
						// Negative generator was used
						auto& generator = gen.gen4neg.at(ii);
						int ii, jj, kk, ll;

						ii = idx.at(generator.indices1.at(0));
						jj = idx.at(generator.indices1.at(1));
						kk = idx.at(generator.indices1.at(2));
						ll = idx.at(generator.indices1.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values1, alpha, f_debug);

						ii = idx.at(generator.indices2.at(0));
						jj = idx.at(generator.indices2.at(1));
						kk = idx.at(generator.indices2.at(2));
						ll = idx.at(generator.indices2.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values2, alpha, f_debug);
					}
				}
			}
		}


		double min_g = constant + C + graph.FindMaxFlow();
		vector<label> xfull(n);
		for (int i = 0; i < n; ++i) {
			xfull[i] = graph.GetLabel(i);
		}

		if (f_debug) {
			std::cout << "Generic cuts\n";
			std::cout << "C=" << C << " min_g=" << min_g << "\n";

			for (int i = 0; i < nVars; ++i) {
				std::cout << xfull[i];
			}
			std::cout << ", ";
			for (int i = 0; i < nVars; ++i) {
				std::cout << xfull[i + nVars];
			}
			std::cout << "\n";
		}

		nlabelled = 0;
		for (int i=0; i<nVars; ++i) {
			bool used = false;
			auto itr = var_used.find(i);
			if (itr != var_used.end()) {
				used = itr->second;
			}

			if (used) {
				x[i]     = xfull[i];
				label yi = xfull[i+nVars];
				if (x[i] == yi) {
					x[i] = -1;
				}
				else {
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


		if (f_debug) {
			//
			// Minimize f_debug with exhaustive search.
			//
			vector<label> x_debug(n + 2, 0), x_debug_opt(n + 2, 0);
			real optimum = f_debug->eval(x_debug);
			while (true) {
				x_debug[0]++;
				int i=0;
				while (x_debug[i]>1) {
					x_debug[i]=0;
					i++;
					if (i == n + 2) {
						break;
					}
					x_debug[i]+=1;
				}
				if (i == n + 2) {
					break;
				}

				real energy = f_debug->eval(x_debug);
				if (energy < optimum) {
					optimum = energy;
					x_debug_opt = x_debug;
				}
			}

			std::cout << "Exhaustive debug\n";
			std::cout << "C=" << C << " min_f_debug=" << constant + C + optimum << "\n";
			for (int i = 0; i < nVars; ++i) {
				std::cout << x_debug_opt[i];
			}
			std::cout << ", ";
			for (int i = 0; i < nVars; ++i) {
				std::cout << x_debug_opt[i + nVars];
			}
			std::cout << "\n";
		}

		return min_g;
	}



	//adds a triple to the existing clique E, check if its already insearted.
	//poss correspond the which variables out of four that is used.
	template<typename real>
	void GeneratorPseudoBoolean<real>::add_triplet(int poss, int ii, int jj, int kk, float* E,  map<triple, vector<float>>& alpha_ijk, int nVars) const
	{
		auto itr1 = alpha_ijk.find(make_triple(ii, jj, kk));
		cout << "ADD_TRIPLET: " << endl;
		cout << "ii:  " << ii <<" " << "jj    " << jj <<" " << "kk :  " << kk << endl;

		//dosent exist!
		if(itr1 != alpha_ijk.end()){
			vector<float> vec = itr1->second;
			float E1[8];  //meaningless cpy
			for(int i = 0; i < 8; ++i) E1[i] = vec[i];

			insert_clique2(poss,E,E1);
			cout <<RED <<"ADDED TRIPLE: "  <<  "ii:  " << ii <<" " << "jj:  " << jj <<" " << "kk:  " << kk << NORMAL <<  endl;
			itr1 = alpha_ijk.erase(itr1);
		}
	}

	template<typename real>
	void GeneratorPseudoBoolean<real>::add_pair(int poss,int ii, int jj, float* E,  map<pair, vector<float> >& alpha_ij, int nVars) const
		{
			auto itr1 = alpha_ij.find(make_pair(ii, jj));
			cout << "ADD_PAIR: " << endl;
			cout << "ii:  " << ii <<" " << "jj    " << jj << endl;

			if(itr1 != alpha_ij.end()){
				vector<float> vec = itr1->second;
				float E1[4];  //meaningless cpy
				for(int i = 0; i < 4; ++i) E1[i] = vec[i];

				insert_clique2(poss,E,E1);
				cout <<RED <<"ADDED TRIPLE: "  <<  "ii:  " << ii <<" " << "jj:  " << jj << NORMAL <<  endl;
				itr1 = alpha_ij.erase(itr1);
			}
		}


		//Version where smaller cliques are inserted into bigger ones.
		template<typename real>
		real GeneratorPseudoBoolean<real>::minimize_version_3(vector<label>& x, int& nlabelled) const
		{
			cout << "#########################"<< endl;
			cout <<"inside minimize_3" << endl;
			cout << "#########################"<< endl;
			//maps with elements if generator already added.
			map<triple, int> inserted3;   
			map<pair, int> inserted2;

			//ASSERT_STR(this->gen.ngen4 == 0, "Degree-4 generators are not yet supported.");

			index nVars = index( x.size() ); // Number of variables.
			int n = 2 * nVars;  // Because we have x and y.
			int num_cliques = 0;


			for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
				const auto& vec = itr->second;
				for (int ii = 0; ii < std::max(gen.ngen4pos, gen.ngen4neg); ++ii) {
					real alpha = vec.at(ii);
					if (alpha > 0) {
						// Add monomials for this generator to the graph
						num_cliques++;
					}
				}
			}

			real C = 0; // Constant in objective function.
			int clique_size = 4;
			int num_cliques_per_node = 2 * num_cliques; // TODO: Fix this. (Is this parameter used by GC?)

			// We add two extra variables in order to be able to add degree-2 cliques
			// as degree-4 cliques.

			typedef PRGC GCType;
			//typedef APGC GCType;

			GCType graph(n + 2,
				2 * num_cliques, // Each generator gives two cliques.
				clique_size,
				num_cliques_per_node);
			int extra1 = n;
			int extra2 = n + 1;

			std::unique_ptr<PseudoBoolean<real>> f_debug;

			// Uncomment this line to also minimize g exhaustively.
			//f_debug.reset(new PseudoBoolean<real>);

			//
			// Degree-1 terms.
			//
			for (auto itr = alphai.begin(); itr != alphai.end(); ++itr) {
				int i = itr->first;
				real alpha = itr->second;

				graph.AddUnaryTerm(i,         0,     alpha);
				graph.AddUnaryTerm(i + nVars, alpha,     0);

				if (f_debug) {
					f_debug->add_clique(i, 0, alpha);
					f_debug->add_clique(i + nVars, alpha, 0);
				}
			}




			// new part, create a new map for all alphaijk where we can hold any index.
			// in this map we only hold values (i,j,k) and a vector with clique energies.
			map<triple,vector<float>> alpha_ijk;
			for(auto itr = alphaijk.begin(); itr != alphaijk.end(); ++itr){
				const triple& ind = itr-> first;
				int i = get_i(ind);
				int j = get_j(ind);
				int k = get_k(ind);
				const auto& vec = itr->second;
				vector<int> idx(6);
				idx.at(0) = i;
				idx.at(1) = j;
				idx.at(2) = k;
				idx.at(3) = i + nVars;
				idx.at(4) = j + nVars;
				idx.at(5) = k + nVars;

				for (int a=0;a<gen.ngen3;++a) {
					real alpha = vec.at(a);
					if (alpha > 0) {
						auto& generator = gen.gen3.at(a);
						int ii, jj, kk;

						//first part of the symmetric generator.
						ii = idx.at(generator.indices1.at(0));
						jj = idx.at(generator.indices1.at(1));
						kk = idx.at(generator.indices1.at(2));

						vector<float> E1;
						for(int b = 0; b< 8; b++) E1.push_back(alpha * generator.values1.at(b));

						alpha_ijk.insert(make_pair(make_triple(ii,jj,kk), E1));

						//symmetric part!
						ii = idx.at(generator.indices2.at(0));
						jj = idx.at(generator.indices2.at(1));
						kk = idx.at(generator.indices2.at(2));
						vector<float> E2;
						for(int b = 0; b< 8; b++) E2.push_back(alpha * generator.values2.at(b));

						alpha_ijk.insert(make_pair(make_triple(ii,jj,kk), E2));
					}
				}
			}


			// new part, create a new map for all alphaij where we can hold any index.
			// in this map we only hold values (i,j) and a vector with clique energies.
			map<pair,vector<float>> alpha_ij;
			for(auto itr = alphaij.begin(); itr != alphaij.end(); ++itr){
				const pair& ind = itr-> first;
				int i = get_i(ind);
				int j = get_j(ind);

				const auto& vec = itr->second;
				vector<int> idx(4);
				idx.at(0) = i;
				idx.at(1) = j;
				idx.at(2) = i + nVars;
				idx.at(3) = j + nVars;

				for (int a=0;a<gen.ngen2;++a) {
					real alpha = vec.at(a);
					if (alpha > 0) {
						auto& generator = gen.gen2.at(a);
						int ii, jj;

						//first part of the symmetric generator.
						ii = idx.at(generator.indices1.at(0));
						jj = idx.at(generator.indices1.at(1));

						vector<float> E1;
						for(int b = 0; b< 4; b++) E1.push_back(alpha * generator.values1.at(b));

						alpha_ij.insert(make_pair(make_pair(ii,jj), E1));

						//symmetric part!
						ii = idx.at(generator.indices2.at(0));
						jj = idx.at(generator.indices2.at(1));

						vector<float> E2;
						for(int b = 0; b< 4; b++) E2.push_back(alpha * generator.values2.at(b));

						alpha_ij.insert(make_pair(make_pair(ii,jj), E2));
					}
				}
			}
			cout << "---Quartic begin!----" << endl;
			for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
				
				const quad& ind = itr->first;
				int i=get_i(ind);
				int j=get_j(ind);
				int k=get_k(ind);
				int l=get_l(ind);
				const auto& vec = itr->second;

				vector<int> idx(8); // Translates from "local" indices to "global"
				idx.at(0) = i; // x variables
				idx.at(1) = j;
				idx.at(2) = k;
				idx.at(3) = l;
				idx.at(4) = i + nVars; // y variables
				idx.at(5) = j + nVars;
				idx.at(6) = k + nVars;
				idx.at(7) = l + nVars;


				// Was a positive or negative generator used?
				bool pos = false;
				auto mitr = posgen4.find(ind);
				if (mitr != posgen4.end()) {
					pos = mitr->second;
				}
				for (int i=0;i<gen.ngen4pos;++i) {
					real alpha = vec.at(i);
					if (alpha > 0) {
						if (pos) 
						{
							// Positive generator was used
							auto& generator = gen.gen4pos.at(i);
							int ii, jj, kk, ll;

							//first part och the symmetric generator.
							ii = idx.at(generator.indices1.at(0));
							jj = idx.at(generator.indices1.at(1));
							kk = idx.at(generator.indices1.at(2));
							ll = idx.at(generator.indices1.at(3));
							cout << "(i,j,k,l): " << ii << "," << jj << "," << kk << "," << ll << "     POSITIVE"  <<endl;

							ii = idx.at(generator.indices2.at(0));
							jj = idx.at(generator.indices2.at(1));
							kk = idx.at(generator.indices2.at(2));
							ll = idx.at(generator.indices2.at(3));
							cout << "(i,j,k,l): " << ii << "," << jj << "," << kk << "," << ll << "     POSITIVE"  <<endl;
						}
						else 
						{
							
							// Positive generator was used
							auto& generator = gen.gen4neg.at(i);
							int ii, jj, kk, ll;

							//first part och the symmetric generator.
							ii = idx.at(generator.indices1.at(0));
							jj = idx.at(generator.indices1.at(1));
							kk = idx.at(generator.indices1.at(2));
							ll = idx.at(generator.indices1.at(3));
							cout << "(i,j,k,l): " << ii << "," << jj << "," << kk << "," << ll << "     NEGATIVE"  <<   endl;

							ii = idx.at(generator.indices2.at(0));
							jj = idx.at(generator.indices2.at(1));
							kk = idx.at(generator.indices2.at(2));
							ll = idx.at(generator.indices2.at(3));
							cout << "(i,j,k,l): " << ii << "," << jj << "," << kk << "," << ll << "     NEGATIVE"  <<endl;
						}
					}
				}
			}



			cout << "---triplets begin!----" << endl;
			for(auto itr = alpha_ijk.begin(); itr != alpha_ijk.end(); ++itr){

				triple t = itr->first;
				cout << "(i,j,k): " << get_i(t) << "," << get_j(t) << "," << get_k(t) << endl;		
			}

			cout << "---Pair begin!----" << endl;
			for(auto itr = alpha_ij.begin(); itr != alpha_ij.end(); ++itr){
				pair t = itr->first;
				cout << "(i,j): " << get_i(t) << "," << get_j(t) << endl;		
			}


			//
			// Go through all alphas which correspond to quartic generators
			//

			for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {

				cout << "alphaijkl.size= " <<alphaijkl.size()  << endl;

				float E[16];
				const quad& ind = itr->first;
				int i=get_i(ind);
				int j=get_j(ind);
				int k=get_k(ind);
				int l=get_l(ind);
				const auto& vec = itr->second;

				vector<int> idx(8); // Translates from "local" indices to "global"
				idx.at(0) = i; // x variables
				idx.at(1) = j;
				idx.at(2) = k;
				idx.at(3) = l;
				idx.at(4) = i + nVars; // y variables
				idx.at(5) = j + nVars;
				idx.at(6) = k + nVars;
				idx.at(7) = l + nVars;

				// Was a positive or negative generator used?
				bool pos = false;
				auto mitr = posgen4.find(ind);
				if (mitr != posgen4.end()) {
					pos = mitr->second;
				}

				ASSERT(gen.ngen4pos == gen.ngen4neg);
				// The code below assumes a clique size of 4.
				ASSERT(clique_size == 4);

				for (int i=0;i<gen.ngen4pos;++i) {
					//cout << "gen.ngen4pos: " << gen.ngen4pos << endl;	
					real alpha = vec.at(i);

					if (alpha > 0) {
						cout << "i: " << i << endl; 

						if (pos) 
						{
							cout << GREEN<< "POSITIVE GENERATOR:" <<NORMAL  << endl;

							// Positive generator was used
							auto& generator = gen.gen4pos.at(i);
							int ii, jj, kk, ll;

							//first part och the symmetric generator.
							ii = idx.at(generator.indices1.at(0));
							jj = idx.at(generator.indices1.at(1));
							kk = idx.at(generator.indices1.at(2));
							ll = idx.at(generator.indices1.at(3));

							//creates the fourth-order klick E
							add_generators_to_clique(alpha, E, generator.values1);

							cout << "ii: " << ii <<" " << "jj: " << jj <<" " << "kk: " << kk <<" " << "ll: " << ll << endl;
							cout << WHITE << "FIRST PART" << NORMAL << endl;


							/////////Adding all positive-first triples!////////
							add_triplet(123,ii, jj, kk, E, alpha_ijk, nVars);
							add_triplet(124,ii, jj, ll, E, alpha_ijk, nVars);
							add_triplet(134,ii, kk, ll, E, alpha_ijk, nVars);
							add_triplet(234,jj, kk, ll, E, alpha_ijk, nVars);

							/////////Adding all positive-first pairs!////////
							add_pair(12,ii,jj,E,alpha_ij,nVars);
							add_pair(13,ii,kk,E,alpha_ij,nVars);
							add_pair(14,ii,ll,E,alpha_ij,nVars);
							add_pair(23,jj,kk,E,alpha_ij,nVars);
							add_pair(24,jj,ll,E,alpha_ij,nVars);
							add_pair(34,kk,ll,E,alpha_ij,nVars);

							int indices1[] = {ii, jj, kk, ll};
							C += make_clique_positive(4, E);
							graph.AddHigherTerm(indices1, E);


							//SYMMETRIC PART 2
							cout << WHITE << "SECOND PART:" << NORMAL << endl;
							ii = idx.at(generator.indices2.at(0));
							jj = idx.at(generator.indices2.at(1));
							kk = idx.at(generator.indices2.at(2));
							ll = idx.at(generator.indices2.at(3));

							//creates the fourth-order klick E
							add_generators_to_clique(alpha, E, generator.values2);

							cout << "ii: " << ii <<" " << "jj: " << jj <<" " << "kk: " << kk <<" " << "ll: " << ll << endl;


							/////////Adding all positive-second triples!////////
							add_triplet(123,ii, jj, kk, E, alpha_ijk, nVars);
							add_triplet(124,ii, jj, ll, E, alpha_ijk, nVars);
							add_triplet(134,ii, kk, ll, E, alpha_ijk, nVars);
							add_triplet(234,jj, kk, ll, E, alpha_ijk, nVars);


							/////////Adding all positive-second pairs!////////
							add_pair(12,ii,jj,E,alpha_ij,nVars);
							add_pair(13,ii,kk,E,alpha_ij,nVars);
							add_pair(14,ii,ll,E,alpha_ij,nVars);
							add_pair(23,jj,kk,E,alpha_ij,nVars);
							add_pair(24,jj,ll,E,alpha_ij,nVars);
							add_pair(34,kk,ll,E,alpha_ij,nVars);


							int indices2[] = {ii, jj, kk, ll};
							C += make_clique_positive(4, E);
							graph.AddHigherTerm(indices2, E);


						}
						else 
						{
							cout << GREEN<< "NEGATIVE GENERATOR" <<NORMAL  << endl;
							cout << WHITE << "FIRST PART" << NORMAL << endl;
							// Negative generator was used
							auto& generator = gen.gen4neg.at(i);
							int ii, jj, kk, ll;

							//first part och the symmetric generator.
							ii = idx.at(generator.indices1.at(0));
							jj = idx.at(generator.indices1.at(1));
							kk = idx.at(generator.indices1.at(2));
							ll = idx.at(generator.indices1.at(3));

							//creates the fourth-order klick E
							add_generators_to_clique(alpha, E, generator.values1);

							cout << "ii: " << ii <<" " << "jj: " << jj <<" " << "kk: " << kk <<" " << "ll: " << ll << endl;



							/////////Adding all negative-first triples!////////
							add_triplet(123,ii, jj, kk, E, alpha_ijk, nVars);
							add_triplet(124,ii, jj, ll, E, alpha_ijk, nVars);
							add_triplet(134,ii, kk, ll, E, alpha_ijk, nVars);
							add_triplet(234,jj, kk, ll, E, alpha_ijk, nVars);

							/////////Adding all negative-first pairs!////////
							add_pair(12,ii,jj,E,alpha_ij,nVars);
							add_pair(13,ii,kk,E,alpha_ij,nVars);
							add_pair(14,ii,ll,E,alpha_ij,nVars);
							add_pair(23,jj,kk,E,alpha_ij,nVars);
							add_pair(24,jj,ll,E,alpha_ij,nVars);
							add_pair(34,kk,ll,E,alpha_ij,nVars);

							int indices[] = {ii, jj, kk, ll};
							C += make_clique_positive(4, E);
							graph.AddHigherTerm(indices, E);



							//SYMMETRIC PART 2
							cout << WHITE << "SECOND PART:" << NORMAL << endl;
							ii = idx.at(generator.indices2.at(0));
							jj = idx.at(generator.indices2.at(1));
							kk = idx.at(generator.indices2.at(2));
							ll = idx.at(generator.indices2.at(3));

							//creates the fourth-order klick E
							add_generators_to_clique(alpha, E, generator.values2);

							cout << "ii: " << ii <<" " << "jj: " << jj <<" " << "kk: " << kk <<" " << "ll: " << ll << endl;


							/////////Adding all negative-second triples!////////
							add_triplet(123,ii, jj, kk, E, alpha_ijk, nVars);
							add_triplet(124,ii, jj, ll, E, alpha_ijk, nVars);
							add_triplet(134,ii, kk, ll, E, alpha_ijk, nVars);
							add_triplet(234,jj, kk, ll, E, alpha_ijk, nVars);


							/////////Adding all negative-second pairs!////////
							add_pair(12,ii,jj,E,alpha_ij,nVars);
							add_pair(13,ii,kk,E,alpha_ij,nVars);
							add_pair(14,ii,ll,E,alpha_ij,nVars);
							add_pair(23,jj,kk,E,alpha_ij,nVars);
							add_pair(24,jj,ll,E,alpha_ij,nVars);
							add_pair(34,kk,ll,E,alpha_ij,nVars);


							int indices2[] = {ii, jj, kk, ll};
							C += make_clique_positive(4, E);
							graph.AddHigherTerm(indices2, E);


						}
					}
				}
			}

			cout << "---triplets left!----" << endl;
			for(auto itr = alpha_ijk.begin(); itr != alpha_ijk.end(); ++itr){

				triple t = itr->first;
				cout << "(i,j,k): " << get_i(t) << "," << get_j(t) << "," << get_k(t) << endl;		
			}
		

			// add the last triplets old style!
			for(auto itr = alpha_ijk.begin(); itr != alpha_ijk.end(); ++itr){
				vector<float> Ei = itr->second;
				float E1[] = {Ei.at(0), Ei.at(1), Ei.at(2), Ei.at(3), Ei.at(4),Ei.at(5), Ei.at(6), Ei.at(7),
					Ei.at(0), Ei.at(1),Ei.at(3), Ei.at(3), Ei.at(4),Ei.at(5), Ei.at(6), Ei.at(7)};
				C+=make_clique_positive(clique_size, E1);
				
				triple t1 = map_back(get_i(itr->first) , get_j(itr->first), get_k(itr->first), nVars);
				int indices1[] = {2*nVars, get_i(t1) , get_j(t1), get_k(t1)};
				cout << "extra1 " << 2*nVars << " ii: " << get_i(itr->first) << " jj: " << get_j(itr->first)  << " kk: " << get_k(itr->first) << endl; 

			//	for(int i = 0; i< 16; ++i) cout << "i: " << i << " : "<< E1[i] << endl;
				graph.AddHigherTerm(indices1, E1);
			}
			
		


			double min_g = constant + C + graph.FindMaxFlow();
			vector<label> xfull(n);
			for (int i = 0; i < n; ++i) {
				xfull[i] = graph.GetLabel(i);
			}

			if (f_debug) {
				std::cout << "Generic cuts\n";
				std::cout << "C=" << C << " min_g=" << min_g << "\n";

				for (int i = 0; i < nVars; ++i) {
					std::cout << xfull[i];
				}
				std::cout << ", ";
				for (int i = 0; i < nVars; ++i) {
					std::cout << xfull[i + nVars];
				}
				std::cout << "\n";
			}

			nlabelled = 0;
			for (int i=0; i<nVars; ++i) {
				bool used = false;
				auto itr = var_used.find(i);
				if (itr != var_used.end()) {
					used = itr->second;
				}

				if (used) {
					x[i]     = xfull[i];
					label yi = xfull[i+nVars];
					if (x[i] == yi) {
						x[i] = -1;
					}
					else {
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


			if (f_debug) {
				//
				// Minimize f_debug with exhaustive search.
				//
				vector<label> x_debug(n + 2, 0), x_debug_opt(n + 2, 0);
				real optimum = f_debug->eval(x_debug);
				while (true) {
					x_debug[0]++;
					int i=0;
					while (x_debug[i]>1) {
						x_debug[i]=0;
						i++;
						if (i == n + 2) {
							break;
						}
						x_debug[i]+=1;
					}
					if (i == n + 2) {
						break;
					}

					real energy = f_debug->eval(x_debug);
					if (energy < optimum) {
						optimum = energy;
						x_debug_opt = x_debug;
					}
				}

				std::cout << "Exhaustive debug\n";
				std::cout << "C=" << C << " min_f_debug=" << constant + C + optimum << "\n";
				for (int i = 0; i < nVars; ++i) {
					std::cout << x_debug_opt[i];
				}
				std::cout << ", ";
				for (int i = 0; i < nVars; ++i) {
					std::cout << x_debug_opt[i + nVars];
				}
				std::cout << "\n";
			}

			return min_g;
		}







		template<typename real>
		real GeneratorPseudoBoolean<real>::minimize_version_1(vector<label>& x, int& nlabelled) const
		{
			index nVars = index( x.size() ); // Number of variables
			index var = 2*nVars; // Current variable

			real C = 0; // Constant in energy function
			Graph<real,real,real> graph(var, 50*var); //TODO: calculate the number of nodes and edges
			graph.add_node(var);

			//
			// Linear generator
			//
			for (auto itr = alphai.begin(); itr != alphai.end(); ++itr) {
				int i = itr->first;
				real alpha = itr->second;

				add_monomial_1_to_graph(C,graph,i,alpha);
				add_monomial_1_to_graph(C,graph,i+nVars,-alpha);
				C+=alpha;
			}

			//
			// Go through all alphas which correspond to quadratic generators
			//
			for (auto itr = alphaij.begin(); itr != alphaij.end(); ++itr) {
				const pair& ind = itr->first;
				int i=get_i(ind);
				int j=get_j(ind);
				const auto& vec = itr->second;
				for (int ii=0;ii<gen.ngen2;++ii) {
					real alpha = vec.at(ii);
					if (alpha > 0) {
						// Add monomials for this generator to the graph
						add_generator_to_graph(nVars, C, graph, alpha, gen.gen2red.at(ii), i,j );
					}
				}
			}

			//
			// Go through all alphas which correspond to cubic generators
			//
			for (auto itr = alphaijk.begin(); itr != alphaijk.end(); ++itr) {
				const triple& ind = itr->first;
				int i=get_i(ind);
				int j=get_j(ind);
				int k=get_k(ind);
				const auto& vec = itr->second;
				for (int ii=0;ii<gen.ngen3;++ii) {
					real alpha = vec.at(ii);
					if (alpha > 0) {
						// Add monomials for this generator to the graph
						add_generator_to_graph(nVars, C, graph, alpha, gen.gen3red.at(ii), i,j,k );
					}
				}
			}

			//
			// Go through all alphas which correspond to quartic generators
			//
			for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
				const quad& ind = itr->first;
				int i=get_i(ind);
				int j=get_j(ind);
				int k=get_k(ind);
				int l=get_l(ind);
				const auto& vec = itr->second;

				// Was a positive or negative generator used?
				bool pos = false;
				auto mitr = posgen4.find(ind);
				if (mitr != posgen4.end()) {
					pos = mitr->second;
				}

				ASSERT(gen.ngen4pos == gen.ngen4neg);

				for (int ii=0;ii<gen.ngen4pos;++ii) {
					real alpha = vec.at(ii);
					if (alpha > 0) {
						// Add monomials for this generator to the graph
						if (pos) {
							// Positive generator was used
							add_generator_to_graph(nVars, C, graph, alpha, gen.gen4redpos.at(ii), i,j,k,l );
						}
						else {
							// Negative generator was used
							add_generator_to_graph(nVars, C, graph, alpha, gen.gen4redneg.at(ii), i,j,k,l );
						}
					}
				}
			}

			// Compute the maximum flow and labeling
			real ming = constant + C + graph.maxflow();

			// Save solution from graph
			//
			int total_num_vars = graph.get_node_num();
			vector<char> xfull(total_num_vars);
			for (int i=0;i<total_num_vars; ++i) {
				xfull[i] = graph.what_segment(i);
			}

			//  try to obtain (0,1) and (1,0) solutions if possible
			//  does not matter that much for random polynomials
			vector<std::pair<int,int> > pairs;
			for (int i=0;i<nVars;++i) {
				pairs.push_back( std::make_pair(i, i+nVars) );
			}
			resolve_different(graph,xfull,pairs);

			// Extract labeling
			nlabelled = 0;
			for (int i=0; i<nVars; ++i) {
				bool used = false;
				auto itr = var_used.find(i);
				if (itr != var_used.end()) {
					used = itr->second;
				}

				if (used) {
					//x[i] = graph.what_segment(i);
					//label yi = graph.what_segment(i+nVars);;
					x[i] = xfull[i];
					label   yi = xfull[i+nVars];
					if (x[i] == yi) {
						x[i] = -1;
					}
					else {
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

			return ming;
		}






	}

#include "pb_instances.inc"
