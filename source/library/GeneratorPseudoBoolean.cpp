
#include <fstream>
#include <typeinfo>

// BK maxflow
#include "graph.h"
// LP solver
#include <coin/ClpSimplex.hpp>

#include "PseudoBoolean.h"
#include "Minimizer.h"
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

	template<typename real> 
	Generators<real>::Generators(const string& filename)
	{
		ifstream fin(filename);

		fin >> ngen2 >> ngen3 >> ngen4pos;
		ASSERT_STR(fin, "Could not read file");

		ngen4=2*ngen4pos; //Half of generators required to be positive!
		ngen4neg = ngen4pos;


		//
		// Read generated forms
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

		//
		// Read aa,bb and cc
		// 

		nentries2 = read_vec(cc,fin,ngen2);
		read_vec(obj2,fin,ngen2);

		ASSERT_STR(fin, "Could not read cc");

		nentries3  = read_vec(bb12,fin,ngen3);
		nentries3 += read_vec(bb13,fin,ngen3);
		nentries3 += read_vec(bb23,fin,ngen3);
		nentries3 += read_vec(bb123,fin,ngen3);

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
