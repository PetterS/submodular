#include <typeinfo>

#include "PseudoBoolean.h"

#include <coin/ClpSimplex.hpp>

#include <fstream>

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

	template<typename real> 
	GeneratorPseudoBoolean<real>::GeneratorPseudoBoolean(string filename)
	{
		ifstream fin(filename);

		fin >> ngen2 >> ngen3 >> ngen4pos;
		ASSERT_STR(fin, "Could not read file");

		ngen4=2*ngen4pos; //Half of generators required to be positive!
		ngen4neg = ngen4pos;

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

//		cout << "nentries234: " << nentries2 << " " << nentries3 << " " << nentries4 << endl;

	}
	//alpha[make_triple(i,j,k)].resize(ngen3);

	template<typename real>
	void GeneratorPseudoBoolean<real>::clear()
	{
		nlpvars = 0;

		alphaij.clear();
		alphaijk.clear();
		alphaijkl.clear();

		indaa.clear();
		indbb.clear();
		indcc.clear();
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


		int nLPVars = int( ngen2*pbf.aij.size() + ngen3*pbf.aijk.size() + ngen4pos*pbf.aijkl.size() );

		int nConstraints = int( pbf.aij.size() + pbf.aijk.size() + pbf.aijkl.size() ); // equality constraints


		// Compute the number of entries in the 
		// constraint matrix
		size_t nEntries = nentries2*pbf.aij.size() + nentries3*pbf.aijk.size() + nentries4*pbf.aijkl.size();
			
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

//			cout << "row: " << row << " col: " << col << " value: " << value << endl;
//			cout << rows.size() << endl;
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

			for (int ii=0;ii<ngen2;++ii){
				int colind=ind+ii;
				add_element(con, colind, cc[ii]); // cc[ii] should always be non-zero, since degree 2 generator
				//Objective function
				cost[colind] = -obj2[ii];
			}
			change_rhs(con, aij);
			con++;

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

			for (int ii=0;ii<ngen3;++ii){
				int colind=ind+ii;
				add_element(con, colind, bb123[ii]); // bb123[ii] should always be non-zero, since degree 3 generator
				if (bb12[ii]!=0){
					add_element(con_12, colind, bb12[ii]);
				}
				if (bb13[ii]!=0){
					add_element(con_13, colind, bb13[ii]);
				}
				if (bb23[ii]!=0){
					add_element(con_23, colind, bb23[ii]);
				}
				//Objective function
				cost[colind] = -obj3[ii];
			}
			change_rhs(con, aijk);
			con++;
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
				for (int ii=0;ii<ngen4pos;++ii){
					int colind=ind+ii;
					add_element(con, colind, aa1234pos[ii]); // aa1234[ii] should always be non-zero, since degree 4 generator
					if (aa12pos[ii]!=0){
						add_element(con_12, colind, aa12pos[ii]);
					}
					if (aa13pos[ii]!=0){
						add_element(con_13, colind, aa13pos[ii]);
					}
					if (aa14pos[ii]!=0){
						add_element(con_14, colind, aa14pos[ii]);
					}
					if (aa23pos[ii]!=0){
						add_element(con_23, colind, aa23pos[ii]);
					}
					if (aa24pos[ii]!=0){
						add_element(con_24, colind, aa24pos[ii]);
					}
					if (aa34pos[ii]!=0){
						add_element(con_34, colind, aa34pos[ii]);
					}
					if (aa123pos[ii]!=0){
						add_element(con_123, colind, aa123pos[ii]);
					}
					if (aa124pos[ii]!=0){
						add_element(con_124, colind, aa124pos[ii]);
					}
					if (aa134pos[ii]!=0){
						add_element(con_134, colind, aa134pos[ii]);
					}
					if (aa234pos[ii]!=0){
						add_element(con_234, colind, aa234pos[ii]);
					}
					//Objective function
					cost[colind] = -obj4pos[ii];
				}
			}
			else{
				for (int ii=0;ii<ngen4neg;++ii){
					int colind=ind+ii;
					add_element(con, colind, aa1234neg[ii]); // aa1234[ii] should always be non-zero, since degree 4 generator
					if (aa12neg[ii]!=0){
						add_element(con_12, colind, aa12neg[ii]);
					}
					if (aa13neg[ii]!=0){
						add_element(con_13, colind, aa13neg[ii]);
					}
					if (aa14neg[ii]!=0){
						add_element(con_14, colind, aa14neg[ii]);
					}
					if (aa23neg[ii]!=0){
						add_element(con_23, colind, aa23neg[ii]);
					}
					if (aa24neg[ii]!=0){
						add_element(con_24, colind, aa24neg[ii]);
					}
					if (aa34neg[ii]!=0){
						add_element(con_34, colind, aa34neg[ii]);
					}
					if (aa123neg[ii]!=0){
						add_element(con_123, colind, aa123neg[ii]);
					}
					if (aa124neg[ii]!=0){
						add_element(con_124, colind, aa124neg[ii]);
					}
					if (aa134neg[ii]!=0){
						add_element(con_134, colind, aa134neg[ii]);
					}
					if (aa234neg[ii]!=0){
						add_element(con_234, colind, aa234neg[ii]);
					}
					//Objective function
					cost[colind] = -obj4neg[ii];
				}
			}
			change_rhs(con, aijkl);
			con++;
		}

		cout << endl;
		cout << rows.size() << endl;
		cout << nEntries << endl;


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
			cout << "ii lpvars cost: " << ii << " " << lpvars[ii] << " " << cost[ii] << endl;
			obj+= (lpvars[ii]*cost[ii]);
		}
		cout << "Objective: " << obj << endl;

/*
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


*/




	} //end of create_lp
}

#include "pb_instances.inc"
