//
//
//
// class PseudoBoolean : represents a pseudo-boolean polynomial of degree <= 4
// 
// class SymmetricPseudoBoolean : represents a symmetric pseudo-boolean polynomial of degree <= 4
//                                with special features, such as keeping a linear index for all
//                                coefficients (for LP solving).
// 
//
#ifndef PETTER_PSEUDOBOOLEAN_H
#define PETTER_PSEUDOBOOLEAN_H


#include <vector>
#include <map>
#include <stdexcept>
#include <iostream>
#include <tuple>

#include <sstream>

#include "graph.h" //maxflow-v3.01

#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
                                    sout << "Error (line " << __LINE__ << " in " << __FILE__ << "): " << #cond; \
                                    throw std::runtime_error(sout.str()); }
#define ASSERT_STR(cond,msg) if (!(cond)) { std::stringstream sout; \
                                            sout << "Error (line " << __LINE__ << " in " << __FILE__ << "): " << (msg); \
                                            throw std::runtime_error(sout.str()); }

namespace Petter
{
	using std::map;
	using std::vector;

	typedef signed short label;
	typedef int index; 

	typedef std::pair<int,int> pair;
	typedef std::pair<int, std::pair<int,int> > triple;
	typedef std::pair<std::pair<int,int>, std::pair<int,int> > quad;

	// Indexing maps with tuples doesn't work very well in VC++ (known bug)
	// otherwise, the following typedefs would have been nice:
	//typedef std::tuple<int,int,int> triple;
	//typedef std::tuple<int,int,int,int> quad;
	

	pair make_pair(int i, int j);
	triple make_triple(int i, int j, int k);
	quad make_quad(int i, int j, int k, int l);
	int get_i(const pair& p);
	int get_j(const pair& p);
	int get_i(const triple& t);
	int get_j(const triple& t);
	int get_k(const triple& t);
	int get_i(const quad& q);
	int get_j(const quad& q);
	int get_k(const quad& q);
	int get_l(const quad& q);

	template<typename real>
	class SymmetricPseudoBoolean;

	template<typename real>
	class PseudoBoolean
	{
	public:
		friend class SymmetricPseudoBoolean<real>;
		void print_helper(std::ostream& out) const;

		PseudoBoolean(); //Creates the zero polynomial
		PseudoBoolean(std::string filename); //Reads f from a file

		void save_to_file(std::string filename); // Saves f to file

		void clear();
		int nvars() const;

		// Add monomials
		void add_monomial(int i, real a);
		void add_monomial(int i, int j, real a);
		void add_monomial(int i, int j, int k, real a);
		void add_monomial(int i, int j, int k, int l, real a);

		// Add entire cliques
		void add_clique(int i, real E0, real E1);
		void add_clique(int i, int j, real E00, real E01, real E10, real E11);
		void add_clique(int i, int j, int k, const vector<real>& E);
		void add_clique(int i, int j, int k, real E000, real E001, real E010, real E011,
											 real E100, real E101, real E110, real E111);
		void add_clique(int i, int j, int k, int l, const vector<real>& E);
		void add_clique(int i, int j, int k, int l, real E0000, real E0001, real E0010, real E0011,
		                                            real E0100, real E0101, real E0110, real E0111,
		                                            real E1000, real E1001, real E1010, real E1011,
		                                            real E1100, real E1101, real E1110, real E1111);

		// Evaluate the function
		real eval(const vector<label>& x) const;

		// Reduce the function given a partial labeling
		void reduce(const vector<label>& x);

		// Minimize using higher-order clique reduction
		real minimize_reduction(vector<label>& x) const;
		real minimize_reduction(vector<label>& x, int& nlabelled) const;

		// Minimizing using a symmetric function g(x,y)
		// NOTE: will change (reduce) *this
		real minimize(vector<label>& x, int& labeled, bool heuristic = false);

	protected:
		/////////////////////////////
		// Polynomial coefficients //
		/////////////////////////////
		real constant;
		map<int, real>    ai;
		map<pair  , real> aij;
		map<triple, real> aijk;
		map<quad  , real> aijkl;
	};

	template<typename real>
	std::ostream& operator<<(std::ostream& out, const PseudoBoolean<real>& pbf);

	template<typename real>
	class SymmetricPseudoBoolean
	{
	public:
		friend class PseudoBoolean<real>;
		void print_helper(std::ostream& out);

		SymmetricPseudoBoolean();
		void clear();

		real eval(const vector<label>& x, const vector<label>& y) const;

		// Create using LP
		void create_lp(const PseudoBoolean<real>& pbf);
		// Create using heuristics
		void create_heuristic(PseudoBoolean<real>& pbf);

		// Minimize; requires submodularity
		real minimize(vector<label>& x) const;
		real minimize(vector<label>& x, int& nlabelled) const;

		void test();
		void test_clear();

	protected:

		/////////////////////////////
		// Polynomial coefficients //
		/////////////////////////////
		real constant;
		map<int, real> bi;
		map<pair, real> bij,cij;
		map<triple, real> bijk,cijk,dijk,eijk;
		map<quad  , real> bijkl,cijkl,dijkl,eijkl,pijkl,qijkl,rijkl,sijkl;

		///////////
		// Index //
		///////////

		template <typename Map, typename K>
		int getindex(Map& m, const K& key) 
		{
			auto itr = m.find( key );
			if ( itr == m.end() ) {
				m[key] = nlpvars;
				return nlpvars++;
			}
			else {
				return itr->second;
			}
		}

		//int ib(int i);
		int ib(int i, int j);
		int ic(int i, int j);
		int ib(int i, int j, int k);
		int ic(int i, int j, int k);
		int id(int i, int j, int k);
		int ie(int i, int j, int k);
		int ib(int i, int j, int k, int l);
		int ic(int i, int j, int k, int l);
		int id(int i, int j, int k, int l);
		int ie(int i, int j, int k, int l);
		int ip(int i, int j, int k, int l);
		int iq(int i, int j, int k, int l);
		int ir(int i, int j, int k, int l);
		int is(int i, int j, int k, int l);

		//////////////////////////////
		// Graph (for minimization) //
		//////////////////////////////
		Graph<real,real,real>* maxflow_graph;

	private:
		//These are not to be modified directly
		int nlpvars;
		map<int, int> indbi;
		map<pair, int> indbij,indcij;
		map<triple, int> indbijk,indcijk,inddijk,indeijk;
		map<quad  , int> indbijkl,indcijkl,inddijkl,indeijkl,indpijkl,indqijkl,indrijkl,indsijkl;
	};

	template<typename real>
	std::ostream& operator<<(std::ostream& out, SymmetricPseudoBoolean<real>& pbf);

}

#endif
