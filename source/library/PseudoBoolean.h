//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// class PseudoBoolean : represents a pseudo-boolean polynomial of degree <= 4
// 
// class SymmetricPseudoBoolean : represents a symmetric pseudo-boolean polynomial of degree <= 4
//                                with special features, such as keeping a linear index for all
//                                coefficients (for LP solving).
// 
// The template argument to the classes specifies the underlying type. The double type is prone to
// rounding errors, although there is code to prevent that. But floating-point number is required
// to use linear programming. If only heuristics is required, int can be used instead.
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
	class GeneratorPseudoBoolean;

	template<typename real>
	class PseudoBoolean
	{
	public:
		template<typename real2> friend class SymmetricPseudoBoolean;
		template<typename real2> friend class GeneratorPseudoBoolean;
		template<typename real2,int degree> friend class Posiform;
		void print_helper(std::ostream& out) const;

		PseudoBoolean(); //Creates the zero polynomial
		PseudoBoolean(std::string filename); //Reads f from a file

		void save_to_file(std::string filename); // Saves f to file

		//Reset the function to the zero functions
		void clear();
		// Returns the maximum index used + 1
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

		// Minimize using higher-order clique reduction (HOCR, Ishikawa 2011)
		real minimize_reduction(vector<label>& x) const;
		real minimize_reduction(vector<label>& x, int& nlabelled) const;

		// Minimize using Fix et al., ICCV 2011
		real minimize_reduction_fixetal(vector<label>& x) const;
		real minimize_reduction_fixetal(vector<label>& x, int& nlabelled) const;
		
		// Minimize using LP relaxation
		// (persistency does not hold in general)
		// mostly for testing purposes
		real minimize_lp(vector<label>& x,bool verbose=false) const;

		// Minimizing using a symmetric, submodular function g(x,y)
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
	std::ostream& operator<<(std::ostream& out, const PseudoBoolean<real>& pbf)
	{
		pbf.print_helper(out);
		return out;
	}

	template<typename real>
	class SymmetricPseudoBoolean
	{
	public:
		template<typename real2> friend class PseudoBoolean;
		void print_helper(std::ostream& out);

		SymmetricPseudoBoolean();
		void clear();

		real eval(const vector<label>& x, const vector<label>& y) const;

		// Create using LP
		template<typename orgreal>
		void create_lp(const PseudoBoolean<orgreal>& pbf);

		// Create using heuristics
		template<typename orgreal>
		void create_heuristic(PseudoBoolean<orgreal>& pbf);

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

	private:
		//These are not to be modified directly
		int nlpvars;
		map<int, int> indbi;
		map<pair, int> indbij,indcij;
		map<triple, int> indbijk,indcijk,inddijk,indeijk;
		map<quad  , int> indbijkl,indcijkl,inddijkl,indeijkl,indpijkl,indqijkl,indrijkl,indsijkl;
	};

	template<typename real>
	std::ostream& operator<<(std::ostream& out, SymmetricPseudoBoolean<real>& pbf)
	{
		pbf.print_helper(out);
		return out;
	}

	template<typename real>
	class GeneratorPseudoBoolean
	{
	public:
		template<typename real2> friend class PseudoBoolean;
		GeneratorPseudoBoolean(std::string filename);
		void clear();

		// Create using LP
		void create_lp(const PseudoBoolean<real>& pbf);

		///////////
		// Index //
		///////////

		template <typename Map>
		int getindex(Map& m, const pair& key) 
		{
			auto itr = m.find( key );
			if ( itr == m.end() ) {
				m[key] = nlpvars;
				int tmp = nlpvars;
				nlpvars+=ngen2;
				return tmp;
			}
			else {
				return itr->second;
			}
		}
		template <typename Map>
		int getindex(Map& m, const triple& key) 
		{
			auto itr = m.find( key );
			if ( itr == m.end() ) {
				m[key] = nlpvars;
				int tmp = nlpvars;
				nlpvars+=ngen3;
				return tmp;
			}
			else {
				return itr->second;
			}
		}
		template <typename Map>
		int getindex(Map& m, const quad& key) 
		{
			auto itr = m.find( key );
			if ( itr == m.end() ) {
				m[key] = nlpvars;
				int tmp = nlpvars;
				nlpvars+=ngen4;
				return tmp;
			}
			else {
				return itr->second;
			}
		}

		int iaa(int i, int j, int k, int l);
		int ibb(int i, int j, int k);
		int icc(int i, int j);

	private:
		int ngen2, ngen3, ngen4, ngen4pos, ngen4neg;
		size_t nentries2, nentries3, nentries4;

		std::vector<int> cc;
		std::vector<int> bb12, bb13, bb23, bb123;
		std::vector<int> aa12pos, aa13pos, aa14pos, aa23pos, aa24pos, aa34pos, aa123pos, aa124pos, aa134pos, aa234pos, aa1234pos;
		std::vector<int> aa12neg, aa13neg, aa14neg, aa23neg, aa24neg, aa34neg, aa123neg, aa124neg, aa134neg, aa234neg, aa1234neg;

		int nlpvars;
		map<pair, int> indcc;
		map<triple, int> indbb;
		map<quad, int> indaa;

		map<pair,   vector<real> > alphaij;
		map<triple, vector<real> > alphaijk;
		map<quad,   vector<real> > alphaijkl;
	};

}

#endif
