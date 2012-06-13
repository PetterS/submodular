/// \mainpage 
///
/// Petter Strandmark 2012
/// petter@maths.lth.se
///
/// \section classes Main Classes
///  Petter::PseudoBoolean : represents a pseudo-boolean polynomial of degree <= 4
/// 
///  Petter::SymmetricPseudoBoolean : represents a symmetric pseudo-boolean polynomial of degree <= 4
///                               with special features, such as keeping a linear index for all
///                                coefficients (for LP solving).
/// 
/// The template argument to the classes specifies the underlying type. The double type is prone to
/// rounding errors, although there is code to prevent that. But floating-point number is required
/// to use linear programming. If only heuristics is required, int can be used instead.
///
///
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

/// The namespace everything resides in
namespace Petter
{
	using std::map;
	using std::vector;

	/// Used to specify optimization method
	enum Method {
		/// Reductions described by Boros and Hammer (2002). For comparisons only---HOCR is better.
		M_reduction,
		/// HOCR (Ishikawa, 2011)
		HOCR, 
		/// Fix et al. (2001)
		Fix,  
		/// Generalized roof duality
		GRD, 
		/// Heuristic GRD
		GRD_heur, 
		/// GRD with generators
		GRD_gen,  
		/// Linear program relaxation
		LP 
	};
	const char* str(Method m);

	/// Represents a boolean label (0 or 1).  Negative values mean undetermined.
	typedef signed short label;
	/// Represents an index 
	typedef int index; 

	/// A pair of indices
	typedef std::pair<int,int> pair;
	/// A 3-tuple of indices
	typedef std::pair<int, std::pair<int,int> > triple;
	/// A 4-tuple of indices
	typedef std::pair<std::pair<int,int>, std::pair<int,int> > quad;

	// Indexing maps with tuples doesn't work very well in VC++ (known bug)
	// otherwise, the following typedefs would have been nice:
	//typedef std::tuple<int,int,int> triple;
	//typedef std::tuple<int,int,int,int> quad;

	/// Creates a pair from two indices
	pair make_pair(int i, int j);
	/// Creates a triple from three indices
	triple make_triple(int i, int j, int k);
	/// Creates a quadruple from four indices
	quad make_quad(int i, int j, int k, int l);
	/// Retrieves the first index
	int get_i(const pair& p);
	/// Retrieves the second index
	int get_j(const pair& p);
	/// Retrieves the first index
	int get_i(const triple& t);
	/// Retrieves the second index
	int get_j(const triple& t);
	/// Retrieves the third index
	int get_k(const triple& t);
	/// Retrieves the first index
	int get_i(const quad& q);
	/// Retrieves the second index
	int get_j(const quad& q);
	/// Retrieves the third index
	int get_k(const quad& q);
	/// Retrieves the fourth index
	int get_l(const quad& q);

	template<typename real>
	class SymmetricPseudoBoolean;

	template<typename real>
	class GeneratorPseudoBoolean;

	//
	/// Represents a pseudo-boolean function of degree <= 4
	//
	/// Use this class to create and hold a pseudo-boolean function.
    /// The function can be saved to disk and minimized using various methods.
	template<typename real>
	class PseudoBoolean
	{
	public:
		template<typename real2> friend class SymmetricPseudoBoolean;
		template<typename real2> friend class GeneratorPseudoBoolean;
		template<typename real2,int degree> friend class Posiform;
        /// \private
		void print_helper(std::ostream& out) const;

		/// Creates the zero function
		PseudoBoolean(); 
		/// Reads a polynomial from file
		//
		/// \param filename the file name from which the contructor should read
		PseudoBoolean(std::string filename); 

		/// Saves the function to file
		//
		/// \param filename the file name the function should write to
		void save_to_file(std::string filename); 

		/// Resets the function to the zero functions
		void clear();
		/// Returns the maximum index used + 1
		int nvars() const;

		/// Adds a monomial of degree 1 to the function
		void add_monomial(int i, real a);
        /// Adds a monomial of degree 2 to the function
		void add_monomial(int i, int j, real a);
        /// Adds a monomial of degree 3 to the function
		void add_monomial(int i, int j, int k, real a);
        /// Adds a monomial of degree 4 to the function
		void add_monomial(int i, int j, int k, int l, real a);

		/// Adds a clique of a single variable to the function
		void add_clique(int i, real E0, real E1);
        /// Adds a clique of two variables to the function
		void add_clique(int i, int j, real E00, real E01, real E10, real E11);
        /// Adds a clique of three variables to the function
		void add_clique(int i, int j, int k, const vector<real>& E);
        /// Adds a clique of three variables to the function
		void add_clique(int i, int j, int k, real E000, real E001, real E010, real E011,
											 real E100, real E101, real E110, real E111);
        /// Adds a clique of four variables to the function
		void add_clique(int i, int j, int k, int l, const vector<real>& E);
        /// Adds a clique of four variables to the function
		void add_clique(int i, int j, int k, int l, real E0000, real E0001, real E0010, real E0011,
		                                            real E0100, real E0101, real E0110, real E0111,
		                                            real E1000, real E1001, real E1010, real E1011,
		                                            real E1100, real E1101, real E1110, real E1111);

		/// Evaluates the function
		real eval(const vector<label>& x) const;

		/// Reduce the function given a partial labeling
		void reduce(const vector<label>& x);

		/// \private
		real minimize_reduction(vector<label>& x) const;
        /// \private
		real minimize_reduction(vector<label>& x, int& nlabelled) const;

		/// \private
		real minimize_reduction_fixetal(vector<label>& x) const;
        /// \private
		real minimize_reduction_fixetal(vector<label>& x, int& nlabelled) const;
		
		/// \private
		real minimize_reduction_M(vector<label>& x, int& nlabelled) const;

		/// Minimize using LP relaxation

		/// Persistency does not hold in general. This function is
		/// mostly for testing purposes
		real minimize_lp(vector<label>& x,bool verbose=false) const;

		/// Minimizing using any method

		/// Same functionality as below.
		real minimize(vector<label>& x, Method method, const char* generators_file=0);
        /// Minimizing using any method

		/// NOTE: might change (reduce) *this
		/// \param x vector with solution. Has to have x.size() == n.
		/// \param labeled The number of variables the function were able to assign.
		/// \param method The minimization method to use
		/// \param generators_file The file from which to read the generators. (optional)
		real minimize(vector<label>& x, int& labeled, Method method = GRD, const char* generators_file=0);

		/// \private
		real minimize(vector<label>& x, int& labeled, bool heuristic = false);
        /// \private
		real minimize_generators(vector<label>& x, int& labeled, bool heuristic = false);

	protected:
		// //////////////////////////
		// Polynomial coefficients //
		// //////////////////////////
		/// constant term in function
		real constant; 
		/// Degree-1 terms
		map<int, real>    ai; 
		/// Degree-2 terms
		map<pair  , real> aij; 
		/// Degree-3 terms
		map<triple, real> aijk; 
		/// Degree-4 terms
		map<quad  , real> aijkl; 
	};

	template<typename real>
	std::ostream& operator<<(std::ostream& out, const PseudoBoolean<real>& pbf)
	{
		pbf.print_helper(out);
		return out;
	}

    /// Holds info about a branch-and-bound run. See branch_and_bound().
	struct BBInfo
	{
		/// (input) What minimization method should be used
		Method method; 
		 /// (output) How many iterations were used
		int iterations;
		/// (output) Total solver time
		double total_time; 
		/// (output) Time spent in pseudo-boolean solver
        double solver_time; 
	};

	/// Minimizes a pseudo-boolean function exactly using branch and bound
	///
	/// \param f the PseudoBoolean function to minimize
	/// \param x vector which holds the solution. Should satisfy x.size()==n
	/// \param bbinfo (optional) a BBInfo struct with extra input and output
	template<typename real>
	real branch_and_bound(const PseudoBoolean<real>& f, vector<label>& x, BBInfo* bbinfo=0);
    
	//
	/// Represents a *symmetric* pseudo-boolean function of degree <= 4
	//
	/// This holds a symmetric pseudo-boolean function; it is typically only used when
    /// minimizing another pseudo-boolean function.
    /// It can only be created from an existing pseudo-boolean function.
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

		void make_submodular(const PseudoBoolean<real>& pbf);

	private:

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


	//
	/// Helper class that holds informations about generators and how to reduce them
	//
	template<typename real>
	class Generators
	{
	public:

		Generators(const std::string& filename);

		// This internal struct holds a monomial of degrees 0 (i,j<0), 1 (j<0) or 2.
		struct Monomial
		{
			Monomial(real coef=0) 
			{
				i = j = -1;
				c = coef;
			}
			Monomial(int i_in, real coef) 
			{
				i = i_in;
				j = -1;
				c = coef;
			}
			Monomial(int i_in, int j_in, real coef) 
			{
				i = i_in;
				j = j_in;
				c = coef;
			}
			void check() const
			{
				ASSERT(i>=0 || j<0);
			}
			int i,j;
			real c;
		};

		int ngen2, ngen3, ngen4, ngen4pos, ngen4neg;
		size_t nentries2, nentries3, nentries4;

		// These variables hold the information needed to build the constraint matrix
		std::vector<int> cc, obj2;
		std::vector<int> bb12, bb13, bb23, bb123, obj3;
		std::vector<int> aa12pos, aa13pos, aa14pos, aa23pos, aa24pos, aa34pos, aa123pos, aa124pos, aa134pos, aa234pos, aa1234pos, obj4pos;
		std::vector<int> aa12neg, aa13neg, aa14neg, aa23neg, aa24neg, aa34neg, aa123neg, aa124neg, aa134neg, aa234neg, aa1234neg, obj4neg;

		// These variables hold the reduced forms of the generators
		std::vector< std::vector< Monomial > > gen2red, gen3red, gen4redpos, gen4redneg;
	};

	//
	/// Represents a *symmetric* submodular pseudo-boolean function of degree <= 4
	//
	/// This holds a symmetric pseudo-boolean function; it is typically only used when
    /// minimizing another pseudo-boolean function. The internal representation is based on
	/// the generators specified in the constructor.
    /// It can only be created from an existing pseudo-boolean function.
	template<typename real>
	class GeneratorPseudoBoolean
	{
	public:

		template<typename real2> friend class PseudoBoolean;
		/// Creates the zero function. 
		//
		/// The internal representation of the function is specified by the generators.
		GeneratorPseudoBoolean(const Generators<real>& generators);
		/// Resets the function to the zero function.
		void clear();

		/// Creates using linear programming
		void create_lp(const PseudoBoolean<real>& pbf);

		/// Minimizes the function. 
		real minimize(vector<label>& x, int& nlabelled) const;

	private:

		template <typename Map>
		int getindex(Map& m, const pair& key) 
		{
			auto itr = m.find( key );
			if ( itr == m.end() ) {
				m[key] = nlpvars;
				int tmp = nlpvars;
				nlpvars+=gen.ngen2;
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
				nlpvars+=gen.ngen3;
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
				nlpvars+=gen.ngen4pos;
				return tmp;
			}
			else {
				return itr->second;
			}
		}

		/// Gives the linear programming variable index for aa_{i,j,k}.
		/// Mostly used internally.
		int iaa(int i, int j, int k, int l);
		/// Gives the linear programming variable index for bb_{i,j,k}.
		/// Mostly used internally.
		int ibb(int i, int j, int k);
		/// Gives the linear programming variable index for bb_{i,j,k}.
		/// Mostly used internally.
		int icc(int i, int j);

		const Generators<real>& gen;

		// LP indices
		int nlpvars;
		map<pair, int> indcc;
		map<triple, int> indbb;
		map<quad, int> indaa;

		// Constant
		real constant;
		// Linear terms
		map<int, real> alphai;
		// Coefficients in front of generators
		map<pair,   vector<real> > alphaij;
		map<triple, vector<real> > alphaijk;
		map<quad,   vector<real> > alphaijkl;
		// Positive or negative generator used
		map<quad, bool> posgen4;
		// Keeps track of variables present
		map<int,bool> var_used;
	};

}

#endif
