//
// Petter Strandmark 2011
// petter@maths.lth.se
//
#ifndef VERTEX_PACKING_PETTER
#define VERTEX_PACKING_PETTER

#include <sstream>
#include <vector>
#include <string>

#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
                                    sout << "Error (line " << __LINE__ << " in " << __FILE__ << "): " << #cond; \
                                    throw std::runtime_error(sout.str()); }
#define ASSERT_STR(cond,msg) if (!(cond)) { std::stringstream sout; \
                                            sout << "Error (line " << __LINE__ << " in " << __FILE__ << "): " << (msg); \
                                            throw std::runtime_error(sout.str()); }

namespace Petter 
{

	template<typename real>
	class VertexPacking
	{
	public:
		VertexPacking(size_t n);

		size_t size() { return w.size(); }

		void set_weight(size_t i, real w);
		void add_weight(size_t i, real w);
		void add_edge(size_t i, size_t j);

		real solve(std::vector<signed char>& x);
		real solve_slower(std::vector<signed char>& x);
		real solve_lp(std::vector<signed char>& x);
		real solve_exhaustive(std::vector<signed char>& x);

		bool is_feasible(const std::vector<signed char>& x) const;

		void load_from_file(std::string file);
		void save_to_file(std::string file);

		void print();

	protected:
		std::vector<real> w;
		std::vector< std::pair<int,int> > E;
	};


	template<typename real>
	class BipartiteVertexPacking
	{
	public:
		BipartiteVertexPacking(size_t m, size_t n);

		void set_left_weight(size_t i, real w);
		void set_right_weight(size_t i, real w);
		void add_edge(size_t i, size_t j);

		real solve(std::vector<signed char>& xm, std::vector<signed char>& xn);
		real solve_exhaustive(std::vector<signed char>& xmsol, std::vector<signed char>& xnsol);

		real energy(const std::vector<signed char>& xmsol, const std::vector<signed char>& xnsol) const;
		bool is_feasible(const std::vector<signed char>& xm, const std::vector<signed char>& xn) const;

		void print();

	protected:
		std::vector<real> wn,wm;
		std::vector< std::pair<int,int> > E;
	};

}

#endif