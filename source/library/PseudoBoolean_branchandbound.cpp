//
// Petter Strandmark 2012
// petter@maths.lth.se
//
// Minimizes a Pseudo-Boolean function using branch and bound
//


#include "PseudoBoolean.h"
#include "Petter-Color.h"
#include "Petter-Timer.h"

#include <iomanip>
#include <queue>
#include <cmath>
#include <cstdio>
using namespace std;

// From main program
void print_x(const std::vector<Petter::label>& x);

namespace Petter
{

	//
	// Subproblem
	//
	// A subproblem consists of a set of fixed variables and 
	// a function to minimize over the remaining
	//
	template<typename real>
	struct Subproblem
	{
		vector<label> fixed;
		PseudoBoolean<real> f;
		real bound;

		Subproblem(const vector<label>& x_in, const PseudoBoolean<real>& f_in) :
			fixed(x_in),
			f(f_in)
		{
			bound = -1000000000;
		}

		bool operator<(const Subproblem<real> right) const
		{
			//int nleft=0;
			//for (int i=0;i<fixed.size();++i) {
			//	if (fixed[i] >= 0) {
			//		nleft++;
			//	}
			//}

			//int nright=0;
			//for (int i=0;i<right.fixed.size();++i) {
			//	if (right.fixed[i] >= 0) {
			//		nright++;
			//	}
			//}

			//return nleft < nright;

			return bound < right.bound;
		}
	};

	template<typename real>
	struct BB_Format_Strings
	{
		static const char* iteration_string;
		static const char* completion_string;
	};


	template<typename real>
	real branch_and_bound(const PseudoBoolean<real>& f, vector<label>& x, BBInfo* bbinfo)
	{
		using namespace std;

		const bool verbose = false;

		size_t n = x.size();

		priority_queue<Subproblem<real> > problems;
		
		// Reset x
		for (size_t i=0;i<x.size();++i) {
			x[i] = -1;
		}
		Subproblem<real> start_problem(x,f);
		problems.push(start_problem);

		real upper_bound = 1000000000;
		vector<label> x_best;
		real f_best = -1;

		double solver_time = 0;
		double total_time  = 0;
		int    total_iterations = 0;
		int    nprints = 0;

		start();
		while (!problems.empty()) {
			total_time += stop(); 
			
			//
			// INFO: Prints progress information
			// 
			if ( abs(total_time - nprints) > 1.0) { // Every second
				nprints++;

				int prefix_length = 6; //64
				int found_prefix = 0;

				real lower_bound = 1000000000;

				vector<label> prefix(prefix_length,0);
				while (true) {
					prefix[0]++;
					int i=0;
					while (prefix[i]>1) {
						prefix[i]=0;
						i++;
						if (i==prefix_length) {
							break;
						}
						prefix[i]+=1;
					}
					if (i==prefix_length) {
						break;
					}

					// Does a problem with this prefix exist in the problem queue?
					bool found = false;
					priority_queue<Subproblem<real> > problems_copy = problems;
					while (!problems_copy.empty()) {
						const Subproblem<real>& problem = problems_copy.top();

						lower_bound = min(lower_bound, problem.bound);

						bool all_correct = true;
						for (int p=0;p<prefix_length;++p) {
							if (prefix.at(p) != problem.fixed.at(p) && problem.fixed.at(p) >= 0) {
								all_correct = false;
								break;
							}
						}
						if (all_correct) {
							found =true;
							break;
						}

						problems_copy.pop();
					}

					if (found) {
						found_prefix++;
					}
				}

				printf(BB_Format_Strings<real>::iteration_string,
				       total_iterations,
				       100*(pow(2.0,prefix_length) - found_prefix)/pow(2.0,prefix_length),
				       lower_bound,
				       upper_bound,
				       problems.size());
			}

			start();



			Subproblem<real> problem = problems.top();
			problems.pop();
			

			//
			// BOUND: Solve problem
			//        Compute a lower bound and persistencies
			//
			if (verbose) {
				cout << "Solving subproblem ";
				print_x(problem.fixed);
			}

			int nlabeled;

			start2();
			if (!bbinfo) {
				problem.bound = problem.f.minimize( problem.fixed, nlabeled, HOCR);
			}
			else if (bbinfo->method == GRD_heur) {

				//bound = problem.f.minimize( problem.fixed, nlabeled, true);

				// Don't iterate --- to save time
				SymmetricPseudoBoolean<real> g;
				g.create_heuristic(problem.f);
				problem.bound = g.minimize(problem.fixed, nlabeled);
			}
			else {
				problem.bound = problem.f.minimize( problem.fixed, nlabeled, bbinfo->method);
			}
			solver_time += stop2();

			total_iterations++;
			problem.f.reduce(problem.fixed);

			if (verbose) {
				cout << " solution = ";
				print_x(problem.fixed);
				cout << " bound = ";
			}

			bool should_branch = false;

			if (nlabeled == int(n) && problem.bound < upper_bound) {
				// Found a new best solution
				upper_bound = problem.bound;
				x_best = problem.fixed;
				f_best = problem.bound;

				if (verbose) {
					cout << WHITE;
				}
			}
			else if (nlabeled < int(n) && problem.bound <= upper_bound)  {
				// Need to branch, because this subproblem
				// is not completely solved
				should_branch = true;

				if (verbose) {
					cout << RED;
				}
			}
			else {
				// This branch can be discarded, because the lower bound 
				// obtained was higher than the global upper bound
				// Alternatively, nlabeled == n and we are done

				if (verbose) {
					cout << GREEN;
				}
			}
			
			if (verbose) {
				cout << problem.bound << NORMAL;
				cout << " upper = " << upper_bound << endl;
			}

			//
			// BRANCH: Replace the subproblem by two new problems
			//         with one more variable fixed
			//
			if (should_branch) {
				
				Subproblem<real> problem2 = problem;

				// Go through all variables
				for (size_t i=0;i<n;++i) {
					// Look for a variable which isn't fixed for this problem
					if (problem.fixed[i] < 0) {

						// Branch into two problems for this variable
						problem.fixed[i] = 0;
						problem.f.reduce(problem.fixed);

						problem2.fixed[i] = 1;
						problem2.f.reduce(problem2.fixed);

						if (verbose) {
							cout << "  New problem: "; print_x(problem.fixed); cout << endl;
							cout << "  New problem: "; print_x(problem2.fixed); cout << endl;
						}

						// Add the two new problems to the queue
						problems.push(problem);
						problems.push(problem2);

						break;
					}
				}
			}
			
			
		}

		total_time += stop();


		printf(BB_Format_Strings<real>::completion_string,
				       total_iterations,
				       f_best,
				       problems.size());

		x = x_best;
		if (bbinfo) {
			bbinfo->total_time = total_time;
			bbinfo->solver_time = solver_time;
			bbinfo->iterations = total_iterations;
		}

		return f_best;
	}



	template<typename real>
	const char* BB_Format_Strings<real>::iteration_string = "  Iteration";
	template<typename real>
	const char* BB_Format_Strings<real>::completion_string = "  Completed";

	template<>
	const char* BB_Format_Strings<double>::iteration_string  = "  Iters=%-8d %2.0f%% of branches done. bounds=(%7.1f, %7.1f) queue size: %d\n";
	template<>
	const char* BB_Format_Strings<double>::completion_string = "  Iters=%-7d 100%% of branches done. optimum=%7.1f\n";

	template<>
	const char* BB_Format_Strings<int>::iteration_string = "  Iters=%-8d %2.0f%% of branches done. bounds=(%6d, %6d) queue size: %d\n";
	template<>
	const char* BB_Format_Strings<int>::completion_string = "  Iters=%-7d 100%% of branches done. optimum=%6d\n";
}


#define INSTANTIATE(c) \
	template c     Petter::branch_and_bound<c>(const Petter::PseudoBoolean<c>& f, vector<label>& x, Petter::BBInfo* bbinfo);

INSTANTIATE(double);
INSTANTIATE(int);
