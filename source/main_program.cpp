//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Demonstration program for generalized roof duality
//

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <new>
#include <cstdlib>
#include <ctime>
#include <map>
#include <list>
#include <random>
#include <functional>
#include <algorithm>
using namespace std;


#include "Petter-Color.h"
#include "PseudoBoolean.h"
#include "Posiform.h"
#include "Minimizer.h" // For tests
#include "Petter-Timer.h"
using namespace Petter;

// Functions running some quick tests
// In: submodular_tests.cpp
template<typename real> void test_pseudoboolean();
void test_minimize();
template<typename real> void test_posiform();

// Function creating a polynomial from command line options
// In: get_polynomial.cpp
void get_polynomial(std::map<std::string,std::string>& cmd_line,
                    Petter::PseudoBoolean<double>* pb);

// Random number generator
namespace {
	#ifdef WIN32
		// Time Stamp Counter gives better seeds than seconds
		// when many small problems are generated consecutively
		std::mt19937 engine((unsigned long)(__rdtsc()&0xffffffff));
	#else
		std::mt19937 engine(unsigned(time(0)));
	#endif
}

void print_x(const std::vector<label>& x)
{
	size_t n = x.size();
	for (size_t i=0;i<n&&i<20;++i) {
		if (x[i]>=0) {
			cout << x[i];
		}
		else {
			cout << '-';
		}
	}
	if (n>20) {
		cout << " ... ";
	}
}

template<typename real>
void print_info(std::string name,
                const std::vector<label>& x,
                real bound, int labeled,
                Petter::Color color, double time=-1)
{
	using namespace std;
	cout << left << setw(16) << name << "f(";
	print_x(x);
	cout << ") = " << color << right << setw(8) << bound << NORMAL;
	cout << ",   labeled : " << color << labeled << NORMAL;
	if (time>=0) {
		cout << ", time : " << color << time << NORMAL;
	}
	cout << endl;
}

//
// Checks that the persistencies in x agrees with some optimal solution
//
void check_persistency( const vector< vector<label> >& optimal_solutions,
                        const vector<label>& x, int reported_labelled=-1)
{
	bool any_ok = false;
	for (auto itr = optimal_solutions.begin(); itr != optimal_solutions.end(); ++itr) {
		bool this_ok = true;
		for (size_t i=0;i<x.size();++i) {
			if (x.at(i) >= 0 && x.at(i) != itr->at(i)) {
				this_ok = false;
			}
		}
		if (this_ok) {
			any_ok = true;
			break;
		}
	}
	ofstream log("percistency.data", ios::app);
	if (!any_ok) {
		// Persistency did not hold for this solution

		
		log <<  0      << '\t'  << endl;// 0




		throw runtime_error("Persistency error");
	}
	else{
		log <<  1      << '\t'  << endl;// 0
	}
	if (reported_labelled >= 0) {
		int labelled = 0;
		for (size_t i=0;i<x.size();++i) {
			if (x.at(i) >= 0) {
				labelled++;
			}
		}

		if (labelled != reported_labelled) {
			throw runtime_error("Incorrect number of labels reported");
		}
	}
}

//
// Checks that the lower bound agrees with the global optimum
//
template<typename real>
void check_bound(real optimum, real lower_bound)
{
	if (lower_bound > optimum) {
		throw runtime_error("Lower bound error");
	}
}
template<>
void check_bound<double>(double optimum, double lower_bound)
{
	// We have to have a loose check, because doubles have
	// rounding errors
	if ( (optimum - lower_bound) / abs(optimum) < -1e-7 ) {
		throw runtime_error("Lower bound error (double)");
	}
}

//
// Checks for equality
//
template<typename T>
bool isequal(T t1, T t2)
{
	return t1 == t2;
}
template<>
bool isequal<double>(double t1, double t2)
{
	if (t1==0 && t2==0) {
		return true;
	}
	return abs(t1-t2)/(abs(t1)+abs(t2)) < 1e-8;
}



template<typename real>
void test_branchandbound(const PseudoBoolean<real>& pb,
                         int n, Method method, ostream& fout,
                         Petter::Color color)
{
	vector<label> x(n,0);
	real val;
	BBInfo bbinfo;

	bbinfo.method = method;
	val = branch_and_bound(pb,x,&bbinfo);

	cout << setw(10) << left << str(method) << color << bbinfo.iterations << NORMAL
	     << " iterations, total_time = " << color << bbinfo.total_time << NORMAL
	     << " s, solver_time = " << color << bbinfo.solver_time << NORMAL << endl;
	fout << setw(10) << left << str(method) <<  bbinfo.iterations
	     << " iterations, total_time = " << bbinfo.total_time
	     << " s, solver_time = " << bbinfo.solver_time << endl;

  if (!isequal(val, pb.eval(x))) {
		throw runtime_error("Branch and bound value inconsistent");
	}
}


//
// Main program
//
int main_program(int num_args, char** args)
{
	using namespace std;
	using namespace Petter;

	// Need doubles for linear programming
	typedef double real;

	if (num_args == 1) {
		// Run some tests
		statusTry("Testing minimization...");
		test_minimize();
		statusTry("Testing pseudo-Boolean functions...");
		test_pseudoboolean<double>();
		//test_pseudoboolean<int>(); // Does not have LP and therefore fails
		statusTry("Testing posiform (double)...");
		//test_posiform<double>();
		statusTry("Testing posiform (int)...");
		//test_posiform<int>();
		statusOK();
		statusTry("Testing graph functions...");
		test_graph_functions<int>();
		statusOK();

		statusTry("Testing generators...");
		Generators<real> generators("generators/generators.txt");
		Petter::GeneratorPseudoBoolean<real> genpb(generators);
		statusOK();
		cerr << "Number of generators : (" << generators.ngen2 << ", " << generators.ngen3 << ", " << generators.ngen4 << ")" << endl;
		statusTry("Testing create lp...");
		PseudoBoolean<double> f("../tests/quartic_paper.txt");
		genpb.create_lp(f);
		statusOK();

		statusTry("Minimizing g...");
		vector<label> x(f.nvars(), -1);
		int nlabelled=-1;
		genpb.minimize(x, nlabelled);
		statusOK();

		cerr << "Possible choices : " << endl;
		cerr << "  " << args[0] << " -m <int> -n <int> -nterms <int>  : runs random examples" << endl;
		cerr << "  " << args[0] << " -m <int> -example                : examples from paper" << endl;
		cerr << "  " << args[0] << " -file <str>                      : read polynomial from file" << endl;
		cerr << "  " << args[0] << " -sat <str>                       : read SAT problem from file" << endl;
		cerr << endl;
		cerr << "    -optimal                         : use linear programming" << endl;
		cerr << "    -generators                      : use generators (with linear programming)" << endl;
		cerr << "        -generators-file             : file with generators to use (optional)" << endl;
		cerr << "    -heuristic                       : use heuristics" << endl;
		cerr << "    -fixetal                         : use reductions from Fix et al." << endl;
		cerr << "    -exhaustive                      : use exhaustive search (n<=30)" << endl;
		cerr << endl;
		cerr << "    -iterate                         : also iterate reduction methods" << endl;
		cerr << endl;
		cerr << "    -packing                         : compute solution using vertex packing" << endl;
		cerr << "    -lprelax                         : compute LP relaxation" << endl;
		cerr << endl;
		cerr << "    -verbose                         : print polynomials" << endl;
		cerr << endl;

		return 0;
	}

	//Command line
	map<string,string> cmd_line;
	//Default parameters
	cmd_line["-m"] = "3";
	cmd_line["-n"] = "0";
	cmd_line["-nterms"] = "0";
	cmd_line["-generators-file"] = "generators/generators.txt";

	//Read command line into map
	for (int i=1;i<num_args;++i) {
		string cmd = args[i];
		if (i<num_args-1 && args[i+1][0]!='-') {
			cmd_line[cmd] = args[i+1];
			++i;
		}
		else {
			cmd_line[cmd] = "";
		}
	}

	bool do_m = false;
	bool do_hocr = true;
	bool do_optimal = false;
	bool do_generators = false;
	bool do_fixetal = false;
	bool do_exhaustive = false;
	bool do_heuristic = false;
	bool do_lprelax = false;
	bool do_packing = false;
	bool iterate_reduction_methods = false;
	if (cmd_line.find("-optimal") != cmd_line.end() || cmd_line.find("-grd") != cmd_line.end()) {
		do_optimal = true;
	}
	if (cmd_line.find("-generators") != cmd_line.end() || cmd_line.find("-grd-gen") != cmd_line.end()) {
		do_generators = true;
	}
	if (cmd_line.find("-fixetal") != cmd_line.end()) {
		do_fixetal = true;
	}
	if (cmd_line.find("-m-reduction") != cmd_line.end() ) {
		do_m = true;
	}
    if (cmd_line.find("-hocr") != cmd_line.end()) {
		do_hocr = true;
	}
	if (cmd_line.find("-heuristic") != cmd_line.end() || cmd_line.find("-grd-heuristic") != cmd_line.end()) {
		do_heuristic = true;
	}
	if (cmd_line.find("-lprelax") != cmd_line.end()) {
		do_lprelax = true;
	}
	if (cmd_line.find("-exhaustive") != cmd_line.end()) {
		do_exhaustive = true;
	}
	if (cmd_line.find("-packing") != cmd_line.end()) {
		do_packing = true;
	}
	if (cmd_line.find("-iterate") != cmd_line.end()) {
		iterate_reduction_methods = true;
	}


	bool submodular = cmd_line.find("-submodular") != cmd_line.end();
	bool verbose = cmd_line.find("-verbose") != cmd_line.end();

	Petter::PseudoBoolean<real> pb;
	get_polynomial(cmd_line, &pb);
	int n = pb.nvars();

	cout << "Polynomial : " << pb << endl;
	cout << WHITE;
	cout << "n = " << n << endl;
	cout << NORMAL;

	// Save to temporary file
	if (std::getenv("TEMP")) {
		string tmp = std::getenv("TEMP");
		pb.save_to_file(tmp + "/pb.txt");
	}


	///////////////////////////////////
	// Solve using different methods //
	///////////////////////////////////

	int m_labeled = -1;
	int hocr_labeled = -1;
	int hocr_itr_labeled = -1;
	int fixetal_labeled = -1;
	int fixetal_itr_labeled = -1;
	int optimal_labeled = -1;
	int generators_labeled = -1;
	int heur_labeled = -1;
	int packing_labeled = -1;

	int generators_labeled_new = -1;
	int generators_labeled_new_10 = -1;

	real m_bound = 100;
	real hocr_bound = 100;
	real hocr_itr_bound = 100;
	real fixetal_bound = 100;
	real fixetal_itr_bound = 100;
	real optimal_bound = 100;
	real generators_bound = 100;
	real heur_bound = 100;
	real lp_bound = 100;
	real packing_bound = 100;
	real packing_itr_bound = 100;

	real generators_bound_new =100;
	real generators_bound_new_10 = 100;

	double m_time = -1;
	double hocr_time = -1;
	double hocr_itr_time = -1;
	double fixetal_time = -1;
	double fixetal_itr_time = -1;
	double optimal_time = -1;
	double generators_time = -1;
	double heur_time = -1;
	double packing_time = -1;

	double generators_time_new = 0;
	double generators_time_new_10 = 0;
	double time_gen10_solve = 0;
	double time_gen_solve = 0;

	if (do_exhaustive && n>30) {
		cout << "Not using exhaustive search for n=" << n << endl;
		do_exhaustive = false;
	}

	if (do_exhaustive) {
		cout << WHITE << "WHITE" << NORMAL << " is global optimum" << endl;
	}
	if (do_packing) {
		cout << BROWN << "BROWN" << NORMAL << " is using vertex packing" << endl;
	}
	if (do_lprelax) {
		cout << NORMAL << "GRAY" << NORMAL << " is an LP relaxation (Rhys form)" << endl;
	}
	if (do_m) {
		cout << DKRED << "DKRED" << NORMAL << " is M-reductions" << endl;
	}
	cout << RED << "RED" << NORMAL << " is HOCR" << endl;
	if (iterate_reduction_methods) {
		cout << DKRED << "DKRED" << NORMAL << " is iterated HOCR" << endl;
	}
	if (do_fixetal) {
		cout << BLUE << "BLUE" << NORMAL << " is Fix et al. from ICCV 2011" << endl;
		if (iterate_reduction_methods) {
			cout << DKBLUE << "DKBLUE" << NORMAL << " is Fix et al. iterated" << endl;
		}
	}
	if (do_optimal) {
		cout << GREEN << "GREEN" << NORMAL << " is generalized roof duality (using LP)" << endl;
	}
	if (do_generators) {
		cout << DKGREEN << "DKGREEN" << NORMAL << " is generalized roof duality with generators (using LP)" << endl;
	}
	if (do_heuristic) {
		cout << YELLOW << "YELLOW" << NORMAL << " is heuristic submodular relaxation" << endl;
	}
	cout << endl;


	// Start testing the different methods.
	// If there is a failure inside this try block, the polynomial
	// will be saved to a temporary file for debugging.
	try {

		//Holds optimal solutions
		vector< vector<label> > optimal_solutions;
		real optimum = -1;

		if (do_exhaustive) {
			ASSERT(n<=30); // Otherwise too big
			cout << "Exhaustive search: " << endl;

			vector<label> x(n,0);
			optimum = pb.eval(x);
			while (true) {
				x[0]++;
				int i=0;
				while (x[i]>1) {
					x[i]=0;
					i++;
					if (i==n) {
						break;
					}
					x[i]+=1;
				}
				if (i==n) {
					break;
				}

				real energy = pb.eval(x);
				if (energy < optimum) {
					optimum = energy;
				}
			}

			for (int i=0;i<n;++i) { x[i] = 0; }

			x[0]=-1;
			while (true) {
				x[0]++;
				int i=0;
				while (x[i]>1) {
					x[i]=0;
					i++;
					if (i==n) {
						break;
					}
					x[i]+=1;
				}
				if (i==n) {
					break;
				}

				real energy = pb.eval(x);
				if (energy == optimum) {
					optimal_solutions.push_back(x);
					print_info("Global minimum",x,energy,n,WHITE);
				}
			}
			cout << endl;
		}

		if (do_packing) {
			packing_bound = -1000000000;
			for (int num_tests=1;num_tests<=10;++num_tests) {

				double bound = 0;
				double first_bound = 0;
				int labeled = 0;
				bool should_continue;
				Petter::PseudoBoolean<real> f = pb;
				vector<label> x(n,0);

				cout << " *** " << endl;
				bool first_time = true;

				packing_time = 0;
				do {
					// Create posiform for -pb
					Posiform<real,4> phi(f,true);

					start();
					// Maximize it
					real bound = -phi.maximize(x);
					int new_labeled=0;
					for (int i=0;i<n;++i) {
						if (x.at(i) >= 0) {
							new_labeled++;
						}
					}

					f.reduce(x);
					packing_time += stop();

					should_continue = new_labeled > labeled;
					labeled = new_labeled;
					if (labeled == n ) {
						//Nothing more to do
						should_continue = false;
					}
					if (first_time) {
						first_bound = bound;
						first_time = false;
					}

					// Print solution
					print_info("Vertex packing",x,bound,new_labeled,BROWN,packing_time);
				} while (should_continue);


				if (first_bound > packing_bound) {
					packing_bound = first_bound;
				}
				if (bound > packing_itr_bound) {
					packing_itr_bound = bound;
				}
				if (labeled > packing_labeled) {
					packing_labeled = labeled;
				}

				// If we know the optimal solution, we can verify
				// that the persistencies are correct
				if (do_exhaustive) {
					check_persistency(optimal_solutions, x, labeled);
					check_bound(optimum, packing_bound);
				}

			}
			if (verbose) {
				cout << "Vertex packing time : " << BROWN << packing_time << NORMAL << endl;
			}
			cout << endl;
		}

		if (do_lprelax) {
			vector<label> x_lp(n);
			lp_bound = pb.minimize_lp(x_lp,verbose);
			print_info("LP relaxation",x_lp,lp_bound,0,NORMAL);
			cout << endl;
			if (do_exhaustive) {
				check_bound(optimum, lp_bound);
			}
		}

		if (do_m) {
			vector<label> x(n,0);
			start();
			m_bound = pb.minimize(x, m_labeled, M_reduction);
			m_time = stop();

			print_info("M-red.",x,m_bound,m_labeled,DKRED,m_time);
			cout << endl;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, m_labeled);
				check_bound(optimum, m_bound);
			}
		}

		Petter::PseudoBoolean<real> f_hocrreduced;

		if (do_hocr) {

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f = pb;
			vector<label> x(n,0);

			const Petter::Color* COL = &RED;

			do {
				iters++;

				int new_labeled = 0;
				start();
				bound = f.minimize_reduction(x,new_labeled);
				double t_minimize = stop();

				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n || !iterate_reduction_methods) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (iters == 1) {
					hocr_bound = bound;
					hocr_labeled = new_labeled;
					hocr_time = t_minimize;
					hocr_itr_time = t_minimize + t_reduce;
					COL = &RED;
				}
				else {
					COL = &DKRED;
					hocr_itr_time += t_minimize + t_reduce;
				}

				print_info("HOCR",x,bound,labeled,*COL,hocr_itr_time);
				if (verbose) {
					cout << "time (minimize) : " << *COL << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << *COL << t_reduce <<  NORMAL << endl;
				}

			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
			}
			cout << endl;

			hocr_itr_bound = bound;
			hocr_itr_labeled = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, hocr_itr_labeled);
				check_bound(optimum, hocr_bound);
				check_bound(optimum, hocr_itr_bound);
			}

			f_hocrreduced = f;
		}


		if (do_fixetal) {

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f = pb;
			vector<label> x(n,0);

			const Petter::Color* COL = &DKBLUE;

			do {
				iters++;

				int new_labeled = 0;
				start();
				bound = f.minimize_reduction_fixetal(x,new_labeled);
				double t_minimize = stop();

				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n || !iterate_reduction_methods) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (iters == 1) {
					fixetal_bound = bound;
					fixetal_labeled = new_labeled;
					fixetal_time = t_minimize;
					fixetal_itr_time = t_minimize + t_reduce;
					COL = &BLUE;
				}
				else {
					COL = &DKBLUE;
					fixetal_itr_time += t_minimize + t_reduce;
				}

				print_info("Fix et al.",x,bound,labeled,*COL,fixetal_itr_time);
				if (verbose) {
					cout << "time (minimize) : " << *COL << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << *COL << t_reduce <<  NORMAL << endl;
				}

			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
			}
			cout << endl;

			fixetal_itr_bound = bound;
			fixetal_itr_labeled = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, fixetal_itr_labeled);
				check_bound(optimum, fixetal_bound);
				check_bound(optimum, fixetal_itr_bound);
			}
		}


		if (do_optimal) {

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f;
			vector<label> x(n,0);

			if (cmd_line.find("-usehocr") != cmd_line.end()) {
				//Start at the HOCR solution (for speed)
				f = f_hocrreduced;
				cout << "USING HOCR REDUCED" << endl;
			}
			else {
				//Default
				f = pb;
			}

			optimal_time = 0;
			do {
				iters++;

				Petter::SymmetricPseudoBoolean<real> spb;
				start();
				spb.create_lp(f);
				double t_create = stop();

				int new_labeled = 0;
				start();
				bound = spb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (verbose) {
					cout << "Relaxation g : " << spb << endl;
				}

				print_info("GRD solution",x,bound,labeled,GREEN);
				if (verbose) {
					cout << "time (create)   : " << GREEN << t_create <<  NORMAL << endl;
					cout << "time (minimize) : " << GREEN << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << GREEN << t_reduce <<  NORMAL << endl;
				}

				optimal_time += t_create + t_minimize + t_reduce;

			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
			}
			cout << endl;

			optimal_bound = bound;
			optimal_labeled = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, optimal_labeled);
				check_bound(optimum, optimal_bound);
			}
		}
		cout << "111111" << endl;
		if (do_generators) {

			

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f;
			vector<label> x(n,0);

			if (cmd_line.find("-usehocr") != cmd_line.end()) {
				//Start at the HOCR solution (for speed)
				f = f_hocrreduced;
				cout << "USING HOCR REDUCED" << endl;
			}
			else {
				//Default
				f = pb;
			}

			Generators<real> generators(cmd_line["-generators-file"]);
			generators_time = 0;
			do {
				iters++;

				//TODO: Reading the file every iteration is not optimal
				GeneratorPseudoBoolean<real> gpb(generators);
				start();
				gpb.create_lp(f);

				double t_create = stop();

				int new_labeled = 0;
				start();
				bound = gpb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (verbose) {
					//cout << "Relaxation g : " << gpb << endl;
				}

				print_info("Gener. solution",x,bound,labeled,DKGREEN);
				if (verbose) {
					cout << "time (create)   : " << DKGREEN << t_create <<  NORMAL << endl;
					cout << "time (minimize) : " << DKGREEN << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << DKGREEN << t_reduce <<  NORMAL << endl;
				}

				generators_time += t_create + t_minimize + t_reduce;
				time_gen_solve += t_minimize + t_reduce;
			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;

				if (abs(bound) > 1e-5 && abs(pb.eval(x) - bound) / abs(bound) > 1e-5) {
					throw runtime_error("For all variables assigned, min g should be equal to min f.");
				}
			}
			cout << endl;

			generators_bound = bound;
			generators_labeled = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {

			

				check_persistency(optimal_solutions, x, generators_labeled);
				check_bound(optimum, generators_bound);
			}
		}
		cout << "111111" << endl;

		//////////////////////////////////
		//new generators!  allgen4.txt////
		//////////////////////////////////

		bool new_generators = true;
		if (new_generators) {

			cout << "-------Generators allgen4.txt---------" << endl;
			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f;
			vector<label> x(n,0);

			if (cmd_line.find("-usehocr") != cmd_line.end()) {
				//Start at the HOCR solution (for speed)
				f = f_hocrreduced;
				cout << "USING HOCR REDUCED" << endl;
			}
			else {
				//Default
				f = pb;
			}

			Generators<real> generators("generators/allgen4.txt");
			generators_time_new = 0;
			do {
				cout << "------------------" << endl;
				cout << f << endl;
				
				iters++;

				//TODO: Reading the file every iteration is not optimal
				GeneratorPseudoBoolean<real> gpb(generators);
			//	gpb.save_to_file("g-error.txt");
				start();
				gpb.create_lp(f);
				double t_create = stop();

				int new_labeled = 0;
				start();
				bound = gpb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (verbose) {
					//cout << "Relaxation g : " << gpb << endl;
				}

				print_info("Gener. solution",x,bound,labeled,DKGREEN);
				if (verbose) {
					cout << "time (create)   : " << DKGREEN << t_create <<  NORMAL << endl;
					cout << "time (minimize) : " << DKGREEN << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << DKGREEN << t_reduce <<  NORMAL << endl;
				}

				generators_time_new += t_create + t_minimize + t_reduce;
				cout << "ITERS = " <<GREEN <<iters << NORMAL << endl;
			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;

				if (abs(bound) > 1e-5 && abs(pb.eval(x) - bound) / abs(bound) > 1e-5) {
					throw runtime_error("For all variables assigned, min g should be equal to min f.");
				}
			}
			cout << endl;

			generators_bound_new = bound;
			generators_labeled_new = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, generators_labeled_new);
				check_bound(optimum, generators_bound);
				
			}
		}
		
		/////////////////////////////////////////
		//Minimize 3!///////////////////////////
		///////////////////////////////////////

		bool minimize_3 = true;
		if (minimize_3) {

			cout << "-------Minimize_3---------" << endl;
			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f;
			vector<label> x(n,0);

			if (cmd_line.find("-usehocr") != cmd_line.end()) {
				//Start at the HOCR solution (for speed)
				f = f_hocrreduced;
				cout << "USING HOCR REDUCED" << endl;
			}
			else {
				//Default
				f = pb;
			}

			Generators<real> generators("generators/allgen4_3.txt");
			generators_time_new = 0;
			do {
				cout << "------------------" << endl;
				cout << f << endl;
				
				iters++;

				//TODO: Reading the file every iteration is not optimal
				GeneratorPseudoBoolean<real> gpb(generators);
			//	gpb.save_to_file("g-error.txt");
				start();
				gpb.create_lp(f);
				double t_create = stop();

				int new_labeled = 0;
				start();
				bound = gpb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (verbose) {
					//cout << "Relaxation g : " << gpb << endl;
				}

				print_info("Gener. solution",x,bound,labeled,DKGREEN);
				if (verbose) {
					cout << "time (create)   : " << DKGREEN << t_create <<  NORMAL << endl;
					cout << "time (minimize) : " << DKGREEN << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << DKGREEN << t_reduce <<  NORMAL << endl;
				}

				generators_time_new += t_create + t_minimize + t_reduce;
				cout << "ITERS = " <<GREEN <<iters << NORMAL << endl;
			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;

				if (abs(bound) > 1e-5 && abs(pb.eval(x) - bound) / abs(bound) > 1e-5) {
					throw runtime_error("For all variables assigned, min g should be equal to min f.");
				}
			}
			cout << endl;

			generators_bound_new = bound;
			generators_labeled_new = labeled;


			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, generators_labeled_new);
				check_bound(optimum, generators_bound);
				
			}
		}






		

		///////////////////////////////////////////////
		//new generators!  generators_10classes.txt////
		///////////////////////////////////////////////
		//bool new_generators = true;



		if (new_generators) {
			cout << "-------Generators 10 classes---------" << endl;
			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f;
			vector<label> x(n,0);

			if (cmd_line.find("-usehocr") != cmd_line.end()) {
				//Start at the HOCR solution (for speed)
				f = f_hocrreduced;
				cout << "USING HOCR REDUCED" << endl;
			}
			else {
				//Default
				f = pb;
			}

			Generators<real> generators("generators/generators_10classes.txt");
			generators_time_new = 0;
			do {
				cout << "------------------" << endl;
				cout << f << endl;
				
				iters++;

				//TODO: Reading the file every iteration is not optimal
				GeneratorPseudoBoolean<real> gpb(generators);
			//	gpb.save_to_file("g-error.txt");
				start();
				gpb.create_lp(f);
				double t_create = stop();

				int new_labeled = 0;
				start();
				bound = gpb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (verbose) {
					//cout << "Relaxation g : " << gpb << endl;
				}

				print_info("Gener. solution",x,bound,labeled,DKGREEN);
				if (verbose) {
					cout << "time (create)   : " << DKGREEN << t_create <<  NORMAL << endl;
					cout << "time (minimize) : " << DKGREEN << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << DKGREEN << t_reduce <<  NORMAL << endl;
				}

				generators_time_new_10 += t_create + t_minimize + t_reduce;
				time_gen10_solve +=t_minimize + t_reduce;
				cout << "ITERS = " <<GREEN <<iters << NORMAL << endl;
			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;

				if (abs(bound) > 1e-5 && abs(pb.eval(x) - bound) / abs(bound) > 1e-5) {
					throw runtime_error("For all variables assigned, min g should be equal to min f.");
				}
			}
			cout << endl;

			generators_bound_new_10 = bound;
			generators_labeled_new_10 = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, generators_labeled_new_10);
				check_bound(optimum, generators_bound);
			}
		}




		if (do_heuristic) {
			//
			// For the heuristic relaxations, we use
			// integer arithmetic
			//
			typedef int heurreal;

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f = pb;
			vector<label> x(n,0);

			heur_time = 0;

			do {
				iters++;

				Petter::SymmetricPseudoBoolean<heurreal> spb;
				start();
				spb.create_heuristic(f);
				double t_create = stop();

				int new_labeled = 0;
				start();
				bound = spb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}


				//Test
				if (n<=1000) {
					vector<label> u(n,0);
					vector<label> v(n,0);
					for (int iter=0; iter<=100; ++iter) {
						for (int i=0;i<n;++i) {
							u[i] = rand()%2;
							v[i] = 1-u[i];
						}
						real fval = f.eval(u);
						heurreal gval = spb.eval(u,v);
						real diff = fval - gval;
						if ( diff * diff > 1e-12 ) {
							cout << "f = " << fval << endl;
							cout << "g = " << gval << endl;
							throw runtime_error("f(x) =/= g(x,bar(x)) for heuristic");
						}
					}
				}

				start();
				f.reduce(x);
				double t_reduce = stop();

				if (verbose) {
					cout << "Relaxation g : " << spb << endl;
				}

				print_info("Heuristics",x,bound,labeled,YELLOW);
				if (verbose) {
					cout << "time (create)   : " << YELLOW << t_create <<  NORMAL << endl;
					cout << "time (minimize) : " << YELLOW << t_minimize <<  NORMAL << endl;
					cout << "time (reduce)   : " << YELLOW << t_reduce <<  NORMAL << endl;
				}


				heur_time += t_create + t_minimize + t_reduce;

			} while (should_continue);

			if (labeled == n) {
				cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
			}

			cout << endl;
			heur_bound = bound;
			heur_labeled = labeled;

			// If we know the optimal solution, we can verify
			// that the persistencies are correct
			if (do_exhaustive) {
				check_persistency(optimal_solutions, x, heur_labeled);
				check_bound(optimum, heur_bound);
			}
		}

		//
		// Branch and bound
		//
		if (cmd_line.find("-bb") != cmd_line.end() ) {

			ofstream fout("bb.log");
			fout << "m=" << cmd_line["-m"] << ", n=" << n << ", nterms=" << cmd_line["-nterms"] << endl;
			cout << "Starting branch and bound." << endl;


			if (do_generators) {
				test_branchandbound(pb, n,  GRD_gen, fout, DKGREEN);
			}
			if (do_optimal) {
				test_branchandbound(pb, n,  GRD, fout, GREEN);
			}
			if (do_heuristic) {
				test_branchandbound(pb, n, GRD_heur, fout, YELLOW);
			}
			if (do_fixetal) {
				test_branchandbound(pb, n, Fix, fout, BLUE);
			}
			if (cmd_line.find("-hocr") != cmd_line.end()) {
				test_branchandbound(pb, n,  HOCR, fout, RED);
			}

		}


	}
	catch (exception& e) {
		statusFailed();

		// Save erroneous polynomial to temporary file
		stringstream filename;
		if (std::getenv("TEMP")) {
			filename << std::getenv("TEMP") << "/";
		}
		int errorid = bind(uniform_int_distribution<int>(0,999), engine)();
		filename << "pb-error-" << errorid << ".txt";
		pb.save_to_file(filename.str());

		// Output filename saved to
		cerr << RED << "Solver error, polynomial saved to " << filename.str() << NORMAL << endl;

		// Write to error log
		stringstream errorlogfile;
		errorlogfile << args[0] << ".errorlog";
			// Get local time
			struct tm *current;
			time_t now;
			time(&now);
			current = localtime(&now);
		ofstream errorfile(errorlogfile.str(), ios::app);
		errorfile << 1900+current->tm_year << "-" << setw(2) << setfill('0')
		          << current->tm_mon+1 << '-' << current->tm_mday << ' '
		          << current->tm_hour << ':' << current->tm_min << ":"
		          << current->tm_sec << "  ";
		errorfile << "Solver error \"" << e.what()
		          << ",\" polynomial saved to " << filename.str() << endl;

		// Rethrow exception
		throw;
	}



	// Write to log file
	if (cmd_line.find("-nolog") == cmd_line.end()) {
		ofstream log("logfile.data", ios::app);
		log << n                << '\t' // 0
			<< 0                << '\t' // 1
			<< 0                << '\t' // 2
			<< hocr_labeled     << '\t' // 3
			<< hocr_itr_labeled << '\t' // 4
			<< optimal_labeled  << '\t' // 5
			<< heur_labeled     << '\t' // 6
			<< hocr_bound       << '\t' // 7
			<< hocr_itr_bound   << '\t' // 8
			<< optimal_bound    << '\t' // 9
			<< heur_bound       << '\t' // 10
			<< hocr_time        << '\t' // 11
			<< hocr_itr_time    << '\t' // 12
			<< optimal_time     << '\t' // 13
			<< heur_time        << '\t' // 14
			<< fixetal_labeled  << '\t' // 15
			<< fixetal_itr_labeled << '\t' // 16
			<< fixetal_bound    << '\t' // 17
			<< fixetal_itr_bound<< '\t' // 18
			<< fixetal_time     << '\t' // 19
			<< fixetal_itr_time << '\t' // 20
			<< generators_labeled<<'\t' // 21
			<< generators_bound << '\t' // 22
			<< generators_time  << '\t' // 23
			<< packing_labeled  <<'\t'  // 24
			<< packing_bound    << '\t' // 25
			<< packing_time     << '\t' // 26
			<< m_labeled        << '\t' // 27
			<< m_bound          << '\t' // 28
			<< m_time           << '\t' // 29

			<< generators_labeled_new<<'\t' // 30
			<< generators_bound_new << '\t' // 31
			<< generators_time_new  << '\t' // 32

			<< generators_labeled_new_10<<'\t' // 33
			<< generators_bound_new_10 << '\t' // 34
			<< generators_time_new_10  << '\t' // 35


			<< time_gen10_solve << '\t'  //36
			<< time_gen_solve << '\t'    //37

            << endl;
	}
	if(((generators_bound_new_10 - generators_bound) < -0.5 )){    //  || ((-generators_bound_new_10 + generators_bound) < -0.0001)){
		cout << "value1: " << (generators_bound_new_10 - generators_bound)<< endl;
		cout << "value2: " << (-generators_bound_new_10 + generators_bound)<< endl;
		cout << "old energy: " << generators_bound <<  "   10_gen energy: " << generators_bound_new_10 <<endl;
		pb.save_to_file("diff_energy.txt");
	}



	if (generators_labeled_new > generators_labeled)  {
		cout << "new: " <<  generators_labeled_new << " old: " << generators_labeled << endl;
		pb.save_to_file("diff_new_better.txt");
	}
	if (generators_labeled_new < generators_labeled)  {
		cout << "new: " <<  generators_labeled_new << " old: " << generators_labeled << endl;
		pb.save_to_file("diff_old_better.txt");
	}
	if (generators_labeled_new_10 < generators_labeled_new)  {
		cout << "new_10: " <<  generators_labeled_new_10 << " old: " << generators_labeled_new << endl;
		pb.save_to_file("diff_10_worse.txt");
	}

	if (generators_labeled_new_10 > generators_labeled_new)  {
		cout << "new_10: " <<  generators_labeled_new_10 << " old: " << generators_labeled_new << endl;
		pb.save_to_file("diff_10_better.txt");
	}


	//cin.get();

	return 0;
}

int main(int argc, char** argv)
{
	using namespace std;
	using namespace Petter;
#ifdef WIN32 //Also defined under Windows-x64
	system("cls");
#else
	system("clear");
#endif

	try {
		return main_program(argc,argv);
	}
	catch (runtime_error& e) {
		statusFailed();
		cerr << RED << "Run-time error : " << e.what() << NORMAL << endl;
		return 100;
	}
	catch (bad_alloc&) {
		statusFailed();
		cerr << RED << "Out of memory." << NORMAL << endl;
		return 101;
	}
	catch (exception& e) {
		statusFailed();
		cerr << RED << "Exception : " << e.what() << NORMAL << endl;
		return 102;
	}
	catch (...) {
		cerr << RED << "Unknown error" << NORMAL << endl;
		return 103;
	}
}



