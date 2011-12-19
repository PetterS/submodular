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
using namespace Petter;


// Random number generator
namespace {
	#ifdef WIN32
		// Time Stamp Counter gives better seeds than seconds
		// when many small problems are generated consecutively
		mt19937 engine(__rdtsc()); 
	#else
		mt19937 engine(unsigned(time(0)));
	#endif
}

// Functions running some quick tests
// In: submodular_tests.cpp
template<typename real> void test_pseudoboolean();
void test_minimize();
template<typename real> void test_posiform();

//Simple routine for conversion of strings
//used for the command line
template <typename T>
T convert_string(const std::string s) 
{
	std::istringstream is(s);
	T result;
	is >> result;
	if (!is) {
		throw std::runtime_error("conversion of \"" + s + "\" failed"); 
	}
	return result;
}

template <typename T>
T absolute(const T t) 
{
	return t < 0 ? -t : t;
}

void print_x(const std::vector<label>& x)
{
	int n = x.size();
	for (int i=0;i<n&&i<20;++i) {
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
void print_info(std::string name, const std::vector<label>& x, real bound, int labeled, Petter::Color color)
{
	using namespace std;
	cout << left << setw(15) << name << "f(";
	print_x(x);
	cout << ") = " << color << right << setw(8) << bound << NORMAL;
	cout << ",   labeled : " << color << labeled << NORMAL << endl;
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
		test_posiform<double>();
		statusTry("Testing posiform (int)...");
		test_posiform<int>();
		statusOK();

		cerr << "Possible choices : " << endl;
		cerr << "  " << args[0] << " -m <int> -n <int> -nterms <int>  : runs random examples" << endl;
		cerr << "  " << args[0] << " -m <int> -example                : examples from paper" << endl;
		cerr << "  " << args[0] << " -file <str>                      : read polynomial from file" << endl;
		cerr << "  " << args[0] << " -sat <str>                       : read SAT problem from file" << endl;
		cerr << endl;
		cerr << "    -fixetal                         : use reductions from Fix et al." << endl;
		cerr << "    -optimal                         : use linear programming" << endl;
		cerr << "    -heuristic                       : use heuristics" << endl;
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


	int m = convert_string<int>(cmd_line["-m"]);
	int n = convert_string<int>(cmd_line["-n"]);
	int nterms = convert_string<int>(cmd_line["-nterms"]);
	int griddim = 0;



	bool do_hocr = true;
	bool do_optimal = false;
	bool do_fixetal = false;
	bool do_exhaustive = false;
	bool do_heuristic = false;
	bool do_lprelax = false;
	bool do_packing = false;
	bool iterate_reduction_methods = false;
	if (cmd_line.find("-optimal") != cmd_line.end()) {
		do_optimal = true;
	}
	if (cmd_line.find("-fixetal") != cmd_line.end()) {
		do_fixetal = true;
	}
	if (cmd_line.find("-heuristic") != cmd_line.end()) {
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

	////////////////////
	// Read from file //
	////////////////////
	if  (cmd_line.find("-file") != cmd_line.end() ) {
		pb = PseudoBoolean<real>(cmd_line["-file"]);
		n = pb.nvars();
	}
	///////////////
	//  Dat file //
	///////////////
	else if (cmd_line.find("-dat") != cmd_line.end() ) {
		cout << "File : " << cmd_line["-dat"] << endl;
		statusTry("Reading file...");
		ifstream fin(cmd_line["-dat"]);
		ASSERT(fin);

		uniform_real_distribution<double> distribution(1.0, 1.000001);
		//auto random_mod = bind(distribution, engine);

		fin >> n;
		for (int i=0;i<n;++i) {
			real E0,E1;
			fin >> E0 >> E1;
			//E0 *= random_mod();
			//E1 *= random_mod();
			pb.add_clique(i,E0,E1);
		}
		ASSERT(fin);
		int np;
		fin >> np;
		for (int c=0;c<np;++c) {
			int i,j;
			real E00,E01,E10,E11;
			fin >> i >> j >> E00 >> E01 >> E10 >> E11;
			i--;
			j--;
			pb.add_clique(i,j,E00,E01,E10,E11);
		}
		ASSERT(fin);
		int nt;
		fin >> nt;
		for (int c=0;c<nt;++c) {
			int i,j,k;
			real E000,E001,E010,E011,E100,E101,E110,E111;
			fin >> i >> j >> k >> E000 >> E001 >> E010 >> E011 >> E100 >> E101 >> E110 >> E111;
			i--;
			j--;
			k--;
			pb.add_clique(i,j,k,E000,E001,E010,E011,E100,E101,E110,E111);
		}
		ASSERT(fin);
		statusOK();
		cout << n << " variables, " << np << " quadratic terms and " << nt << " cubic terms." << endl;
	}
	/////////
	// SAT //
	/////////
	else if (cmd_line.find("-sat") != cmd_line.end() ) {
		ifstream fin(cmd_line["-sat"]);
		ASSERT(fin);
		string line;
		char cmd;
		do {
			fin >> cmd;
			if (cmd == 'c')
				getline(fin, line);
		} while (cmd=='c');
		ASSERT( cmd == 'p');
		string cnf;
		int nclauses;
		fin >> cnf >> n >> nclauses;
		ASSERT( cnf == "cnf");
		cout << "n        = " << n << endl;
		cout << "nclauses = " << nclauses << endl;

		typedef std::pair<int, bool> variable;
		vector< vector<variable> > clauses;

		m = 0;

		for (int c=1; c<=nclauses; ++c) {
			vector<variable> vars;
			int i;
			while (true) {
				fin >> i;
				if (i == 0) {
					break;
				}
				vars.push_back( make_pair(abs(i)-1, i < 0) );
			}

			sort(vars.begin(), vars.end());

			clauses.push_back( vars );
			
			m = max(m, (int)vars.size());
		}

		if (verbose) {
			for (auto itr=clauses.begin(); itr != clauses.end(); ++itr) {
				vector<variable>& vars = *itr;

				for (auto itr2=vars.begin(); itr2 != vars.end(); ++itr2) {
					int  i   = itr2->first;
					bool neg = itr2->second;
					if (neg) {
						cout << "not(";
					}
					cout << i;
					if (neg) {
						cout << ")";
					}
					cout << " ";
				}
				cout << endl;
			}
		}

		for (auto itr=clauses.begin(); itr != clauses.end(); ++itr) {
			vector<variable>& vars = *itr;
			ASSERT(vars.size() == 2 || vars.size() == 3); // Only 3-SAT or lower for now

			if (vars.size() == 3) {
				// If the clause is "xi or not(xj) or not(xk)"
				// then all assignments of (xi,xj,xk) except
				// (0,1,1) are OK.

				int eind = 0;
				if (vars.at(0).second) {
					eind += 4;
				}
				if (vars.at(1).second) {
					eind += 2;
				}
				if (vars.at(2).second) {
					eind += 1;
				}
				vector<real> E(8,0);
				E.at(eind) = 1;

				int i = vars.at(0).first;
				int j = vars.at(1).first;
				int k = vars.at(2).first;
				pb.add_clique(i,j,k, E);
			}
			else if (vars.size() == 2) {
				int i = vars.at(0).first;
				int j = vars.at(1).first;
				bool i_neg = vars.at(0).second;
				bool j_neg = vars.at(1).second;

				if (i_neg && j_neg) {
					pb.add_clique(i,j, 0,0,0,1);
				}
				else if (i_neg && !j_neg) {
					pb.add_clique(i,j, 0,0,1,0);
				}
				else if (!i_neg && j_neg) {
					pb.add_clique(i,j, 0,1,0,0);
				}
				else if (!i_neg && !j_neg) {
					pb.add_clique(i,j, 1,0,0,0);
				}
			}
			else {
				throw runtime_error("Invalid clause length");
			}
		}
	}
	/////////////////////////
	// Examples from paper //
	/////////////////////////
	else if (cmd_line.find("-example") != cmd_line.end() ) {
		
		if (m == 3) {
			n = 3;

			pb.add_monomial(0, -2.0);
			pb.add_monomial(1,  1.0);
			pb.add_monomial(2, -1.0);

			pb.add_monomial(0,1, 4.0);
			pb.add_monomial(0,2, 4.0);
			pb.add_monomial(1,2, -2.0);

			pb.add_monomial(0,1,2, -2.0);
			cout << endl << YELLOW << "Cubic example" << NORMAL << endl;
		}
		else if (m == 4) {
			n = 4;
			pb.add_monomial(0,  1.0);
			pb.add_monomial(2,  2.0);
			pb.add_monomial(3, -1.0);

			pb.add_monomial(0,3, 1.0);
			pb.add_monomial(1,2, 1.0);
			pb.add_monomial(2,3, -1.0);

			pb.add_monomial(0,1,2,3, 2.0);

			cout << endl << YELLOW << "Quartic example" << NORMAL << endl;
			cout << "Polynomial : " << pb << endl << endl;
		}
		else {
			throw runtime_error("m != {3,4}");
		}
	}
	////////////////////////////////
	// Generate random polynomial //
	////////////////////////////////
	else {
		ASSERT(n>=m);
		ASSERT(m>=2);
		ASSERT(nterms>0);

		statusTry("Generating random polynomial...");

		auto random_coef  = bind(uniform_int_distribution<int>(-100,100), engine);
		real upper = 100;
		if (submodular) {
			upper = 0;
		}
		auto random_coef2 = bind(uniform_int_distribution<int>(-100,upper), engine);
		auto random_index = bind(uniform_int_distribution<int>(0,n-1), engine);

		map< quad, bool> exists;
		map< int, bool > used;

		auto add_term = [&nterms,&used,&exists,&pb,&random_coef,&random_coef2,m,n](int i, int j, int k, int l) 
		{
			if (m<4) {
				l = 3*n;
			}
			if (m<3) {
				k = 2*n;
			}

			if (nterms > 0 && !exists[ make_quad(i,j,k,l) ]) {
				exists[ make_quad(i,j,k,l) ] = true;
				used[i] = true;
				used[j] = true;
				if (m >= 3) {
					used[k] = true;
				}
				if (m >= 4) {
					used[l] = true;
				}

				pb.add_monomial(i, random_coef());
				pb.add_monomial(j, random_coef());
				pb.add_monomial(i,j, random_coef2());
				if (m >= 3) {
					pb.add_monomial(k, random_coef());
					pb.add_monomial(i,k, random_coef2());
					pb.add_monomial(j,k, random_coef2());
					pb.add_monomial(i,j,k, random_coef2());
				}
				if (m >= 4) {
					pb.add_monomial(l, random_coef());
					pb.add_monomial(i,l, random_coef2());
					pb.add_monomial(j,l, random_coef2());
					pb.add_monomial(k,l, random_coef2());
					pb.add_monomial(i,j,l, random_coef2());
					pb.add_monomial(i,k,l, random_coef2());
					pb.add_monomial(j,k,l, random_coef2());
					pb.add_monomial(i,j,k,l, random_coef2());
				}

				//cout << i << ' ' << j << ' ' << k << ' ' << l << endl;

				nterms--;
			}
		}; 

		//First make sure every variable is used once
		for (int i=0; i<=n-m-1 && nterms>0; i+=m) {
			add_term(i,i+1,i+2,i+3);
		}
		int i = n-m;
		add_term(i,i+1,i+2,i+3);

		while (nterms > 0) {
			int i,j,k,l;
			do {
				i = random_index();
				j = random_index();

				if (m>=3) {
					k = random_index();
				}
				else {
					k = 2*n;
				}

				if (m>=4) {
					l = random_index();
				}
				else {
					l = 3*n;
				}
			} while (i>=j || i>=k || j>=k || k>=l);

			add_term(i,j,k,l);			
		}

		// Save to temporary file
		if (std::getenv("TEMP")) {
			string tmp = std::getenv("TEMP");
			pb.save_to_file(tmp + "/pb.txt");
		}

		statusOK();
		cout << "Number of used variables : " << used.size() << endl;
	}


	cout << "Polynomial : " << pb << endl;
	cout << WHITE;
	cout << "n = " << n << endl;
	cout << "m = " << m << endl << endl;
	cout << NORMAL;


	

	///////////////////////////////////
	// Solve using different methods //
	///////////////////////////////////

	

	int hocr_labeled = -1;
	int hocr_itr_labeled = -1;
	int fixetal_labeled = -1;
	int fixetal_itr_labeled = -1;
	int optimal_labeled = -1;
	int heur_labeled = -1;
	int packing_labeled = -1;

	real hocr_bound = 100;
	real hocr_itr_bound = 100;
	real fixetal_bound = 100;
	real fixetal_itr_bound = 100;
	real optimal_bound = 100;
	real heur_bound = 100;
	real lp_bound = 100;
	real packing_bound = 100;
	real packing_itr_bound = 100;

	double hocr_time = -1;
	double hocr_itr_time = -1;
	double fixetal_time = -1;
	double fixetal_itr_time = -1;
	double optimal_time = -1;
	double heur_time = -1;
  double packing_time = -1;

	if (do_exhaustive) {
		cout << WHITE << "WHITE" << NORMAL << " is global optimum" << endl;
	}
	if (do_packing) {
		cout << BROWN << "BROWN" << NORMAL << " is using vertex packing" << endl;
	}
	if (do_lprelax) {
		cout << NORMAL << "GRAY" << NORMAL << " is an LP relaxation (Rhys form)" << endl;
	}
	cout << DKRED << "DKRED" << NORMAL << " is HOCR" << endl;
	if (iterate_reduction_methods) {
		cout << RED << "RED" << NORMAL << " is iterated HOCR" << endl;
	}
	if (do_fixetal) {
		cout << DKBLUE << "DKBLUE" << NORMAL << " is Fix et al. from ICCV 2011" << endl;
		if (iterate_reduction_methods) {
			cout << BLUE << "BLUE" << NORMAL << " is Fix et al. iterated" << endl;
		}
	}
	if (do_optimal) {
		cout << GREEN << "GREEN" << NORMAL << " is optimal generalized roof duality (using LP)" << endl;
	}
	if (do_heuristic) {
		cout << YELLOW << "YELLOW" << NORMAL << " is heuristic submodular relaxation" << endl;
	}
	cout << endl;


  // For timing
	clock_t t_raw;
	auto start = [&t_raw]() { t_raw = clock(); };
	auto stop  = [&t_raw]() -> double { return double(clock()-t_raw) / double(CLOCKS_PER_SEC); };
	
	try {

		//Holds optimal solutions
		vector< vector<label> > optimal_solutions;

		if (do_exhaustive) {
			ASSERT(n<=30); // Otherwise too big
			cout << "Exhaustive search: " << endl;

			vector<label> x(n,0);
			real best_energy = pb.eval(x);
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
				if (energy < best_energy) {
					best_energy = energy;
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
				if (energy == best_energy) {

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
					print_info("Vertex packing",x,bound,new_labeled,BROWN);
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
					bool any_ok = false;
					for (auto itr = optimal_solutions.begin(); itr != optimal_solutions.end(); ++itr) {
						bool this_ok = true;
						for (int i=0;i<n;++i) {
							if (x.at(i) >= 0 && x.at(i) != itr->at(i)) {
								this_ok = false;
							}
						}
						if (this_ok) {
							any_ok = true;
							break;
						}
					}
					if (!any_ok) {
						// Persistency did not hold for this solution

						pb.save_to_file("packing-error.txt");
						throw runtime_error("Vertex packing solution without persistency");
					}
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
		}

	
		vector<label> x(n,0),x1(n,0),x2(n,0);

		

		Petter::PseudoBoolean<real> f_hocrreduced;

		if (do_hocr) {

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f = pb;

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
					COL = &DKRED;
				}
				else {
					COL = &RED;
					hocr_itr_time += t_minimize + t_reduce;
				}

				print_info("HOCR",x,bound,labeled,*COL);
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

			f_hocrreduced = f;
		}


		if (do_fixetal) {

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f = pb;

			const Petter::Color* COL = &BLUE;

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
					COL = &DKBLUE;
				}
				else {
					COL = &BLUE;
					fixetal_itr_time += t_minimize + t_reduce;
				}

				print_info("Fix et al.",x,bound,labeled,*COL);
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

			f_hocrreduced = f;
		}

		
		if (do_optimal) {

			int iters = 0;
			double bound = 0;
			int labeled = 0;
			bool should_continue;
			Petter::PseudoBoolean<real> f;

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
						heurreal fval = f.eval(u);
						heurreal gval = spb.eval(u,v);
						if ( absolute(fval-gval) > 1e-6 ) {
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
		}
	}
	catch (...) {
		// Save erroneous polynomial to temporary file
		if (std::getenv("TEMP")) {
			string tmp = std::getenv("TEMP");
			pb.save_to_file(tmp + "/pb-error.txt");
		}
		// Rethrow exception
		throw;
	}

	// Write to log file
	if (cmd_line.find("-nolog") == cmd_line.end()) {
		ofstream log("logfile.data", ios::app);
		log << n                << '\t' // 0 
			<< nterms           << '\t' // 1
			<< griddim          << '\t' // 2
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
			<< fixetal_itr_time << endl;// 20
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



