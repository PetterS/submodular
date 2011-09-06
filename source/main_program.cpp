
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
using namespace Petter;


// Random number generator
namespace {
	mt19937 engine(unsigned(time(0)&0xffffffff));
}

// Function running some quick tests
// In: submodular_tests.cpp
void test_pseudoboolean();


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

//
// Main program
//
int main_program(int num_args, char** args)
{
	using namespace std;
	using namespace Petter;

	if (num_args == 1) {
		// Run some tests
		statusTry("Testing...");
		test_pseudoboolean();
		statusOK();

		cerr << "Possible choices : " << endl;
		cerr << "  " << args[0] << " -m <int> -n <int> -nterms <int>  : runs random examples" << endl;
		cerr << "  " << args[0] << " -m <int> -example                : examples from paper" << endl;
		cerr << "  " << args[0] << " -file <str>                      : read polynomial from file" << endl;
		cerr << "  " << args[0] << " -sat <str>                       : read SAT problem from file" << endl;
		cerr << endl;
		cerr << "    -lp                              : use linear programming" << endl;
		cerr << "    -heuristic                       : use heuristics" << endl;
		cerr << endl;
		cerr << "    -verbose                         : print polynomials" << endl;
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

	Petter::PseudoBoolean pb;

	////////////////////
	// Read from file //
	////////////////////
	if  (cmd_line.find("-file") != cmd_line.end() ) {
		pb = PseudoBoolean(cmd_line["-file"]);
		n = pb.nvars();
	}
	///////////////
	//  Dat file //
	///////////////
	else if (cmd_line.find("-dat") != cmd_line.end() ) {
		statusTry("Reading file...");
		ifstream fin(cmd_line["-dat"]);
		ASSERT(fin);
		fin >> n;
		for (int i=0;i<n;++i) {
			real E0,E1;
			fin >> E0 >> E1;
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

		if (cmd_line.find("-verbose") != cmd_line.end()) {
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

		cout << endl << YELLOW << "Random polynomial" << NORMAL << endl;

		uniform_int_distribution<int> distribution(-100, 100);
		auto random_value = bind(distribution, engine);

		map< quad, bool> exists;
		for (int t=1; t <= nterms; ++t) {
			int i,j,k,l;
			do {
				i = rand()%n;
				j = rand()%n;

				if (m>=3) {
					k = rand()%n;
				}
				else {
					k = 2*n;
				}

				if (m>=4) {
					l = rand()%n;
				}
				else {
					l = 3*n;
				}
			} while (i>=j || i>=k || j>=k || k>=l || exists[ make_quad(i,j,k,l) ] );
			exists[ make_quad(i,j,k,l) ] = true;


			pb.add_monomial(i, random_value());
			pb.add_monomial(j, random_value());
			pb.add_monomial(i,j, random_value());
			if (m >= 3) {
				pb.add_monomial(k, random_value());
				pb.add_monomial(i,k, random_value());
				pb.add_monomial(j,k, random_value());
				pb.add_monomial(i,j,k, random_value());
			}
			if (m >= 4) {
				pb.add_monomial(l, random_value());
				pb.add_monomial(i,l, random_value());
				pb.add_monomial(j,l, random_value());
				pb.add_monomial(k,l, random_value());
				pb.add_monomial(i,j,l, random_value());
				pb.add_monomial(i,k,l, random_value());
				pb.add_monomial(j,k,l, random_value());
				pb.add_monomial(i,j,k,l, random_value());
			}
		}

		// Save to temporary file
		string tmp = std::getenv("TEMP");
		pb.save_to_file(tmp + "/pb.txt");

	}

	const int print_limit = 10;

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
	int lp_labeled = -1;
	int heur_labeled = -1;

	Petter::real hocr_bound = 100;
	Petter::real hocr_itr_bound = 100;
	Petter::real lp_bound = 100;
	Petter::real heur_bound = 100;

	cout << DKRED << "DKRED" << NORMAL << " is HOCR" << endl;
	cout << RED << "RED" << NORMAL << " is iterated HOCR" << endl;
	if (cmd_line.find("-lp") != cmd_line.end()) {
		cout << GREEN << "GREEN" << NORMAL << " is LP optimal relaxation" << endl;
	}
	if (cmd_line.find("-heuristic") != cmd_line.end()) {
		cout << YELLOW << "YELLOW" << NORMAL << " is heuristic relaxation" << endl;
	}
	cout << endl;

	

	// For timing
	clock_t t_raw;
	auto start = [&t_raw]() { t_raw = clock(); };
	auto stop  = [&t_raw]() -> double { return double(clock()-t_raw) / double(CLOCKS_PER_SEC); };

	vector<label> x(n,0),x1(n,0),x2(n,0);

	int iters = 0;
	double bound = 0;
	int labeled = 0;
	bool should_continue;

	bool run_hocr = true;
	if (run_hocr) {

		Petter::PseudoBoolean f = pb;

		const Petter::Color* COL = &RED;

		do {
			iters++;

			int new_labeled = 0;
			start();
			bound = f.minimize_reduction(x,new_labeled);
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

			if (iters == 1) {
				hocr_bound = bound;
				hocr_labeled = new_labeled;
				COL = &DKRED;
			}
			else {
				COL = &RED;
			}

			cout << "labeled : " << *COL << labeled << NORMAL << endl;
			cout << "f_bound : " << *COL << bound << NORMAL << endl;
			cout << "time (minimize) : " << *COL << t_minimize <<  NORMAL << endl;
			cout << "time (reduce)   : " << *COL << t_reduce <<  NORMAL << endl;
			cout << endl;

		} while (should_continue);

		if (labeled == n) {
			cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
			cout << endl;
		}

		hocr_itr_bound = bound;
		hocr_itr_labeled = labeled;
	}

		
	if (cmd_line.find("-lp") != cmd_line.end()) {

		Petter::PseudoBoolean f   = pb;

		iters = 0;
		labeled = 0;

		do {
			iters++;

			Petter::SymmetricPseudoBoolean spb;
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

			if (cmd_line.find("-verbose") != cmd_line.end()) {
				cout << "Relaxation g : " << spb << endl;
			}

			cout << "Labeled : "<< GREEN << labeled << NORMAL << endl;
			cout << "Bound   : " << GREEN << bound << NORMAL << endl;
			cout << "time (create)   : " << GREEN << t_create <<  NORMAL << endl;
			cout << "time (minimize) : " << GREEN << t_minimize <<  NORMAL << endl;
			cout << "time (reduce)   : " << GREEN << t_reduce <<  NORMAL << endl;
			cout << endl;

		} while (should_continue);

		if (labeled == n) {
			cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
			cout << endl;
		}

		lp_bound = bound;
		lp_labeled = labeled;
	}

		
	if (cmd_line.find("-heuristic") != cmd_line.end()) {

		Petter::PseudoBoolean f = pb;

		iters = 0;
		labeled = 0;

		do {
			iters++;

			Petter::SymmetricPseudoBoolean spb;
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
			vector<label> u(n,0);
			vector<label> v(n,0);
			for (int iter=0; iter<=100; ++iter) {
				for (int i=0;i<n;++i) {
					u[i] = rand()%2;
					v[i] = 1-u[i];
				}
				double fval = f.eval(u);
				double gval = spb.eval(u,v);
				if ( absolute(fval-gval) > 1e-6 ) {
					cout << "f = " << fval << endl;
					cout << "g = " << gval << endl;
					throw runtime_error("f(x) =/= g(x,bar(x)) for heuristic");
				}
			}

			start();
			f.reduce(x);
			double t_reduce = stop();

			if (cmd_line.find("-verbose") != cmd_line.end()) {
				cout << "Relaxation g : " << spb << endl;
			}

			cout << "Labeled : "<< YELLOW << labeled << NORMAL << endl;
			cout << "Bound   : " << YELLOW << bound << NORMAL << endl;
			cout << "time (create)   : " << YELLOW << t_create <<  NORMAL << endl;
			cout << "time (minimize) : " << YELLOW << t_minimize <<  NORMAL << endl;
			cout << "time (reduce)   : " << YELLOW << t_reduce <<  NORMAL << endl;	
			cout << endl;

		} while (should_continue);

		if (labeled == n) {
			cout << "Global minimum : " << WHITE << pb.eval(x) << NORMAL << endl;
		}

		heur_bound = bound;
		heur_labeled = labeled;
	}
		
	// Write to log file
	ofstream log("logfile.data", ios::app);
	log << n                << '\t' // 0 
		<< nterms           << '\t' // 1
		<< griddim          << '\t' // 2
		<< hocr_labeled     << '\t' // 3
		<< hocr_itr_labeled << '\t' // 4
		<< lp_labeled       << '\t' // 5
		<< heur_labeled     << '\t' // 6
		<< hocr_bound       << '\t' // 7
		<< hocr_itr_bound   << '\t' // 8
		<< lp_bound         << '\t' // 9
		<< heur_bound      << endl; // 10

	//cin.get();

	return 0;

#if 0


		////////////////////
		// Solve using LP //
		////////////////////

		int labeled, reduction_labeled;
		bool should_continue = false;
		double bound, reduction_bound;

		reduction_bound = pb.minimize_reduction(x,reduction_labeled);

		cout << "x = (";
		for (int i=0;i<n-1 && i<=print_limit;++i) {
			if (i < n-2 && i == print_limit) {
				cout << " ... ";
			}
			else {
				cout << x.at(i) << ", ";
			}
		}
		cout << x.at(n-1) << ")" << endl;
		cout << "Labeled : " << RED << reduction_labeled << NORMAL<< endl;
		cout << "Bound   : "<< RED << reduction_bound << endl << NORMAL<< endl;

		labeled = 0;
		int iters = 0;
		do {
			iters++;

			Petter::SymmetricPseudoBoolean spb;
			spb.create_lp(pb);

			cout << "Submodular relaxation : " << spb << endl;

			int new_labeled = 0;
			bound = spb.minimize(x, new_labeled);
			should_continue = new_labeled > labeled;
			labeled = new_labeled;
			if (labeled == n) {
				//Nothing more to do
				should_continue = false;
			}

			cout << "x = (";
			for (int i=0;i<n-1 && i<=print_limit;++i) {
				if (i < n-2 && i == print_limit) {
					cout << " ... ";
				}
				else {
					cout << x.at(i) << ", ";
				}
			}
			cout << x.at(n-1) << ")" << endl;
			cout << "Labeled : "<< GREEN << labeled << NORMAL << endl;
			cout << "Bound   : " << GREEN << bound << NORMAL << endl;

			pb.reduce(x);

		} while (should_continue);


		// Write to log file
		ofstream log("logfile.data", ios::app);
		log << n << '\t' << nterms << '\t' << reduction_labeled << '\t' << labeled << '\t' << griddim << '\t' << reduction_bound << '\t' << bound << '\t' << iters << endl;

		return 0;
	
#endif

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




