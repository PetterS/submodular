
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


/*
void f() 
{
	typedef std::tuple<int,int> Pair;
	std::map<Pair,int> m;
	Pair p(1,2);
	m[p] = 1;
}
*/



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
		cerr << "  " << args[0] << " -heuristic                       : use heuristics as well" << endl;
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
		}

		//for (auto itr=clauses.begin(); itr != clauses.end(); ++itr) {
		//	vector<variable>& vars = *itr;

		//	for (auto itr2=vars.begin(); itr2 != vars.end(); ++itr2) {
		//		int  i   = itr2->first;
		//		bool neg = itr2->second;
		//		if (neg) {
		//			cout << "not(";
		//		}
		//		cout << i;
		//		if (neg) {
		//			cout << ")";
		//		}
		//		cout << " ";
		//	}
		//	cout << endl;
		//}


		m = 3; // Only 3-SAT for now

		for (auto itr=clauses.begin(); itr != clauses.end(); ++itr) {
			vector<variable>& vars = *itr;
			ASSERT(vars.size() == unsigned(m)); // Only 3-SAT for now

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
	vector<label> x(n);



	////////////////////////////
	// Solve using heuristics //
	////////////////////////////
	if  (cmd_line.find("-heuristic") != cmd_line.end() ) {

		cout << DKRED << "DKRED" << NORMAL << " is HOCR" << endl;
		cout << RED << "RED" << NORMAL << " is iterated HOCR" << endl;
		if (cmd_line.find("-nolp") == cmd_line.end()) {
			cout << GREEN << "GREEN" << NORMAL << " is LP optimal relaxation" << endl;
		}
		cout << YELLOW << "YELLOW" << NORMAL << " is heuristic relaxation" << endl;
		cout << endl;

		Petter::PseudoBoolean pb2 = pb;
		Petter::PseudoBoolean f   = pb;

		// For timing
		clock_t t_raw;
		auto start = [&t_raw]() { t_raw = clock(); };
		auto stop  = [&t_raw]() -> double { return double(clock()-t_raw) / double(CLOCKS_PER_SEC); };

		vector<label> x(n,0),x1(n,0),x2(n,0);
		int iters = 0;
		int reduction_labeled = 0;
		double f_bound = 0;
		double f1_bound = 0;
		int labeled0 = 0;
		bool should_continue;

		const Petter::Color* COL = &RED;

		do {
			iters++;

			int new_labeled = 0;
			start();
			f1_bound = pb.minimize_reduction(x,new_labeled);
			double t_minimize = stop();

			should_continue = new_labeled > labeled0;
			labeled0 = new_labeled;
			if (labeled0 == n) {
				//Nothing more to do
				should_continue = false;
			}

			start();
			pb.reduce(x);
			double t_reduce = stop();

			if (iters == 1) {
				f_bound = f1_bound;
				reduction_labeled = new_labeled;
				COL = &DKRED;
			}
			else {
				COL = &RED;
			}

			cout << "labeled : " << *COL << labeled0 << NORMAL << endl;
			cout << "f_bound : " << *COL << f_bound << NORMAL << endl;
			cout << "time (minimize) : " << *COL << t_minimize <<  NORMAL << endl;
			cout << "time (reduce)   : " << *COL << t_reduce <<  NORMAL << endl;
			cout << endl;

		} while (should_continue);

		if (labeled0 == n) {
			cout << "Global minimum : " << WHITE << f.eval(x) << NORMAL << endl;
			cout << endl;
		}

		int labeled1 = 0;
		iters = 0;
		Petter::real g1_bound;
		if (cmd_line.find("-nolp") == cmd_line.end()) {
			do {
				iters++;

				Petter::SymmetricPseudoBoolean spb;
				start();
				spb.create_lp(pb);
				double t_create = stop();

				int new_labeled = 0;
				start();
				g1_bound = spb.minimize(x, new_labeled);
				double t_minimize = stop();
				should_continue = new_labeled > labeled1;
				labeled1 = new_labeled;
				if (labeled1 == n) {
					//Nothing more to do
					should_continue = false;
				}

				start();
				pb.reduce(x);
				double t_reduce = stop();

				cout << "Labeled : "<< GREEN << labeled1 << NORMAL << endl;
				cout << "Bound   : " << GREEN << g1_bound << NORMAL << endl;
				cout << "time (create)   : " << GREEN << t_create <<  NORMAL << endl;
				cout << "time (minimize) : " << GREEN << t_minimize <<  NORMAL << endl;
				cout << "time (reduce)   : " << GREEN << t_reduce <<  NORMAL << endl;
				cout << endl;

			} while (should_continue);

			if (labeled1 == n) {
				cout << "Global minimum : " << WHITE << f.eval(x) << NORMAL << endl;
				cout << endl;
			}
		}

		

		int labeled2 = 0;
		iters = 0;
		should_continue = false;
		Petter::real g2_bound;
		do {
			iters++;

			Petter::SymmetricPseudoBoolean spb;
			start();
			spb.create_heuristic(pb2);
			double t_create = stop();

			int new_labeled = 0;
			start();
			g2_bound = spb.minimize(x, new_labeled);
			double t_minimize = stop();
			should_continue = new_labeled > labeled2;
			labeled2 = new_labeled;
			if (labeled2 == n) {
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
				double f = pb2.eval(u);
				double g = spb.eval(u,v);
				if ( absolute(f-g) > 1e-6 ) {
					cout << "f = " << f << endl;
					cout << "g = " << g << endl;
					throw runtime_error("f(x) =/= g(x,bar(x)) for heuristic");
				}
			}

			start();
			pb2.reduce(x);
			double t_reduce = stop();

			cout << "Labeled : "<< YELLOW << labeled2 << NORMAL << endl;
			cout << "Bound   : " << YELLOW << g2_bound << NORMAL << endl;
			cout << "time (create)   : " << YELLOW << t_create <<  NORMAL << endl;
			cout << "time (minimize) : " << YELLOW << t_minimize <<  NORMAL << endl;
			cout << "time (reduce)   : " << YELLOW << t_reduce <<  NORMAL << endl;	
			cout << endl;

		} while (should_continue);

		if (labeled2 == n) {
			cout << "Global minimum : " << WHITE << f.eval(x) << NORMAL << endl;
		}
		
		// Write to log file
		ofstream log("logfile_heuristic.data", ios::app);
		log << n << '\t' << nterms << '\t' << griddim << '\t'
			<< reduction_labeled << '\t' << labeled1 << '\t' << labeled2 << '\t'
			<< f_bound           << '\t' << g1_bound << '\t' << g2_bound  << endl;

		//cin.get();

		return 0;
	}
	else  {


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
	}
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
	srand(unsigned(time(0)));

	try {
		int ret;
		
		//for (int iter=0; iter <= 100; ++iter) {
			ret = main_program(argc,argv);
		//}

		return ret;
	}
	catch (runtime_error& e) {
		statusFailed();
		cerr << RED << "Run-time error : " << e.what() << NORMAL << endl;
	}
	catch (bad_alloc&) {
		statusFailed();
		cerr << RED << "Out of memory." << NORMAL << endl;
	}
	catch (exception& e) {
		statusFailed();
		cerr << RED << "Exception : " << e.what() << NORMAL << endl;
	}
	catch (...) {
		cerr << RED << "Unknown error" << NORMAL << endl;
	}
	return 1;
}




