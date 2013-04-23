//
// Petter Strandmark 2012
// petter@maths.lth.se
//
// get_polynomial creates a polynomial from the command line.
//


#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include "Petter-Color.h"
#include "PseudoBoolean.h"

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

//Simple routine for conversion of strings
//used for the command line
template <typename T>
T convert_string(const std::string& s)
{
	std::istringstream is(s);
	T result;
	is >> result;
	if (!is) {
		throw std::runtime_error("conversion of \"" + s + "\" failed");
	}
	return result;
}

void get_polynomial(std::map<std::string,std::string>& cmd_line,
                    Petter::PseudoBoolean<double>* pb)
{
	using namespace std;
	using namespace Petter;

	typedef double real;

	int m = convert_string<int>(cmd_line["-m"]);
	int n = convert_string<int>(cmd_line["-n"]);
	int nterms = convert_string<int>(cmd_line["-nterms"]);
	int griddim = 0;
	bool submodular = cmd_line.find("-submodular") != cmd_line.end();
	bool verbose = cmd_line.find("-verbose") != cmd_line.end();

	////////////////////
	// Read from file //
	////////////////////
	if  (cmd_line.find("-file") != cmd_line.end() ) {
		*pb = PseudoBoolean<real>(cmd_line["-file"]);
		n = pb->nvars();
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
			pb->add_clique(i,E0,E1);
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
			pb->add_clique(i,j,E00,E01,E10,E11);
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
			pb->add_clique(i,j,k,E000,E001,E010,E011,E100,E101,E110,E111);
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
				pb->add_clique(i,j,k, E);
			}
			else if (vars.size() == 2) {
				int i = vars.at(0).first;
				int j = vars.at(1).first;
				bool i_neg = vars.at(0).second;
				bool j_neg = vars.at(1).second;

				if (i_neg && j_neg) {
					pb->add_clique(i,j, 0,0,0,1);
				}
				else if (i_neg && !j_neg) {
					pb->add_clique(i,j, 0,0,1,0);
				}
				else if (!i_neg && j_neg) {
					pb->add_clique(i,j, 0,1,0,0);
				}
				else if (!i_neg && !j_neg) {
					pb->add_clique(i,j, 1,0,0,0);
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

			pb->add_monomial(0, -2.0);
			pb->add_monomial(1,  1.0);
			pb->add_monomial(2, -1.0);

			pb->add_monomial(0,1, 4.0);
			pb->add_monomial(0,2, 4.0);
			pb->add_monomial(1,2, -2.0);

			pb->add_monomial(0,1,2, -2.0);
			cout << endl << YELLOW << "Cubic example" << NORMAL << endl;
		}
		else if (m == 4) {
			n = 4;
			pb->add_monomial(0,  1.0);
			pb->add_monomial(2,  2.0);
			pb->add_monomial(3, -1.0);

			pb->add_monomial(0,3, 1.0);
			pb->add_monomial(1,2, 1.0);
			pb->add_monomial(2,3, -1.0);

			pb->add_monomial(0,1,2,3, 2.0);

			cout << endl << YELLOW << "Quartic example" << NORMAL << endl;
			cout << "Polynomial : " << pb << endl << endl;
		}
		else {
			throw runtime_error("m != {3,4}");
		}
	}
	////////////////////////////////
	// Read previously generated  //
	////////////////////////////////
	else if (cmd_line.find("-prev") != cmd_line.end()) {
		string tmp = std::getenv("TEMP");
		*pb = PseudoBoolean<real>(tmp + "/pb.txt");
		n = pb->nvars();
		m = 4; //TODO
	}
	else if (cmd_line.find("-kolmogorov") != cmd_line.end()) {
		// Example from Kolmogorov's 2012 DAM paper
		real E[16] = {3, 2, 4, 10, 2, 12, 13, 12, 1, 3, 0, 12, 7, 10, 12, 14};
		vector<real> Evec(E, E+16);
		pb->add_clique(0,1,2,3, Evec);
		n = 4;
		m = 4;
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
		auto random_index = bind(uniform_int_distribution<int>(0,n-1), engine);

		map< quad, bool> exists;
		map< int, bool > used;

		auto add_term = [&nterms,&used,&exists,&pb,&random_coef,m,n]
		                (int i, int j, int k, int l)
		{
			if (m<4) {
				l = 3*n;
			}
			if (m<3) {
				k = 2*n;
			}

			if (nterms > 0 && !exists[ Petter::make_quad(i,j,k,l) ]) {
				exists[ Petter::make_quad(i,j,k,l) ] = true;
				used[i] = true;
				used[j] = true;
				if (m >= 3) {
					used[k] = true;
				}
				if (m >= 4) {
					used[l] = true;
				}

				pb->add_monomial(i, random_coef());
				pb->add_monomial(j, random_coef());
				pb->add_monomial(i,j, random_coef());
				if (m >= 3) {
					pb->add_monomial(k, random_coef());
					pb->add_monomial(i,k, random_coef());
					pb->add_monomial(j,k, random_coef());
					pb->add_monomial(i,j,k, random_coef());
				}
				if (m >= 4) {
					pb->add_monomial(l, random_coef());
					pb->add_monomial(i,l, random_coef());
					pb->add_monomial(j,l, random_coef());
					pb->add_monomial(k,l, random_coef());
					pb->add_monomial(i,j,l, random_coef());
					pb->add_monomial(i,k,l, random_coef());
					pb->add_monomial(j,k,l, random_coef());
					pb->add_monomial(i,j,k,l, random_coef());
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

		statusOK();
		cout << "Number of used variables : " << used.size() << endl;
	}

	if (submodular) {
		cout << "Making polynomial submodular...\n";
		pb->make_submodular();
		if (n <= 10) {
			statusTry("Checking submodularity...");
			if (!pb->is_submodular()) {
				throw runtime_error("Polynomial is not submodular");
			}
			statusOK();
		}
	}
}

