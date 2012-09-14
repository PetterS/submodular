//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Complete workflow for minimizing f(x) using generalized roof duality
//


#include "PseudoBoolean.h"

#include <cstdlib>

namespace Petter
{

	template<typename real>
	real PseudoBoolean<real>::minimize(vector<label>& x, Method method, const char* generators_file)
	{
		int tmp;
		return minimize(x,tmp,method,generators_file);
	}

	template<typename real>
	real PseudoBoolean<real>::minimize(vector<label>& x, int& labeled, Method method, const char* generators_file)
	{
		if (method==HOCR) {
			return minimize_reduction(x,labeled);
		}
		else if (method==Fix) {
			return minimize_reduction_fixetal(x,labeled);
		}
		else if (method==LP) {
			labeled=0;
			return minimize_lp(x);
		}
		else if (method==M_reduction) {
			return minimize_reduction_M(x,labeled);
		}

		bool should_continue;
		real bound = -1;
		labeled = 0;
		int n = int( x.size() );

		Generators<real>* generators = 0;
		if (method == GRD_gen) {
			if (generators_file) {
				generators = new Generators<real>(generators_file);
			}
			else if (getenv("GENERATORS")) {
				generators = new Generators<real>(std::getenv("GENERATORS"));
			}
			else {
				generators = new Generators<real>("generators/generators.txt");
			}
		}

		do {

			int new_labeled = 0;

			if (method==GRD || method==GRD_heur) {
				// Create symmetric relaxation
				SymmetricPseudoBoolean<real> spb;
				if (method==GRD_heur) {
					spb.create_heuristic(*this);
				}
				else {
					spb.create_lp(*this);
				}

				// Minimize relaxation
				bound = spb.minimize(x, new_labeled);
			}
			else if (method == GRD_gen) {
				// Create symmetric relaxation
				GeneratorPseudoBoolean<real> spb(*generators);
				spb.create_lp(*this);
				// Minimize relaxation
				bound = spb.minimize(x, new_labeled);
			}

			// If we have more persistencies, continue
			should_continue = new_labeled > labeled;
			labeled = new_labeled;
			if (labeled == n) {
				//Nothing more to do
				should_continue = false;
			}

			// Reduce this function
			reduce(x);

		} while (should_continue);


		if (generators) {
			delete generators;
		}

		return bound;

	}



	template<typename real>
	real PseudoBoolean<real>::minimize(vector<label>& x, int& labeled, bool heuristic)
	{
		bool should_continue;
		real bound;
		labeled = 0;
		int n = int( x.size() );


		do {

			
			if (heuristic) {
				// Create symmetric relaxation
				SymmetricPseudoBoolean<real> spb;

				spb.create_heuristic(*this);

				// Minimize relaxation
				int new_labeled = 0;
				bound = spb.minimize(x, new_labeled);

				// If we have more persistencies, continue
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				// Reduce this function
				reduce(x);
			}
			else {
				// Create symmetric relaxation, double needed
				SymmetricPseudoBoolean<double> spb;

				spb.create_lp(*this);

				// Minimize relaxation
				int new_labeled = 0;
				bound = spb.minimize(x, new_labeled);

				// If we have more persistencies, continue
				should_continue = new_labeled > labeled;
				labeled = new_labeled;
				if (labeled == n) {
					//Nothing more to do
					should_continue = false;
				}

				// Reduce this function
				reduce(x);
			}

			

		} while (should_continue);

		return bound;
	}


	template<typename real>
	real PseudoBoolean<real>::minimize_generators(vector<label>& x, int& labeled, bool heuristic)
	{
		bool should_continue;
		real bound;
		labeled = 0;
		int n = int( x.size() );

		ASSERT(!heuristic);

		// TODO: should be provided by user
		Generators<real> generators("generators/generators.txt");

		do {
			// Create symmetric relaxation
			GeneratorPseudoBoolean<real> spb(generators);
			spb.create_lp(*this);

			// Minimize relaxation
			int new_labeled = 0;
			bound = spb.minimize(x, new_labeled);

			// If we have more persistencies, continue
			should_continue = new_labeled > labeled;
			labeled = new_labeled;
			if (labeled == n) {
				//Nothing more to do
				should_continue = false;
			}

			// Reduce this function
			reduce(x);

		} while (should_continue);

		return bound;
	}
}

#include "pb_instances.inc"
