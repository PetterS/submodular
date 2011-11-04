//
// Petter Strandmark 2011
// petter@maths.lth.se
//
// Complete workflow for minimizing f(x) using generalized roof duality
//


#include "PseudoBoolean.h"


namespace Petter
{

	template<typename real>
	real PseudoBoolean<real>::minimize(vector<label>& x, int& labeled, bool heuristic)
	{
		bool should_continue;
		real bound;
		labeled = 0;
		int n = int( x.size() );

		do {

			// Create symmetric relaxation
			SymmetricPseudoBoolean<real> spb;
			if (heuristic) {
				spb.create_heuristic(*this);
			}
			else {
				spb.create_lp(*this);
			}

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
