real GeneratorPseudoBoolean<real>::minimize_version_2(vector<label>& x, int& nlabelled) const
	{
		//ASSERT_STR(this->gen.ngen4 == 0, "Degree-4 generators are not yet supported.");

		index nVars = index( x.size() ); // Number of variables.
		int n = 2 * nVars;  // Because we have x and y.
		int num_cliques = 0;

		for (auto itr = alphaij.begin(); itr != alphaij.end(); ++itr) {
			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen2;++ii) {
				real alpha = vec.at(ii);
				if (alpha > 0) {
					num_cliques++;
				}
			}
		}

		for (auto itr = alphaijk.begin(); itr != alphaijk.end(); ++itr) {
			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen3;++ii) {
				real alpha = vec.at(ii);
				if (alpha > 0) {
					// Add monomials for this generator to the graph
					num_cliques++;
				}
			}
		}

		for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
			const auto& vec = itr->second;
			for (int ii = 0; ii < std::max(gen.ngen4pos, gen.ngen4neg); ++ii) {
				real alpha = vec.at(ii);
				if (alpha > 0) {
					// Add monomials for this generator to the graph
					num_cliques++;
				}
			}
		}

		
		real C = 0; // Constant in objective function.
		int clique_size = 4;
		int num_cliques_per_node = 2 * num_cliques; // TODO: Fix this. (Is this parameter used by GC?)

		// We add two extra variables in order to be able to add degree-2 cliques
		// as degree-4 cliques.

		typedef PRGC GCType;
		//typedef APGC GCType;

		GCType graph(n + 2,
		             2 * num_cliques, // Each generator gives two cliques.
		             clique_size,
		             num_cliques_per_node);
		int extra1 = n;
		int extra2 = n + 1;

		std::unique_ptr<PseudoBoolean<real>> f_debug;
		
		// Uncomment this line to also minimize g exhaustively.
		//f_debug.reset(new PseudoBoolean<real>);

		//
		// Degree-1 terms.
		//
		for (auto itr = alphai.begin(); itr != alphai.end(); ++itr) {
			int i = itr->first;
			real alpha = itr->second;

			graph.AddUnaryTerm(i,         0,     alpha);
			graph.AddUnaryTerm(i + nVars, alpha,     0);

			if (f_debug) {
				f_debug->add_clique(i, 0, alpha);
				f_debug->add_clique(i + nVars, alpha, 0);
			}
		}

		//
		// Go through all alphas which correspond to quadratic generators
		//
		for (auto itr = alphaij.begin(); itr != alphaij.end(); ++itr) {
			const pair& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);

			vector<int> idx(4); // Translates from "local" indices to "global"
			idx[0] = i; // x variables
			idx[1] = j;
			idx[2] = i + nVars; // y variables
			idx[3] = j + nVars;

			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen2;++ii) {
				real alpha = vec.at(ii);
				auto& generator = gen.gen2.at(ii);
				if (alpha > 0) {
					// Add cliques for this generator to the graph
					{
						float E1[]= {alpha * generator.values1.at(0), // E0000
						             alpha * generator.values1.at(1), // E0001
						             alpha * generator.values1.at(2), // E0010
						             alpha * generator.values1.at(3), // E0011
						             alpha * generator.values1.at(0), // E0100
						             alpha * generator.values1.at(1), // E0101
						             alpha * generator.values1.at(2), // E0110
						             alpha * generator.values1.at(3), // E0111
						             alpha * generator.values1.at(0), // E1000
						             alpha * generator.values1.at(1), // E1001
						             alpha * generator.values1.at(2), // E1010
						             alpha * generator.values1.at(3), // E1011
						             alpha * generator.values1.at(0), // E1100
						             alpha * generator.values1.at(1), // E1101
						             alpha * generator.values1.at(2), // E1110
						             alpha * generator.values1.at(3)};// E1111
						int indices1[] = {extra1,
						                  extra2,
						                  idx.at(generator.indices1.at(0)),
						                  idx.at(generator.indices1.at(1))};
						C += make_clique_positive(clique_size, E1);
						graph.AddHigherTerm(indices1, E1);

						if (f_debug) {
							std::vector<real> Ev(E1, E1+16);
							f_debug->add_clique(indices1[0], indices1[1], indices1[2], indices1[3], Ev);
						}
					}

					{
						float E2[]= {alpha * generator.values2.at(0), // E0000
						             alpha * generator.values2.at(1), // E0001
						             alpha * generator.values2.at(2), // E0010
						             alpha * generator.values2.at(3), // E0011
						             alpha * generator.values2.at(0), // E0100
						             alpha * generator.values2.at(1), // E0101
						             alpha * generator.values2.at(2), // E0110
						             alpha * generator.values2.at(3), // E0111
						             alpha * generator.values2.at(0), // E1000
						             alpha * generator.values2.at(1), // E1001
						             alpha * generator.values2.at(2), // E1010
						             alpha * generator.values2.at(3), // E1011
						             alpha * generator.values2.at(0), // E1100
						             alpha * generator.values2.at(1), // E1101
						             alpha * generator.values2.at(2), // E1110
						             alpha * generator.values2.at(3)};// E1111
						int indices2[] = {extra1,
						                  extra2,
						                  idx.at(generator.indices2.at(0)),
						                  idx.at(generator.indices2.at(1))};
						C += make_clique_positive(clique_size, E2);
						graph.AddHigherTerm(indices2, E2);

						if (f_debug) {
							std::vector<real> Ev(E2, E2+16);
							f_debug->add_clique(indices2[0], indices2[1], indices2[2], indices2[3], Ev);
						}
					}
				}
			}
		}

		//
		// Go through all alphas which correspond to cubic generators
		//
		for (auto itr = alphaijk.begin(); itr != alphaijk.end(); ++itr) {
			const triple& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);

			vector<int> idx(6); // Translates from "local" indices to "global"
			idx.at(0) = i; // x variables
			idx.at(1) = j;
			idx.at(2) = k;
			idx.at(3) = i + nVars; // y variables
			idx.at(4) = j + nVars;
			idx.at(5) = k + nVars;

			const auto& vec = itr->second;
			for (int ii=0;ii<gen.ngen3;++ii) {
				real alpha = vec.at(ii);
				auto& generator = gen.gen3.at(ii);
				if (alpha > 0) {
					// Add cliques for this generator to the graph

					{
						float E1[]= {alpha * generator.values1.at(0), // E0000
						             alpha * generator.values1.at(1), // E0001
						             alpha * generator.values1.at(2), // E0010
						             alpha * generator.values1.at(3), // E0011
						             alpha * generator.values1.at(4), // E0100
						             alpha * generator.values1.at(5), // E0101
						             alpha * generator.values1.at(6), // E0110
						             alpha * generator.values1.at(7), // E0111
						             alpha * generator.values1.at(0), // E1000
						             alpha * generator.values1.at(1), // E1001
						             alpha * generator.values1.at(2), // E1010
						             alpha * generator.values1.at(3), // E1011
						             alpha * generator.values1.at(4), // E1100
						             alpha * generator.values1.at(5), // E1101
						             alpha * generator.values1.at(6), // E1110
						             alpha * generator.values1.at(7)};// E1111
						int indices1[] = {extra1,
						                  idx.at(generator.indices1.at(0)),
										  idx.at(generator.indices1.at(1)),
										  idx.at(generator.indices1.at(2))};
						C += make_clique_positive(clique_size, E1);
						graph.AddHigherTerm(indices1, E1);

						if (f_debug) {
							std::vector<real> Ev(E1, E1+16);
							f_debug->add_clique(indices1[0], indices1[1], indices1[2], indices1[3], Ev);
						}
					}

					{
						float E2[]= {alpha * generator.values2.at(0), // E0000
						             alpha * generator.values2.at(1), // E0001
						             alpha * generator.values2.at(2), // E0010
						             alpha * generator.values2.at(3), // E0011
						             alpha * generator.values2.at(4), // E0100
						             alpha * generator.values2.at(5), // E0101
						             alpha * generator.values2.at(6), // E0110
						             alpha * generator.values2.at(7), // E0111
						             alpha * generator.values2.at(0), // E1000
						             alpha * generator.values2.at(1), // E1001
						             alpha * generator.values2.at(2), // E1010
						             alpha * generator.values2.at(3), // E1011
						             alpha * generator.values2.at(4), // E1100
						             alpha * generator.values2.at(5), // E1101
						             alpha * generator.values2.at(6), // E1110
						             alpha * generator.values2.at(7)};// E1111
						int indices2[] = {extra1,
						                  idx.at(generator.indices2.at(0)),
						                  idx.at(generator.indices2.at(1)),
						                  idx.at(generator.indices2.at(2))};
						C += make_clique_positive(clique_size, E2);
						graph.AddHigherTerm(indices2, E2);

						if (f_debug) {
							std::vector<real> Ev(E2, E2+16);
							f_debug->add_clique(indices2[0], indices2[1], indices2[2], indices2[3], Ev);
						}
					}
				}
			}
		}

		//
		// Go through all alphas which correspond to quartic generators
		//

		for (auto itr = alphaijkl.begin(); itr != alphaijkl.end(); ++itr) {
			const quad& ind = itr->first;
			int i=get_i(ind);
			int j=get_j(ind);
			int k=get_k(ind);
			int l=get_l(ind);
			const auto& vec = itr->second;

			vector<int> idx(8); // Translates from "local" indices to "global"
			idx.at(0) = i; // x variables
			idx.at(1) = j;
			idx.at(2) = k;
			idx.at(3) = l;
			idx.at(4) = i + nVars; // y variables
			idx.at(5) = j + nVars;
			idx.at(6) = k + nVars;
			idx.at(7) = l + nVars;

			// Was a positive or negative generator used?
			bool pos = false;
			auto mitr = posgen4.find(ind);
			if (mitr != posgen4.end()) {
				pos = mitr->second;
			}

			ASSERT(gen.ngen4pos == gen.ngen4neg);
			// The code below assumes a clique size of 4.
			ASSERT(clique_size == 4);

			for (int ii=0;ii<gen.ngen4pos;++ii) {
				real alpha = vec.at(ii);

				if (alpha > 0) {

					// Add monomials for this generator to the graph
					if (pos) {
						// Positive generator was used
						auto& generator = gen.gen4pos.at(ii);
						int ii, jj, kk, ll;

						ii = idx.at(generator.indices1.at(0));
						jj = idx.at(generator.indices1.at(1));
						kk = idx.at(generator.indices1.at(2));
						ll = idx.at(generator.indices1.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values1, alpha, f_debug);

						ii = idx.at(generator.indices2.at(0));
						jj = idx.at(generator.indices2.at(1));
						kk = idx.at(generator.indices2.at(2));
						ll = idx.at(generator.indices2.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values2, alpha, f_debug);
					}
					else {
						// Negative generator was used
						auto& generator = gen.gen4neg.at(ii);
						int ii, jj, kk, ll;

						ii = idx.at(generator.indices1.at(0));
						jj = idx.at(generator.indices1.at(1));
						kk = idx.at(generator.indices1.at(2));
						ll = idx.at(generator.indices1.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values1, alpha, f_debug);

						ii = idx.at(generator.indices2.at(0));
						jj = idx.at(generator.indices2.at(1));
						kk = idx.at(generator.indices2.at(2));
						ll = idx.at(generator.indices2.at(3));
						add_generator_to_graph(&graph, &C, ii, jj ,kk ,ll, generator.values2, alpha, f_debug);
					}
				}
			}
		}


		double min_g = constant + C + graph.FindMaxFlow();
		vector<label> xfull(n);
		for (int i = 0; i < n; ++i) {
			xfull[i] = graph.GetLabel(i);
		}

		if (f_debug) {
			std::cout << "Generic cuts\n";
			std::cout << "C=" << C << " min_g=" << min_g << "\n";

			for (int i = 0; i < nVars; ++i) {
				std::cout << xfull[i];
			}
			std::cout << ", ";
			for (int i = 0; i < nVars; ++i) {
				std::cout << xfull[i + nVars];
			}
			std::cout << "\n";
		}

		nlabelled = 0;
		for (int i=0; i<nVars; ++i) {
			bool used = false;
			auto itr = var_used.find(i);
			if (itr != var_used.end()) {
				used = itr->second;
			}

			if (used) {
			    x[i]     = xfull[i];
				label yi = xfull[i+nVars];
				if (x[i] == yi) {
					x[i] = -1;
				}
				else {
					nlabelled++;
				}
			}
			else {
				// This variable is not part of the polynomial,
				// therefore labelled
				if (x[i]<0) {
					x[i]=0;
				}
				nlabelled++;
			}
		}


		if (f_debug) {
			//
			// Minimize f_debug with exhaustive search.
			//
			vector<label> x_debug(n + 2, 0), x_debug_opt(n + 2, 0);
			real optimum = f_debug->eval(x_debug);
			while (true) {
				x_debug[0]++;
				int i=0;
				while (x_debug[i]>1) {
					x_debug[i]=0;
					i++;
					if (i == n + 2) {
						break;
					}
					x_debug[i]+=1;
				}
				if (i == n + 2) {
					break;
				}

				real energy = f_debug->eval(x_debug);
				if (energy < optimum) {
					optimum = energy;
					x_debug_opt = x_debug;
				}
			}

			std::cout << "Exhaustive debug\n";
			std::cout << "C=" << C << " min_f_debug=" << constant + C + optimum << "\n";
			for (int i = 0; i < nVars; ++i) {
				std::cout << x_debug_opt[i];
			}
			std::cout << ", ";
			for (int i = 0; i < nVars; ++i) {
				std::cout << x_debug_opt[i + nVars];
			}
			std::cout << "\n";
		}

		return min_g;
	}