#include "stdafx.h"
#include "GadgetGraph.h"

GadgetGraph::GadgetGraph(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node)
{
	DEBUG_STATEMENT(GADGET_GRAPH_DEBUG_VERBOSE = true);
	DEBUG_STATEMENT(PRINT_ENABLED = true);

	NUM_NODES = num_nodes;
	GG_CLIQUE_SIZE = clique_size;
	D_GG_CLIQUE_SIZE = 2*GG_CLIQUE_SIZE;
	NUM_CLIQUES_PER_NODE = num_cliques_per_node;
	
	ITER_NUM = 0;
	TIME = 0;

	CODE_ALL_ONES_WORD = (CODE_TYPE)((1 << GG_CLIQUE_SIZE)-1);
	CODE_ALL_ONES_DWORD = (CODE_TYPE)((1 << D_GG_CLIQUE_SIZE)-1);
	INDEX_CODE_ONE = new CODE_TYPE[D_GG_CLIQUE_SIZE];
	INDEX_CODE_ZERO = new CODE_TYPE[D_GG_CLIQUE_SIZE];
	for(int i=0; i<D_GG_CLIQUE_SIZE; ++i)
	{
		CODE_TYPE bit_code = (1 << i);
		INDEX_CODE_ONE[i] = bit_code;
		INDEX_CODE_ZERO[i] = (bit_code ^ CODE_ALL_ONES_DWORD);
	}

	NUM_CONSTRAINTS = (1 << GG_CLIQUE_SIZE);
	CONSTRAINT_INFO = new ConstraintInfo[NUM_CONSTRAINTS];

	NUM_CONTAINING_CONSTRAINTS = NUM_CONSTRAINTS/2;

	CONTAINING_CONSTRAINTS = new int*[D_GG_CLIQUE_SIZE];
	for(int i=0; i<D_GG_CLIQUE_SIZE; ++i)
		CONTAINING_CONSTRAINTS[i] = new int[NUM_CONTAINING_CONSTRAINTS];

	NON_CONTAINING_CONSTRAINTS = new int*[D_GG_CLIQUE_SIZE];
	for(int i=0; i<D_GG_CLIQUE_SIZE; ++i)
		NON_CONTAINING_CONSTRAINTS[i] = new int[NUM_CONTAINING_CONSTRAINTS];
}

GadgetGraph::~GadgetGraph()
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
		delete[] CONTAINING_CONSTRAINTS[i];
	delete[] CONTAINING_CONSTRAINTS;

	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
		delete[] NON_CONTAINING_CONSTRAINTS[i];
	delete[] NON_CONTAINING_CONSTRAINTS;

	delete[] CONSTRAINT_INFO;

	delete[] INDEX_CODE_ONE;
	delete[] INDEX_CODE_ZERO;
}