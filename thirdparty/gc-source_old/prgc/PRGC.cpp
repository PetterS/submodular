#include "StdAfx.h"
#include "prgc/PRGC.h"
#include "prgc/PRGCDebug.h"
#include "prgc/PRGCDump.h"

PRGC::PRGC(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node)
:GadgetGraph(num_nodes, num_cliques, clique_size, num_cliques_per_node), _flow(0)
{
	InitializeConfigVariables();

	_nodes = new PRGCNode[NUM_NODES];

	_cliques = new PRGCClique[num_cliques];
	_numCliques = 0;

#ifdef CHECK_FOR_ERRORS
	for(int j=0, counter=0; j<IMG_HEIGHT; ++j)
		for(int i=0;i<IMG_WIDTH; ++i, ++counter)
		{
			DEBUG_STATEMENT(_nodes[counter]._x = i);
			DEBUG_STATEMENT(_nodes[counter]._y = j);
		}
#endif
}

PRGC::~PRGC(void)
{
	delete[] _nodes;
	delete[] _cliques;
}

void PRGC::InitializeConfigVariables()
{
	PRGC_MAX_DIST_LABEL = NUM_NODES;

	std::vector<int> num_cont_constraints_vec, num_non_cont_constraints_vec;
	for(int i=0; i<GG_CLIQUE_SIZE; ++i){
		num_cont_constraints_vec.push_back(0);
		num_non_cont_constraints_vec.push_back(0);
	}
	for(int i=0; i<NUM_CONSTRAINTS; ++i)
	{
		ConstraintInfo& info = CONSTRAINT_INFO[i];
		
		info._labeling = i;

		int num_contained_nodes = 0;

		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
		{
			if(IS_BIT_0(info._labeling, j))
			{
				int& num_non_cont_constraints = num_non_cont_constraints_vec[j];
				NON_CONTAINING_CONSTRAINTS[j][num_non_cont_constraints++] = i;
			}else
			{
				int& num_cont_constraints = num_cont_constraints_vec[j];
				CONTAINING_CONSTRAINTS[j][num_cont_constraints++] = i;

				info._containedNodes[num_contained_nodes++] = j;
			}
		}
	}
}

void PRGC::AddUnaryTerm(int index, float e_0_b, float e_1_a)
{
	//PRINTF("Unary %d => a=%f, b=%f\n", index, e_1_a, e_0_b);

	float min_val = MIN(e_0_b, e_1_a);
	_flow += min_val;

	e_0_b -= min_val;
	e_1_a -= min_val;

	PRGCNode& node = _nodes[index];
	node._excess = e_1_a-e_0_b;
}

void PRGC::AddHigherTerm(int* node_indexes, float* labeling_cost)
{
	PRGCClique& clique = _cliques[_numCliques];

	DEBUG_STATEMENT(clique._x = _nodes[node_indexes[0]]._x);
	DEBUG_STATEMENT(clique._y = _nodes[node_indexes[0]]._y);
	DEBUG_STATEMENT(clique._index = _numCliques);

	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		int node_index = node_indexes[i];
		PRGCNode& node = _nodes[node_index];

		int index = GG_CLIQUE_SIZE-1-i;
		clique.InsertNode(&node, index);
	}

	memcpy(clique._constraintSlacks, labeling_cost, NUM_CONSTRAINTS*sizeof(float));

	++_numCliques;
}

double PRGC::FindMaxFlow()
{
	DEBUG_STATEMENT(PRINTF("\nInitial flow = %f\n\n", _flow));

	//PRINT_ENABLED = true;
	//DEBUG_STATEMENT(PRINTF("\n\nState before Reparametrization\n"); DumpState());

	for(int i=0; i<_numCliques; ++i)
	{
		PRGCClique& clique = _cliques[i];
		clique.Reparametrize(_excessNodes, _activeNodes, _flow);
		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
			clique.RecalculateTightEdgeCodes(j);
		clique.RecalculateSaturatedEdgeCode();
		//DEBUG_STATEMENT(PRINTF("Pushing aux node %d,%d,%d in active list\n", aux->_x, clique->_y, clique->_index));
	}

	for(int i=0; i<NUM_NODES; ++i)
	{
		PRGCNode& node = _nodes[i];
		if(node._excess > EPSILON)
		{
			node._inExcessList = true;
			_excessNodes.push_back(&node);
		}		
	}

	//DEBUG_STATEMENT(PRINTF("\n\nState after Reparametrization\n"); DumpState(); PRINTF("\n\n"));

	//_flow = 0;
	ITER_NUM = 0;

	//DEBUG_STATEMENT(GraphDebug::ValidateAll(_nodes, NUM_NODES, _cliques, _numCliques));
	//DEBUG_STATEMENT(if(DEBUG_START_ITER >= 0) PRINT_ENABLED = false);

	while(true)
	{
		//if((ITER_NUM % 1000) == 0)
		//	printf("[%d,%f] ", ITER_NUM,_flow);

		++ITER_NUM;

		//DEBUG_STATEMENT(if(ITER_NUM >= DEBUG_START_ITER) PRINT_ENABLED = true);
		DEBUG_STATEMENT(PRINTF("\n********* Iteration %I64d flow = %f ********** \n", ITER_NUM, (float)_flow));

		if((ITER_NUM % PRGC_RELABEL_ALL_AFTER_EVERY) == 1){
			RelabelAll();
			DEBUG_STATEMENT(PRGCDebug::ValidateAll(_nodes, NUM_NODES, _cliques, _numCliques));
		}

		if(!_activeNodes.empty()){
			DEBUG_STATEMENT(PRINTF("expanding graph\n"));
			ExpandGraph();
			DEBUG_STATEMENT(PRGCDebug::ValidateAll(_nodes, NUM_NODES, _cliques, _numCliques));
		}

		if(_excessNodes.empty())
		{
			DEBUG_STATEMENT(PRINTF("No excess node available, flow maximal, breaking...\n"));
			DEBUG_STATEMENT(PRGCDebug::ValidateAll(_nodes, NUM_NODES, _cliques, _numCliques));
			break;
		}

		PRGCNode* node = _excessNodes.front();
		_excessNodes.pop_front();
		node->_inExcessList = false;
		PushFlow(node);

		DEBUG_STATEMENT(PRGCDebug::ValidateAll(_nodes, NUM_NODES, _cliques, _numCliques));
	}

	DEBUG_STATEMENT(PRINTF("Doing relabel for one last time, for setting sink and src reachable nodes\n"));
	RelabelAll();
	DEBUG_STATEMENT(PRINTF("\n\n"));

	//DEBUG_STATEMENT(PRINTF("\n\nLast State \n"); DumpState());

	return _flow;
}

void PRGC::DumpState()
{
	DEBUG_STATEMENT(PRGCDump::DumpState(_nodes, NUM_NODES, _cliques, _numCliques));
}

