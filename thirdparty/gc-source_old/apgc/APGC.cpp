#include "StdAfx.h"
#include "apgc/APGC.h"
#include "apgc/APGCDebug.h"
#include "apgc/APGCDump.h"

APGC::APGC(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node)
:GadgetGraph(num_nodes, num_cliques, clique_size, num_cliques_per_node), _flow(0)
{
	InitializeConfigVariables();

	_joinEdge = new APGCEdge();
	_nodes = new APGCPixNode[NUM_NODES];
	_auxNodes = new APGCAuxNode[num_cliques];

#ifdef CHECK_FOR_ERRORS
	for(int j=0, counter=0; j<IMG_HEIGHT; ++j)
		for(int i=0;i<IMG_WIDTH; ++i, ++counter)
		{
			DEBUG_STATEMENT(_nodes[counter]._x = i);
			DEBUG_STATEMENT(_nodes[counter]._y = j);
			DEBUG_STATEMENT(_nodes[counter]._nodeIndex = counter);
		}
#endif
	_numAuxNodes = 0;
}

APGC::~APGC(void)
{
	delete _joinEdge;

	delete[] _nodes;
	delete[] _auxNodes;
}

void APGC::InitializeConfigVariables()
{
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
		if(num_contained_nodes != GG_CLIQUE_SIZE)
			info._containedNodes[num_contained_nodes] = -1;
	}
}

void APGC::AddUnaryTerm(int index, float e_0_b, float e_1_a)
{
	//DEBUG_STATEMENT(PRINTF("Unary %d => a=%f, b=%f\n", index, e_1_a, e_0_b));

	float min_val = MIN(e_0_b, e_1_a);
	_flow += min_val;

	e_0_b -= min_val;
	e_1_a -= min_val;

	APGCPixNode& node = _nodes[index];
	node._excess = e_1_a-e_0_b;
}

void APGC::AddHigherTerm(int* node_indexes, float* labeling_cost)
{
	APGCAuxNode& clique = _auxNodes[_numAuxNodes];

	DEBUG_STATEMENT(clique._x = _nodes[node_indexes[0]]._x);
	DEBUG_STATEMENT(clique._y = _nodes[node_indexes[0]]._y);
	clique._index = _numAuxNodes;

	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		int node_index = node_indexes[i];
		APGCPixNode& node = _nodes[node_index];

		int node_index_in_clique = GG_CLIQUE_SIZE-1-i;
		clique.InsertNode(&node, node_index_in_clique);
	}

	memcpy(clique._constraintSlacks, labeling_cost, NUM_CONSTRAINTS*sizeof(float));

	++_numAuxNodes;
}

void APGC::InitializeStrcuturesBeforeFlow()
{
	_lastPathLength = 1;
	_maxAssignedDistLabel = 0; _activeNodes.push_back(APGCNodeList()); _orphanNodes.push_back(APGCNodeList()); 

	for(int i=0; i<_numAuxNodes; ++i)
	{
		APGCAuxNode& aux_node = _auxNodes[i];
		aux_node.Reparametrize(_flow);
		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
			aux_node.RecalculateTightEdgeCodes(j);
		aux_node.RecalculateSaturatedEdgeCode();
	}

	for(int i=0; i<NUM_NODES; ++i)
	{
		APGCPixNode* pix_node = _nodes + i;

		if(pix_node->_excess < -EPSILON)
		{
			pix_node->_tree = TREE_SINK;
			pix_node->_distanceLabel = TERM_DIST;
			DeclareActive(pix_node);
			for(int j=0; j<pix_node->_numAuxNodes; ++j)
			{
				APGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[j];
				APGCAuxNode* aux_node = aux_info._auxNode;
				aux_node->_sinkLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
			}
		}
		else if(pix_node->_excess > EPSILON)
		{
			pix_node->_tree = TREE_SRC;
			pix_node->_distanceLabel = TERM_DIST;
			DeclareActive(pix_node);
			for(int j=0; j<pix_node->_numAuxNodes; ++j)
			{
				APGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[j];
				APGCAuxNode* aux_node = aux_info._auxNode;
				aux_node->_srcLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
			}
		}
	}
	//PRINTF("\n\nState after Reparametrization\n"); DumpState(); PRINTF("\n\n");

	ITER_NUM = 0;
	TIME = 0;
	_joinEdge->_auxNode = NULL;
}

double APGC::FindMaxFlow()
{
	PRINT_ENABLED = GADGET_GRAPH_DEBUG_VERBOSE;

	//PRINTF("\nInitial Flow = %.1f\n", _flow);
	//PRINTF("\n\nState before Reparametrization\n"); DumpState();

	InitializeStrcuturesBeforeFlow();

	DEBUG_STATEMENT(APGCDebug::ValidateAll(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));

	while(true)
	{
		++ITER_NUM;

		PRINTF("\n********* Iteration %d flow %.1f ********** \n", (int)ITER_NUM, _flow);

		PRINTF("\nExpanding graph\n");
		_joinEdge->_auxNode = NULL;
		ExpandGraph();
		DEBUG_STATEMENT(APGCDebug::ValidateAll(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));
		if(_joinEdge->_auxNode == NULL){
			DEBUG_STATEMENT(for(int i=0; i<NUM_NODES; ++i) MY_ASSERT(!_nodes[i]._active));
			DEBUG_STATEMENT(for(int i=0; i<_numAuxNodes; ++i) MY_ASSERT(!_auxNodes[i]._active));
			PRINTF("No path between src and sink..breaking\n\n");
			break;
		}

		float flow_augmented = AugmentFlow();
		_flow += flow_augmented;
		ProcessTightConstraints();
		ProcessOrphans();
		DEBUG_STATEMENT(APGCDebug::ValidateAll(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));
	}

	DumpState();

	return _flow;
}

void APGC::DumpState()
{
	DEBUG_STATEMENT(APGCDump::DumpState(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));
}

