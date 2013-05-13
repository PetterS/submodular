#include "StdAfx.h"
#include "bkgc/BKGC.h"
#include "bkgc/BKGCSinkOrphanRootProcessor.h"
#include "bkgc/BKGCSrcOrphanRootProcessor.h"
#include "test/CliquePotential.h"
#include "bkgc/BKGCDataStructures.h"
#include "bkgc/BKGCDebug.h"
#include "bkgc/BKGCDump.h"

BKGC::BKGC(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node)
:GadgetGraph(num_nodes, num_cliques, clique_size, num_cliques_per_node)
{
	InitializeConfigVariables();
	EPSILON = 0.1f;

	_joinEdge = new BKGCEdge();
	_nodes = new BKGCPixNode[num_nodes];
	_auxNodes = new BKGCAuxNode[num_cliques];

#ifdef CHECK_FOR_ERRORS
	for(int j=0, counter=0; j<IMG_HEIGHT; ++j)
		for(int i=0;i<IMG_WIDTH; ++i, ++counter){
			DEBUG_STATEMENT(_nodes[counter]._x = i);
			DEBUG_STATEMENT(_nodes[counter]._y = j);
		}
#endif

	_activeNodes = new std::list<BKGCNode*>();
	_sinkTreeOrphanRootProcessor = new BKGCSinkOrphanRootProcessor();
	_srcTreeOrphanRootProcessor = new BKGCSrcOrphanRootProcessor();
	_numAuxNodes = 0;
	_flow = 0;

	BKGC_TIGHT_CONSTRAINTS = &_tightConstraints;
	BKGC_ACTIVE_NODES = _activeNodes;
}

BKGC::~BKGC(void)
{
	delete _joinEdge;
	delete _sinkTreeOrphanRootProcessor;
	delete _srcTreeOrphanRootProcessor;
	delete _activeNodes;

	delete[] _nodes;
	delete[] _auxNodes;
}

void BKGC::InitializeConfigVariables()
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

int BKGC::GetLabel(int node_id){
	if(_nodes[node_id]._tree == TREE_SRC)
		return LABEL_SOURCE_0;
	else 
		return LABEL_SINK_1;
}

void BKGC::AddHigherTerm(int* node_indexes, float* labeling_cost)
{
	BKGCAuxNode& clique = _auxNodes[_numAuxNodes];

	DEBUG_STATEMENT(clique._x = _nodes[node_indexes[0]]._x);
	DEBUG_STATEMENT(clique._y = _nodes[node_indexes[0]]._y);
	
	clique._index = _numAuxNodes;

	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		int node_index = node_indexes[i];
		BKGCPixNode& node = _nodes[node_index];

		int node_index_in_clique = GG_CLIQUE_SIZE-1-i;
		clique.InsertNode(&node, node_index_in_clique);
	}

	memcpy(clique._constraintSlacks, labeling_cost, NUM_CONSTRAINTS*sizeof(float));

	++_numAuxNodes;
}

void BKGC::AddUnaryTerm(int index, float e_0_b, float e_1_a)
{
	BKGCPixNode& node = _nodes[index];
	_flow += MIN(e_1_a , e_0_b);
	node._excess = e_1_a-e_0_b;
}

void BKGC::InitializeStucturesBeforeFlow()
{
	for(int i=0; i<_numAuxNodes; ++i)
	{
		BKGCAuxNode& aux_node = _auxNodes[i];
		aux_node.Reparametrize(_flow);
		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
			aux_node.RecalculateTightEdgeCodes(j);
		aux_node.RecalculateSaturatedEdgeCode();
	}

	for(int i=0; i<NUM_NODES; ++i)
	{
		BKGCPixNode& node = _nodes[i];

		if(EQUAL(node._excess,0)){
			node._tree = TREE_NONE;
			continue;
		}
		_activeNodes->push_back(&node);
		node._active = true;
		if(LESS_THAN(node._excess,0)){
			node._tree = TREE_SINK;
			for(int i=0; i<node._numAuxNodes; ++i){
				BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[i];
				aux_info._auxNode->_sinkLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
			}
		}else{
			node._tree = TREE_SRC;
			for(int i=0; i<node._numAuxNodes; ++i){
				BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[i];
				aux_info._auxNode->_srcLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
			}
		}
	}

	//PRINTF("\n\nState after Reparametrization\n"); DumpState(); PRINTF("\n\n");

	ITER_NUM = 0;
	TIME = 0;
	_joinEdge->_auxNode = NULL;
}

double BKGC::FindMaxFlow()
{
	bool print_enabled_backup = PRINT_ENABLED;
	PRINT_ENABLED = GADGET_GRAPH_DEBUG_VERBOSE;

	InitializeStucturesBeforeFlow();

	DEBUG_STATEMENT(BKGCDebug::ValidateAll(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));

	while(true)
	{
		++ITER_NUM;
		DEBUG_STATEMENT(PRINTF("\n********* Iteration %d, %f ********** \n", (int)ITER_NUM, _flow));

		if(_joinEdge->_auxNode != NULL)
		{
			if(_joinEdge->_sendingNodeIndex == -1)
			{
				BKGCPixNode* sink_side_node = _joinEdge->_auxNode->_nodeInfo[_joinEdge->_receivingNodeIndex]._node;
				if(_joinEdge->_auxNode->_parentIndex == -1 || (sink_side_node->_excess >= -EPSILON && sink_side_node->_parentAuxNodeIndex == -1))
					_joinEdge->_auxNode = NULL;
			}else
			{
				BKGCPixNode* src_side_node = _joinEdge->_auxNode->_nodeInfo[_joinEdge->_sendingNodeIndex]._node;
				if(_joinEdge->_receivingNodeIndex == -1)
				{
					if(_joinEdge->_auxNode->_parentIndex == -1 || (src_side_node->_excess <= EPSILON && src_side_node->_parentAuxNodeIndex == -1))
						_joinEdge->_auxNode = NULL;		
				}else
				{
					BKGCPixNode* sink_side_node = _joinEdge->_auxNode->_nodeInfo[_joinEdge->_receivingNodeIndex]._node;
					if((sink_side_node->_excess >= -EPSILON && sink_side_node->_parentAuxNodeIndex == -1) || (src_side_node->_excess <= EPSILON && src_side_node->_parentAuxNodeIndex == -1))
						_joinEdge->_auxNode = NULL;		
				}
			}
		}

		if(_joinEdge->_auxNode == NULL)
		{
			DEBUG_STATEMENT( PRINTF("expanding graph\n" ));
			ProcessActives();
			DEBUG_STATEMENT(BKGCDebug::ValidateAll(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));
			if(_joinEdge->_auxNode == NULL){
				DEBUG_STATEMENT(PRINTF("No path between src and sink..breaking\n\n"));
				break;
			}
		}
		float flow_augmented = AugmentFlow();
		//if(flow_augmented < EPSILON/GG_CLIQUE_SIZE){
		//	printf("less flow augmentation %f in iter %d\n", flow_augmented, ITER_NUM);
		//	exit(1);
		//}
		_flow += flow_augmented;
		ProcessTightConstraints();
		DEBUG_STATEMENT(BKGCDebug::ValidateAll(_nodes, NUM_NODES, _auxNodes, _numAuxNodes));
	}

	PRINT_ENABLED = print_enabled_backup;

	return _flow;
}

#ifdef CHECK_FOR_ERRORS
void BKGC::DumpState() 
{ 
	BKGCDump::DumpState(_nodes, NUM_NODES, _auxNodes, _numAuxNodes);
}
#endif
