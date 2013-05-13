#include "StdAfx.h"
#include "prgc/PRGC.h"
#include "prgc/PRGCDump.h"

void PRGC::ExpandGraph()
{
	while(!_activeNodes.empty())
	{
		PRGCNode* pix_node = _activeNodes.front();
		pix_node->_active = false;
		_activeNodes.pop_front();
		if(pix_node->_distanceLabel == DIST_DEFAULT)
			continue;
		ExpandNode(pix_node);
	}
}

void PRGC::RelabelAll()
{
	DEBUG_STATEMENT(PRINTF("\nRelabeling all \n\n"));

	_activeNodes.clear();
	for(int i=0; i<_numCliques; ++i)
	{
		PRGCClique* clique = _cliques + i;
		clique->_labeledNodesCode = 0;
	}
	for(int i=0; i<NUM_NODES; ++i)
	{
		PRGCNode* node = _nodes + i;
		node->_active = false;
		if(node->_excess < -EPSILON)
			node->SetAsSink(_activeNodes);
		else
			node->_distanceLabel = DIST_DEFAULT;
	}
	ExpandGraph();
}

void PRGC::ExpandNode(PRGCNode* rec_node)
{
	DEBUG_STATEMENT(PRINTF("Processing active node "); PRGCDump::PrintNode(rec_node); PRINTLN);

	DIST_TYPE child_dist = rec_node->_distanceLabel + 1;
	for(int i=0; i<rec_node->_numCliques; ++i)
	{
		PRGCEdgeInfo& rec_node_edge_info = *(rec_node->_edgeInfo[i]);
		PRGCClique* clique = rec_node_edge_info._clique;
		int rec_node_index = rec_node_edge_info._indexInClique;

		CODE_TYPE possible_join_nodes = (clique->_labeledNodesCode ^ CODE_ALL_ONES_WORD) 
										& rec_node_edge_info._tightEdgeIntersectionCode;

		CODE_TYPE unlabeled_nodes = clique->_labeledNodesCode ^ CODE_ALL_ONES_WORD;
		CODE_TYPE nodes_i_can_give_label = unlabeled_nodes & rec_node_edge_info._tightEdgeIntersectionCode;
		while(nodes_i_can_give_label != 0)
		{
			uchar j = FIRST_NON_ZERO(nodes_i_can_give_label);
			nodes_i_can_give_label &= INDEX_CODE_ZERO[j];
			PRGCEdgeInfo& send_node_info = clique->_edgeInfo[j];
			PRGCNode* send_node = send_node_info._node;
			if(send_node->_distanceLabel != DIST_DEFAULT)
				continue;

			//DEBUG_COND(ITER_NUM == 1 && send_node->_x == 0 && send_node->_y == 0 && send_node->_nodeIndex == 0);

			for(int k=0; k<send_node->_numCliques; ++k){
				PRGCEdgeInfo& ai = *(send_node->_edgeInfo[k]);
				ai._clique->_labeledNodesCode |= INDEX_CODE_ONE[ai._indexInClique];
			}

			send_node->_distanceLabel = child_dist;
			DEBUG_STATEMENT(PRINTF("Giving distance label %d to ", child_dist); PRGCDump::PrintNode(send_node); PRINTLN); 
			if(!send_node->_active){
				send_node->_active = true;
				_activeNodes.push_back(send_node);
				DEBUG_STATEMENT(PRINTF("Pushing pix node %d,%d in active list\n", send_node->_x, send_node->_y));
			}
			if(send_node->_excess > EPSILON && !send_node->_inExcessList){
				send_node->_inExcessList = true;
				_excessNodes.push_front(send_node);
			}
		}
	}
}

void PRGC::PushFlow(PRGCNode* send_node)
{
	if(send_node->_excess <= EPSILON)
		return;
	if(send_node->_distanceLabel == DIST_DEFAULT)
		return;

	DEBUG_STATEMENT(PRINTF("Processing excess node "); PRGCDump::PrintNode(send_node); PRINTLN);

	//DEBUG_COND(ITER_NUM > 20);

	DIST_TYPE my_dist = send_node->_distanceLabel;
	int min_child_label = INT_MAX;
	for(int i=0; i<send_node->_numCliques; ++i)
	{
		PRGCEdgeInfo& send_node_edge_info = *(send_node->_edgeInfo[i]);
		PRGCClique* clique = send_node_edge_info._clique;
		int send_node_index = send_node_edge_info._indexInClique;

		int not_me = (1 << send_node_index) ^ CODE_ALL_ONES_WORD;
		CODE_TYPE possible_rec_nodes = clique->_labeledNodesCode & not_me;

		while(possible_rec_nodes != 0)
		{
			int rec_node_index = FIRST_NON_ZERO(possible_rec_nodes);
			possible_rec_nodes &= INDEX_CODE_ZERO[rec_node_index];
			PRGCEdgeInfo& rec_node_info = clique->_edgeInfo[rec_node_index];
			if((rec_node_info._tightEdgeIntersectionCode & INDEX_CODE_ONE[send_node_index]) == 0)
				continue;
			PRGCNode* rec_node = rec_node_info._node;
			if(rec_node->_distanceLabel == DIST_DEFAULT)
				continue;
			if(rec_node->_distanceLabel >= my_dist)
			{
				if(min_child_label > rec_node->_distanceLabel)
					min_child_label = rec_node->_distanceLabel;
				continue;
			}
			float cap = FindResidual(clique, send_node_index, rec_node_index);
			float flow_to_send = MIN(cap, send_node->_excess);
			DEBUG_STATEMENT(PRINTF("Pushing flow %f from ", flow_to_send); PRGCDump::PrintNode(send_node); PRINTF(" to "); PRGCDump::PrintNode(rec_node); PRINTLN);
			MY_ASSERT(flow_to_send > EPSILON);
			if(rec_node->_excess < -EPSILON)
				_flow += MIN(-rec_node->_excess, flow_to_send);
			AdjustConstraintSlacksAfterSendFlow(clique, send_node_index, rec_node_index, flow_to_send);
		
			send_node->ReduceExcess(flow_to_send);
			rec_node->IncreaseExcess(flow_to_send, _excessNodes, _activeNodes);

			if(send_node->_excess <= EPSILON)
				break;
		}
		if(send_node->_excess <= EPSILON)
			break;
	}

	if(send_node->_excess <= EPSILON)
		return;

	if(min_child_label == INT_MAX)
	{
		DEBUG_STATEMENT(PRINTF("No path to sink available for "); PRGCDump::PrintNode(send_node); PRINTF(" Relabeling to default(infinity)\n"));
		send_node->DeclareOrphan();
		return;
	}

	if(min_child_label >= PRGC_MAX_DIST_LABEL)
	{
		DEBUG_STATEMENT(PRINTF("My nearest neighbor is num_nodes away. No path to sink available for me - "); PRGCDump::PrintNode(send_node); PRINTF(" Relabeling to default(infinity)\n"));
		send_node->DeclareOrphan();
		return;
	}

	DEBUG_STATEMENT(PRINTF("Relabeling "); PRGCDump::PrintNode(send_node); PRINTF(" from %d to %d\n", send_node->_distanceLabel, min_child_label+1));
	send_node->_distanceLabel = min_child_label + 1;
	send_node->_inExcessList = true;
	_excessNodes.push_front(send_node);
}

float PRGC::FindResidual(PRGCClique* clique, int send_node_index, int rec_node_index)
{
	int* constraints_list = CONTAINING_CONSTRAINTS[rec_node_index];
	CODE_TYPE send_node_code = INDEX_CODE_ONE[send_node_index];

	float min_slack = PRGC_CAP_INFINITY;

	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints_list[i];
		CODE_TYPE constraint_labeling = CONSTRAINT_LABELING(constraint_index);
		if((constraint_labeling & send_node_code) != 0)
			continue;
		float slack = clique->_constraintSlacks[constraint_index];
		if(min_slack > slack)
			min_slack = slack;
	}

	return min_slack;
}

void PRGC::AdjustConstraintSlacksAfterSendFlow(PRGCClique* clique, int send_node_index, int rec_node_index, float flow)
{
	int* constraints_list; 

	constraints_list = CONTAINING_CONSTRAINTS[rec_node_index];
	CODE_TYPE send_node_code = INDEX_CODE_ONE[send_node_index];

	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints_list[i];
		CODE_TYPE constraint_labeling = CONSTRAINT_LABELING(constraint_index);
		if((constraint_labeling & send_node_code) != 0)
			continue;
		float& slack = clique->_constraintSlacks[constraint_index];
		MY_ASSERT(slack > EPSILON);
		slack -= flow;
		if(slack > EPSILON)
			continue;
		clique->DeclareConstraintTight(constraint_index, constraint_labeling);
	}

	if(send_node_index == -1)
		return;

	constraints_list = CONTAINING_CONSTRAINTS[send_node_index];
	CODE_TYPE rec_node_code = INDEX_CODE_ONE[rec_node_index];

	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints_list[i];
		CODE_TYPE constraint_labeling = CONSTRAINT_LABELING(constraint_index);
		if((constraint_labeling & rec_node_code) != 0)
			continue;
		float& slack = clique->_constraintSlacks[constraint_index];
		bool previously_tight = (slack <= EPSILON);
		slack += flow;
		if(!previously_tight || slack <= EPSILON)
			continue;
		clique->DeclareConstraintNotTight(constraint_index, constraint_labeling, _activeNodes);
	}
}

