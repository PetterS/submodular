#include "StdAfx.h"
#include "prgc/PRGCDataStructures.h"
#include "prgc/PRGCDump.h"

void PRGCNode::ReduceExcess(float flow_to_send)
{
	_excess -= flow_to_send;
}

void PRGCNode::IncreaseExcess(float flow_to_send, PRGCPNodeList& excess_list, PRGCPNodeList& active_list)
{
	_excess += flow_to_send;

	if(_distanceLabel == TERM_DIST && _excess > -EPSILON){
		//sink saturated here, do something to find new distance label for it and then break
		FindNewDistanceLabel();
	}

	if(!_inExcessList && _excess > EPSILON){
		_inExcessList = true;
		excess_list.push_back(this);
	}
}

void PRGCNode::FindNewDistanceLabel()
{
	int min_child_label = INT_MAX;
	for(int i=0; i<_numCliques; ++i)
	{
		PRGCEdgeInfo& send_node_edge_info = *(_edgeInfo[i]);
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
			if(min_child_label > rec_node->_distanceLabel)
				min_child_label = rec_node->_distanceLabel;
		}
	}

	if(min_child_label == INT_MAX)
	{
		DEBUG_STATEMENT(PRINTF("No path to sink available for "); PRGCDump::PrintNode(this); PRINTF(" Relabeling to default(infinity)\n"));
		DeclareOrphan();
		return;
	}

	if(min_child_label >= PRGC_MAX_DIST_LABEL)
	{
		DEBUG_STATEMENT(PRINTF("My nearest neighbor is MAX_DIST_LABEL away. No path to sink available for me - "));
		DEBUG_STATEMENT(PRGCDump::PrintNode(this); PRINTF(" Relabeling to default(infinity)\n"));
		DeclareOrphan();
		return;
	}

	DEBUG_STATEMENT(PRINTF("Relabeling "); PRGCDump::PrintNode(this); PRINTF(" from %d to %d\n", _distanceLabel, min_child_label+1));
	_distanceLabel = min_child_label + 1;
}

void PRGCNode::SetAsSink(PRGCPNodeList& active_list)
{
	if(!_active){
		_active = true;
		active_list.push_back(this);
	}
	_distanceLabel = TERM_DIST;
	for(int j=0; j<_numCliques; ++j)
	{
		PRGCEdgeInfo* edge_info = _edgeInfo[j];
		PRGCClique* clique = edge_info->_clique;
		clique->_labeledNodesCode |= INDEX_CODE_ONE[edge_info->_indexInClique];
	}
}

void PRGCNode::DeclareOrphan()
{
	_distanceLabel = DIST_DEFAULT;
	for(int j=0; j<_numCliques; ++j)
	{
		PRGCEdgeInfo& edge_info = *(_edgeInfo[j]);
		PRGCClique* clique = edge_info._clique;
		clique->_labeledNodesCode &= INDEX_CODE_ZERO[edge_info._indexInClique];
	}
}

void PRGCClique::RecalculateTightEdgeCodes(int index)
{
	PRGCEdgeInfo& edge_info = _edgeInfo[index];
	edge_info._tightEdgeIntersectionCode = CODE_ALL_ONES_WORD;
	edge_info._tightEdgeUnionCode = 0;

	int *constraints = CONTAINING_CONSTRAINTS[index];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints[i];
		if(_constraintSlacks[constraint_index] > EPSILON)
			continue;
		ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];
		CODE_TYPE constraint_labeling = constraint_info._labeling;
		edge_info._tightEdgeIntersectionCode &= constraint_labeling;
		edge_info._tightEdgeUnionCode |= constraint_labeling;
	}
}

void PRGCClique::RecalculateSaturatedEdgeCode()
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		PRGCEdgeInfo& edge_info = _edgeInfo[i];
		if(edge_info._tightEdgeUnionCode != 0)
			_saturatedEdgeCode |= (1 << i);
	}
}

void PRGCClique::DeclareConstraintNotTight(int constraint_index, CODE_TYPE constraint_labeling,
										   PRGCPNodeList& active_list)
{		
	DEBUG_STATEMENT(PRINTF("Flow constraint %d-%s of ", constraint_index, BINARY(constraint_labeling)));
	DEBUG_STATEMENT(PRGCDump::PrintClique(this); PRINTF(" became relaxed \n"));

	//DEBUG_COND(ITER_NUM == 576);

	bool mark_actives = false;
	if(((_labeledNodesCode ^ CODE_ALL_ONES_WORD) | constraint_labeling) != constraint_labeling)
		mark_actives = true;

	while(constraint_labeling != 0)
	{
		uchar j = FIRST_NON_ZERO(constraint_labeling);
		constraint_labeling &= INDEX_CODE_ZERO[j];
		PRGCEdgeInfo& edge_info = _edgeInfo[j];
		RecalculateTightEdgeCodes(j);
		if(edge_info._tightEdgeUnionCode == 0){
			_saturatedEdgeCode &= INDEX_CODE_ZERO[j];
		}
		if(mark_actives){
			if(!edge_info._node->_active){
				edge_info._node->_active = true;
				active_list.push_back(edge_info._node);
			}
		}
	}
}

void PRGCClique::DeclareConstraintTight(int constraint_index, CODE_TYPE constraint_labeling)
{
	DEBUG_STATEMENT(PRINTF("Flow constraint %d-%s of ", constraint_index, BINARY(constraint_labeling)));
	DEBUG_STATEMENT(PRGCDump::PrintClique(this); PRINTF(" became tight \n"));

	//DEBUG_COND(ITER_NUM == 6);
	CODE_TYPE temp_labeling = constraint_labeling;
	while(temp_labeling != 0)
	{
		int j = FIRST_NON_ZERO(temp_labeling);
		temp_labeling &= INDEX_CODE_ZERO[j];
		PRGCEdgeInfo& edge_info = _edgeInfo[j];
		bool slack_earlier = (edge_info._tightEdgeUnionCode == 0);
		edge_info._tightEdgeIntersectionCode &= constraint_labeling;
		edge_info._tightEdgeUnionCode |= constraint_labeling;
		if(slack_earlier)
			_saturatedEdgeCode |= INDEX_CODE_ONE[j];
	}
}

void PRGCClique::Reparametrize(PRGCPNodeList& excess_list, PRGCPNodeList& active_list, double& flow)
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)		
	{
		if(_constraintSlacks[NUM_CONSTRAINTS-1] <= EPSILON)
			break;

		float min_slack = FindMinContainingSlack(i);
		if(min_slack <= EPSILON)
			continue;

		PRGCNode* pix_node = _edgeInfo[i]._node;

		//DEBUG_STATEMENT(PRINTF("Reparamterizing node "); PRGCDump::PrintNode(pix_node); PRINTF(" by %f\n", min_slack));

		if(pix_node->_excess < 0)
			flow += MIN(min_slack, -pix_node->_excess);

		pix_node->_excess += min_slack;

		int *constraints = CONTAINING_CONSTRAINTS[i];
		for(int j=0; j<NUM_CONTAINING_CONSTRAINTS; ++j)
		{
			int constraint_index = constraints[j];
			float& slack = _constraintSlacks[constraint_index];
			slack -= min_slack;
		}
	}

	for(int i=0; i<GG_CLIQUE_SIZE; ++i)		
	{
		if(_constraintSlacks[0] <= EPSILON)
			break;

		float min_slack = FindMinNonContainingSlack(i);
		if(min_slack <= EPSILON)
			continue;

		PRGCNode* pix_node = _edgeInfo[i]._node;

		//DEBUG_STATEMENT(PRINTF("Reparamterizing node "); PRGCDump::PrintNode(pix_node); PRINTF(" by %f\n", min_slack));

		if(pix_node->_excess > 0)
			flow += MIN(min_slack, pix_node->_excess);

		pix_node->_excess -= min_slack;

		int *constraints = NON_CONTAINING_CONSTRAINTS[i];
		for(int j=0; j<NUM_CONTAINING_CONSTRAINTS; ++j)
		{
			int constraint_index = constraints[j];
			float& slack = _constraintSlacks[constraint_index];
			slack -= min_slack;
		}
	}
}

float PRGCClique::FindMinContainingSlack( int node_index )
{
	int *constraints = CONTAINING_CONSTRAINTS[node_index];
	float min_slack = 1E10;
	for(int j=0; j<NUM_CONTAINING_CONSTRAINTS; ++j)
	{
		int constraint_index = constraints[j];
		float slack = _constraintSlacks[constraint_index];
		if(slack < min_slack)
			min_slack = slack;
	}
	return min_slack;
}

float PRGCClique::FindMinNonContainingSlack( int node_index )
{
	int *constraints = NON_CONTAINING_CONSTRAINTS[node_index];
	float min_slack = 1E10;
	for(int j=0; j<NUM_CONTAINING_CONSTRAINTS; ++j)
	{
		int constraint_index = constraints[j];
		float slack = _constraintSlacks[constraint_index];
		if(slack < min_slack)
			min_slack = slack;
	}
	return min_slack;
}

