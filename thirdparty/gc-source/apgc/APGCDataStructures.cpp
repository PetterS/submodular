#include "StdAfx.h"
#include "apgc/APGCDataStructures.h"
#include "apgc/APGC.h"
#include "apgc/APGCDump.h"
#include "apgc/APGCDump.h"

bool APGCPixNode::HasPathToRoot()
{
	if(_pathCheckTimeStamp == TIME)
		return _pathCheckValidPath;

	_pathCheckTimeStamp = TIME;

	if(_excess < -EPSILON){
		_pathCheckValidPath = true;
		return true;
	}

	_pathCheckValidPath = false;

	if(_parentAuxNodeIndex == -1)
		return false;

	APGCAuxNode* aux_node = _auxNodeInfo[_parentAuxNodeIndex]._auxNode;

	APGCPixNode* parent = NULL;
	if(_parentPixNodeIndex == -1)
	{
		_pathCheckValidPath = aux_node->HasPathToRoot();		
		return _pathCheckValidPath;
	}else
	{
		parent = aux_node->_nodeInfo[_parentPixNodeIndex]._node;
		_pathCheckValidPath = parent->HasPathToRoot();		
		return _pathCheckValidPath;
	}
}

bool APGCAuxNode::HasPathToRoot()
{
	if(_pathCheckTimeStamp == TIME)
		return _pathCheckValidPath;

	_pathCheckTimeStamp = TIME;
	_pathCheckValidPath = false;

	if(_parentPixNodeIndex == -1)
		return false;

	APGCPixNode* parent = _nodeInfo[_parentPixNodeIndex]._node;
	_pathCheckValidPath = parent->HasPathToRoot();		
	return _pathCheckValidPath;
}

void APGCAuxNode::AugmentFlow(int sending_node, int receiving_node, float flow)
{
	CODE_TYPE sending_node_code = INDEX_CODE_ONE[sending_node];
	CODE_TYPE receiving_node_code = INDEX_CODE_ONE[receiving_node];

	//DEBUG_COND(ITER_NUM == 7646);

	int *constraints;

	MY_ASSERT(sending_node != -1 && receiving_node != -1);

	constraints = CONTAINING_CONSTRAINTS[sending_node];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints[i];
		ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];
		CODE_TYPE constraint_labeling = constraint_info._labeling;
		if((constraint_labeling & receiving_node_code) != 0)
			continue;
		float& slack = _constraintSlacks[constraint_index];
		bool previously_tight = (slack <= EPSILON);
		slack += flow;
		if(!previously_tight || slack <= EPSILON)
			continue;
		DeclareConstraintNotTight(constraint_index, constraint_labeling);
	}

	constraints = CONTAINING_CONSTRAINTS[receiving_node];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints[i];
		ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];
		CODE_TYPE constraint_labeling = constraint_info._labeling;
		if((constraint_labeling & sending_node_code) != 0)
			continue;
		float& slack = _constraintSlacks[constraint_index];
		bool previously_tight = (slack <= EPSILON);
		slack -= flow;
		if(previously_tight || slack > EPSILON)
			continue;
		DeclareConstraintTight(constraint_index, constraint_labeling);
	}
}

void APGCAuxNode::DeclareConstraintTight(int constraint_index, CODE_TYPE constraint_labeling)
{
	DEBUG_STATEMENT(PRINTF("Flow constraint %d-%s of ", constraint_index, BINARY(constraint_labeling)));
	DEBUG_STATEMENT(APGCDump::PrintNode(this); PRINTF(" became tight \n"));

	//DEBUG_COND(ITER_NUM == 6);
	CODE_TYPE temp_labeling = constraint_labeling;
	while(temp_labeling != 0)
	{
		uchar j = FIRST_NON_ZERO(temp_labeling);
		temp_labeling &= (1 << j) ^ CODE_ALL_ONES_WORD;
		APGCPixNodeInfo& node_info = _nodeInfo[j];
		++node_info._numTightConstraints;
		node_info._tightEdgeIntersectionCode &= constraint_labeling;
		node_info._tightEdgeUnionCode |= constraint_labeling;
		if(node_info._numTightConstraints == 1)
			_saturatedEdgeCode |= (1 << j);
	}
	TIGHT_CONSTRAINTS.push_back(std::pair<APGCAuxNode*, int>(this, constraint_index));
}

void APGCAuxNode::DeclareConstraintNotTight(int constraint_index, CODE_TYPE constraint_labeling)
{		
	DEBUG_STATEMENT(PRINTF("Flow constraint %d-%s of ", constraint_index, BINARY(constraint_labeling)));
	DEBUG_STATEMENT(APGCDump::PrintNode(this); PRINTF(" became relaxed \n"));
	//DEBUG_COND(ITER_NUM == 21);

	while(constraint_labeling != 0)
	{
		uchar j = FIRST_NON_ZERO(constraint_labeling);
		constraint_labeling &= (1 << j) ^ CODE_ALL_ONES_WORD;
		APGCPixNodeInfo& node_info = _nodeInfo[j];
		RecalculateTightEdgeCodes(j);
		if(node_info._numTightConstraints == 0)
			_saturatedEdgeCode &= ((1 << j) ^ CODE_ALL_ONES_WORD);
	}
}

void APGCAuxNode::RecalculateTightEdgeCodes(int index)
{
	APGCPixNodeInfo& node_info = _nodeInfo[index];
	node_info._tightEdgeIntersectionCode = CODE_ALL_ONES_WORD;
	node_info._tightEdgeUnionCode = 0;
	node_info._numTightConstraints = 0;

	int *constraints = CONTAINING_CONSTRAINTS[index];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints[i];
		if(constraint_index == 0 || constraint_index == CODE_ALL_ONES_WORD)
			continue;
		if(_constraintSlacks[constraint_index] > EPSILON)
			continue;
		node_info._numTightConstraints++;
		ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];
		CODE_TYPE constraint_labeling = constraint_info._labeling;
		node_info._tightEdgeIntersectionCode &= constraint_labeling;
		node_info._tightEdgeUnionCode |= constraint_labeling;
	}
}

void APGCAuxNode::Reparametrize(double& flow)
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)		
	{
		if(_constraintSlacks[NUM_CONSTRAINTS-1] <= EPSILON)
			break;

		float min_slack = FindMinContainingSlack(i);
		if(min_slack <= EPSILON)
			continue;

		APGCPixNode* pix_node = _nodeInfo[i]._node;

		//DEBUG_STATEMENT(PRINTF("Reparamterizing node "); APGCDump::PrintNode(pix_node); PRINTF(" by %f\n", min_slack));

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

		APGCPixNode* pix_node = _nodeInfo[i]._node;

		//DEBUG_STATEMENT(PRINTF("Reparamterizing node "); APGCDump::PrintNode(pix_node); PRINTF(" by %f\n", min_slack));

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

	MY_ASSERT(_constraintSlacks[NUM_CONSTRAINTS-1] <= EPSILON);
	MY_ASSERT(_constraintSlacks[0] <= EPSILON);
}

float APGCAuxNode::FindMinContainingSlack( int node_index )
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

float APGCAuxNode::FindMinNonContainingSlack( int node_index )
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

void APGCAuxNode::RecalculateSaturatedEdgeCode()
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		APGCPixNodeInfo& node_info = _nodeInfo[i];
		if(node_info._numTightConstraints > 0)
			_saturatedEdgeCode |= (1 << i);
	}
}
