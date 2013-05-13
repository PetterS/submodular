#include "stdafx.h"
#include "bkgc/BKGCDataStructures.h"
#include "bkgc/BKGC.h"
#include "bkgc/BKGCDump.h"

void BKGCAuxNode::DeclareActives(std::list<BKGCNode*> *activeNodesList)
{
	int unsaturated = _saturatedEdgeCode ^ CODE_ALL_ONES_WORD;
	int sink_labeled_unsaturated = _sinkLabeledNodesCode & unsaturated;
	int non_children = _childrenCode ^ CODE_ALL_ONES_WORD;
	int	active_code = (sink_labeled_unsaturated | _srcLabeledNodesCode) & non_children;
	while(active_code !=0)
	{
		int i = FIRST_NON_ZERO_WORD[active_code];
		active_code &= (1 << i) ^ CODE_ALL_ONES_WORD;
		BKGCPixNode* pix_node = _nodeInfo[i]._node;
		if(!pix_node->_active){
			pix_node->_active = true;
			activeNodesList->push_back(pix_node);
		}
	}
}

void BKGCPixNode::DeclareActives(std::list<BKGCNode*> *activeNodesList, int aux_index)
{
	BKGCAuxNodeInfo& aux_info = _auxNodeInfo[aux_index];
	BKGCAuxNode* aux_node = aux_info._auxNode;
	int my_index = aux_info._pixNodeIndexInAuxNode;
	BKGCPixNodeInfo& my_info = aux_node->_nodeInfo[my_index];

	if(aux_node->_parentIndex != -1 && 
		((aux_node->_tree == TREE_SRC && my_info._numTightConstraints == 0) || aux_node->_tree==TREE_SINK))
	{
		if(!aux_node->_active){
			aux_node->_active = true;
			activeNodesList->push_back(aux_node);
		}
	}

	if(my_info._numTightConstraints > 0)
	{
		int src_labeled_isec_partners = aux_node->_srcLabeledNodesCode & my_info._tightEdgeIntersectionCode;
		int sink_labeled_uni_partners = aux_node->_sinkLabeledNodesCode & my_info._tightEdgeUnionCode;
		int non_children = aux_info._childrenCode ^ CODE_ALL_ONES_WORD;
		int my_code = (1 << my_index) ;
		int not_me = my_code ^ CODE_ALL_ONES_WORD;
		int active_code = (src_labeled_isec_partners | sink_labeled_uni_partners) 
			& non_children & not_me;
		while(active_code != 0)
		{
			int k = FIRST_NON_ZERO_WORD[active_code];
			active_code &= (1 << k) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo parent_node_info = aux_node->_nodeInfo[k];
			BKGCPixNode* parent_node = parent_node_info._node;
			//MY_ASSERT(parent_node_info._numTightConstraints > 0);
			MY_ASSERT(parent_node->_tree != TREE_NONE);
			if(parent_node->_tree == TREE_SINK && 
				(parent_node_info._tightEdgeIntersectionCode & my_code)==0)
				continue;
			if(!parent_node->_active){
				parent_node->_active = true;
				activeNodesList->push_back(parent_node);
			}
		}
	}

}

bool BKGCPixNode::HasPathToRoot()
{
	if(_pathCheckTimeStamp == TIME)
		return _pathCheckValidPath;

	_pathCheckTimeStamp = TIME;

	if(fabs(_excess) > EPSILON){
		_pathCheckValidPath = true;
		return true;
	}

	_pathCheckValidPath = false;

	if(_parentAuxNodeIndex == -1)
		return false;

	BKGCAuxNode* aux_node = _auxNodeInfo[_parentAuxNodeIndex]._auxNode;

	BKGCPixNode* parent = NULL;
	if(_parentPixNodeIndex == -1)
	{
		if(aux_node->_parentIndex == -1)
			return false;
		if(aux_node->_nodeInfo[aux_node->_parentIndex]._numTightConstraints > 0)
			return false;
		parent = aux_node->_nodeInfo[aux_node->_parentIndex]._node;
	}else
	{
		parent = aux_node->_nodeInfo[_parentPixNodeIndex]._node;
	}
	_pathCheckValidPath = parent->HasPathToRoot();		
	return _pathCheckValidPath;
}

void BKGCAuxNode::RecalculateTightEdgeCodes(int index)
{
	BKGCPixNodeInfo& node_info = _nodeInfo[index];
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

void BKGCAuxNode::AugmentFlow(int sending_node, int receiving_node, float flow)
{
	MY_ASSERT(sending_node != -1 && receiving_node != -1);
	CODE_TYPE sending_node_code = INDEX_CODE_ONE[sending_node];
	CODE_TYPE receiving_node_code = INDEX_CODE_ONE[receiving_node];
	
	std::vector<BKGCPixNode*> pix_node_for_which_to_mark_actives;
	int *constraints;

	//DEBUG_COND(ITER_NUM == 7646);

	constraints = CONTAINING_CONSTRAINTS[sending_node];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints[i];
		if(constraint_index == 0 || constraint_index == CODE_ALL_ONES_WORD)
			continue;
		int constraint_labeling = CONSTRAINT_LABELING(constraint_index);
		if((constraint_labeling & receiving_node_code) != 0)
			continue;
		float& slack = _constraintSlacks[constraint_index];
		bool previously_tight = (slack <= EPSILON);
		slack += flow;
		if(!previously_tight || slack <= EPSILON)
			continue;
		DEBUG_STATEMENT(PRINTF("Flow constraint %d of ",constraint_labeling));
		DEBUG_STATEMENT(BKGCDump::PrintNode(this); PRINTF(" became relaxed \n"));
		//DEBUG_COND(ITER_NUM == 21);

		if(this->_tree == TREE_SRC && !this->_active){
			this->_active = true;
			BKGC_ACTIVE_NODES->push_back(this);
		}
		while(constraint_labeling != 0)
		{
			uchar j = FIRST_NON_ZERO_WORD[constraint_labeling];
			constraint_labeling &= (1 << j) ^ CODE_ALL_ONES_WORD;

			RecalculateTightEdgeCodes(j);
			BKGCPixNodeInfo& node_info = _nodeInfo[j];
			if(node_info._node->_tree == TREE_SINK){
				if(!node_info._node->_active){
					node_info._node->_active = true;
					BKGC_ACTIVE_NODES->push_back(node_info._node);
				}
			}else if(node_info._node->_tree == TREE_NONE){
				pix_node_for_which_to_mark_actives.push_back(node_info._node);
			}

			if(node_info._numTightConstraints == 0)
				_saturatedEdgeCode &= ((1 << j) ^ CODE_ALL_ONES_WORD);
		}
	}

	BKGCPixNodeInfo& receiving_node_info = _nodeInfo[receiving_node];
	constraints = CONTAINING_CONSTRAINTS[receiving_node];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints[i];
		if(constraint_index == 0 || constraint_index == CODE_ALL_ONES_WORD)
			continue;
		int constraint_labeling = CONSTRAINT_LABELING(constraint_index);
		if((constraint_labeling & sending_node_code) != 0)
			continue;
		float& slack = _constraintSlacks[constraint_index];
		bool previously_tight = slack <= EPSILON;
		slack -= flow;
		if(previously_tight || slack >= EPSILON)
			continue;

		DEBUG_STATEMENT(PRINTF("Flow constraint %d of ",constraint_labeling));
		DEBUG_STATEMENT(BKGCDump::PrintNode(this); PRINTF(" became tight \n"));
		//DEBUG_COND(ITER_NUM == 6);
		int temp_labeling = constraint_labeling;
		while(temp_labeling != 0)
		{
			uchar j = FIRST_NON_ZERO_WORD[temp_labeling];
			temp_labeling &= (1 << j) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& node_info = _nodeInfo[j];
			++node_info._numTightConstraints;
			if(!node_info._node->_active){
				node_info._node->_active = true;
				BKGC_ACTIVE_NODES->push_back(node_info._node);
			}
			if(node_info._numTightConstraints == 1 && node_info._node->_tree == TREE_NONE)
				pix_node_for_which_to_mark_actives.push_back(node_info._node);
			node_info._tightEdgeIntersectionCode &= constraint_labeling;
			node_info._tightEdgeUnionCode |= constraint_labeling;
			if(node_info._numTightConstraints == 1)
				 _saturatedEdgeCode |= (1 << j);
		}
		BKGC_TIGHT_CONSTRAINTS->push_back(std::pair<BKGCAuxNode*, int>(this, constraint_index));
	}

	int size = (int)pix_node_for_which_to_mark_actives.size();
	for(int i=0; i<size; ++i){
		BKGCPixNode* pix_node = pix_node_for_which_to_mark_actives[i];
		for(int j=0; j<pix_node->_numAuxNodes; ++j)
			pix_node->DeclareActives(BKGC_ACTIVE_NODES, j);
	}
}

void BKGCAuxNode::Reparametrize(double& flow)
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)		
	{
		if(_constraintSlacks[NUM_CONSTRAINTS-1] <= EPSILON)
			break;

		float min_slack = FindMinContainingSlack(i);
		if(min_slack <= EPSILON)
			continue;

		BKGCPixNode* pix_node = _nodeInfo[i]._node;

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

		BKGCPixNode* pix_node = _nodeInfo[i]._node;

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

float BKGCAuxNode::FindMinContainingSlack( int node_index )
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

float BKGCAuxNode::FindMinNonContainingSlack( int node_index )
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

void BKGCAuxNode::RecalculateSaturatedEdgeCode()
{
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		BKGCPixNodeInfo& node_info = _nodeInfo[i];
		if(node_info._numTightConstraints > 0)
			_saturatedEdgeCode |= (1 << i);
	}
}
