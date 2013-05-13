#include "StdAfx.h"
#include "bkgc/BKGCDebug.h"
#include "bkgc/BKGC.h"
#include "test/CliquePotential.h"

#ifdef CHECK_FOR_ERRORS

void BKGCDebug::ValidateAll(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes)
{
	ValidateConstraints(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);
	ValidateSaturatedEdges(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);
	ValidateChildren(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);
	ValidatePaths(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);	
}

void BKGCDebug::ValidateConstraints(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		BKGCPixNode& node = pix_nodes[i];

		for(int j=0; j<node._numAuxNodes; ++j)
		{
			BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[j];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			BKGCPixNodeInfo& node_info = aux_info._auxNode->_nodeInfo[aux_info._pixNodeIndexInAuxNode];
			uchar num_tight_constraints = 0;
			int union_code = 0, intersection_code = CODE_ALL_ONES_WORD;

			int* constraints = CONTAINING_CONSTRAINTS[aux_info._pixNodeIndexInAuxNode];
			//DEBUG_COND(ITER_NUM >= 7645 && node._x == )
			for(int l=0; l<NUM_CONTAINING_CONSTRAINTS; ++l)
			{
				int constraint_index = constraints[l];
				if(constraint_index == 0 || constraint_index == CODE_ALL_ONES_WORD)
					continue;
				int constraint_labeling = CONSTRAINT_LABELING(constraint_index);
				float slack = aux_node->_constraintSlacks[constraint_index];
				MY_ASSERT(slack >= 0);
				if(slack <= EPSILON)
				{
					++num_tight_constraints;
					union_code |= constraint_labeling;
					intersection_code &= constraint_labeling;
				}
			}
			MY_ASSERT(node_info._tightEdgeIntersectionCode == intersection_code);
			MY_ASSERT(node_info._tightEdgeUnionCode == union_code);
			MY_ASSERT(node_info._numTightConstraints == num_tight_constraints);
		}
	}
}

void BKGCDebug::ValidateSaturatedEdges(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_aux_nodes; ++i)
	{
		BKGCAuxNode& aux_node = aux_nodes[i];
		for(int i=0; i<GG_CLIQUE_SIZE; ++i)
		{
			BKGCPixNodeInfo& node_info = aux_node._nodeInfo[i];
			BKGCPixNode* pix_node = node_info._node;
			if(node_info._numTightConstraints > 0)
				MY_ASSERT((aux_node._saturatedEdgeCode & (1 << i)) != 0);
		}
		int saturated_code = aux_node._saturatedEdgeCode;
		while(saturated_code != 0)
		{
			int j = FIRST_NON_ZERO_WORD[saturated_code];
			saturated_code &= (1 << j) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& node_info = aux_node._nodeInfo[j];
			MY_ASSERT(node_info._numTightConstraints > 0);
		}
	}
}

void BKGCDebug::ValidateChildren(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_aux_nodes; ++i)
	{
		BKGCAuxNode& aux_node = aux_nodes[i];
		if(aux_node._childrenCode != 0)
			MY_ASSERT(aux_node._parentIndex != -1);
		int children_code = aux_node._childrenCode;
		while(children_code != 0)
		{
			int j = FIRST_NON_ZERO_WORD[children_code];
			children_code &= (1 << j) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& pix_node_info = aux_node._nodeInfo[j];
			BKGCPixNode* pix_node = pix_node_info._node;
			MY_ASSERT(pix_node->_parentAuxNodeIndex == pix_node_info._cliqueIndexInNode);
			MY_ASSERT(pix_node->_parentPixNodeIndex == -1);
			MY_ASSERT(pix_node->_tree == aux_node._tree);
		}
	}

	for(int i=0; i<num_pix_nodes; ++i)
	{
		BKGCPixNode& node = pix_nodes[i];
		if(node._excess < -EPSILON)
			MY_ASSERT(node._tree == TREE_SINK);
		else if(node._excess > EPSILON)
			MY_ASSERT(node._tree == TREE_SRC);
		else{ 
			if(node._parentAuxNodeIndex == -1)
				MY_ASSERT(node._tree == TREE_NONE);
			else
				MY_ASSERT(node._tree == TREE_SRC || node._tree == TREE_SINK);
		}

		//DEBUG_COND(ITER_NUM == 13 && node._x == 0 && node._y == 2);

		if(node._parentAuxNodeIndex != -1)
		{
			BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[node._parentAuxNodeIndex];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			if(node._parentPixNodeIndex != -1){
				BKGCPixNodeInfo& parent_pix_info = aux_node->_nodeInfo[node._parentPixNodeIndex];
				BKGCAuxNodeInfo& parent_aux_info = parent_pix_info._node->_auxNodeInfo[parent_pix_info._cliqueIndexInNode];
				MY_ASSERT((parent_aux_info._childrenCode & (1 << aux_info._pixNodeIndexInAuxNode)) != 0);
			}else
			{
				MY_ASSERT((aux_node->_childrenCode & (1 << aux_info._pixNodeIndexInAuxNode)) != 0);
			}
		}

		for(int j=0; j<node._numAuxNodes; ++j)
		{
			BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[j];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			int my_index = aux_info._pixNodeIndexInAuxNode;
			if(node._tree == TREE_SRC)
				MY_ASSERT((aux_node->_srcLabeledNodesCode & (1<<my_index))!= 0);
			else if(node._tree == TREE_SINK)
				MY_ASSERT((aux_node->_sinkLabeledNodesCode & (1<<my_index))!= 0);
			else{
				MY_ASSERT((aux_node->_sinkLabeledNodesCode & (1<<my_index))== 0);
				MY_ASSERT((aux_node->_srcLabeledNodesCode & (1<<my_index))== 0);
				MY_ASSERT(aux_info._childrenCode == 0);
				MY_ASSERT(aux_node->_parentIndex != my_index);
			}
			//if(aux_info._childrenCode != 0)
			//	MY_ASSERT(node._parentAuxNodeIndex != j);
			if(aux_node->_parentIndex == my_index){
				if(node._parentPixNodeIndex==-1)
					MY_ASSERT(node._parentAuxNodeIndex != j);
			}
			int temp_children_code = aux_info._childrenCode;
			while(temp_children_code != 0)
			{
				int k = FIRST_NON_ZERO_WORD[temp_children_code];
				temp_children_code &= (1 << k) ^ CODE_ALL_ONES_WORD;
				BKGCPixNode* child_node = aux_info._auxNode->_nodeInfo[k]._node;
				MY_ASSERT(child_node->_parentPixNodeIndex == my_index);
				MY_ASSERT(child_node->_tree == node._tree);
			}
		}
	}
}

void BKGCDebug::ValidatePaths(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		BKGCPixNode& node = pix_nodes[i];
		if(fabs(node._excess) > EPSILON)
			continue;
		//DEBUG_COND(ITER_NUM == 13 && node._x == 0 && node._y == 0);
		if(node._parentAuxNodeIndex == -1)
		{
			for(int j=0; j<node._numAuxNodes; ++j)
			{
				BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[j];
				BKGCAuxNode* aux_node = aux_info._auxNode;
				int my_index = aux_info._pixNodeIndexInAuxNode;
				BKGCPixNodeInfo& node_info = aux_node->_nodeInfo[my_index];
				if(node_info._numTightConstraints == 0){
					MY_ASSERT(aux_node->_parentIndex == -1 || aux_node->_active);
					continue;
				}
				if(aux_node->_tree == TREE_SINK){
					MY_ASSERT(aux_node->_parentIndex == -1 || aux_node->_active);
					continue;
				}
				ValidateNoP2PPath(node_info, my_index, aux_node);
			}
		}else
		{
			BKGCAuxNodeInfo& aux_info = node._auxNodeInfo[node._parentAuxNodeIndex];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			int my_index = aux_info._pixNodeIndexInAuxNode;
			BKGCPixNodeInfo& node_info = aux_node->_nodeInfo[my_index];
			MY_ASSERT(node._tree != TREE_NONE);
			if(node._parentPixNodeIndex == -1)
			{
				MY_ASSERT(aux_node->_parentIndex != -1);
				MY_ASSERT((aux_node->_childrenCode & (1 << aux_info._pixNodeIndexInAuxNode)) != 0);
				MY_ASSERT(aux_node->_tree == node._tree);
				//MY_ASSERT(aux_info._childrenCode == 0);
				if(aux_node->_tree == TREE_SRC)
					MY_ASSERT(node_info._numTightConstraints == 0);
			}
			else
			{
				int parent_index = node._parentPixNodeIndex;
				BKGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[parent_index];
				BKGCPixNode* parent_node = parent_node_info._node;
				if(node._tree == TREE_SINK)
				{
					MY_ASSERT(aux_node->_parentIndex != my_index);
					//MY_ASSERT(aux_node->_parentIndex != parent_index);
					MY_ASSERT((parent_node_info._tightEdgeIntersectionCode & (1 << my_index)) != 0);
					MY_ASSERT(parent_node->_excess < 0 || parent_node->_parentAuxNodeIndex != -1);
				}else
				{
					MY_ASSERT((aux_node->_childrenCode & (1 << my_index)) == 0);
					//MY_ASSERT((aux_node->_childrenCode & (1 << parent_index)) == 0);
					MY_ASSERT((node_info._tightEdgeIntersectionCode & (1 << parent_index)) != 0);
					MY_ASSERT(parent_node->_excess > 0 || parent_node->_parentAuxNodeIndex != -1);
				}
			}
		}
	}
}

void BKGCDebug::ValidateNoP2PPath(BKGCPixNodeInfo& node_info, int my_index, BKGCAuxNode* aux_node)
{
	for(int k=0; k<GG_CLIQUE_SIZE; ++k)
	{
		BKGCPixNodeInfo& possible_parent_info = aux_node->_nodeInfo[k];
		BKGCPixNode* possible_parent = possible_parent_info._node;
		if(possible_parent->_active)
			continue;
		if(possible_parent->_tree == TREE_SRC)
			MY_ASSERT((node_info._tightEdgeIntersectionCode & (1 << k)) == 0);
		else if(possible_parent->_tree == TREE_SINK)
		{
			if(possible_parent->_parentAuxNodeIndex != -1)
			{
				BKGCAuxNode* parent_aux_node = possible_parent->_auxNodeInfo[possible_parent->_parentAuxNodeIndex]._auxNode;
				if(parent_aux_node == aux_node)
					continue;
			}
			MY_ASSERT((possible_parent_info._tightEdgeIntersectionCode & (1 << my_index)) == 0);
		}
	}		
}
#endif
