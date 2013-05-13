#include "StdAfx.h"
#include "apgc/APGCDebug.h"
#include "apgc/APGC.h"

#ifdef CHECK_FOR_ERRORS

void APGCDebug::ValidateAll(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	++TIME;

	ValidateConstraints(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);
	ValidateSaturatedEdges(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);
	ValidateChildren(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);
	ValidatePaths(pix_nodes, num_pix_nodes, aux_nodes, num_aux_nodes);	
	ValidateDistanceLabel(pix_nodes, num_pix_nodes);
	ValidateDistanceLabel(aux_nodes, num_aux_nodes);
}

void APGCDebug::ValidateConstraints(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		APGCPixNode& node = pix_nodes[i];

		for(int j=0; j<node._numAuxNodes; ++j)
		{
			APGCAuxNodeInfo& aux_info = node._auxNodeInfo[j];
			APGCAuxNode* aux_node = aux_info._auxNode;
			APGCPixNodeInfo& node_info = aux_info._auxNode->_nodeInfo[aux_info._pixNodeIndexInAuxNode];
			uchar num_tight_constraints = 0;
			CODE_TYPE union_code = 0, intersection_code = CODE_ALL_ONES_WORD;
			//DEBUG_COND(ITER_NUM >= 7645 && node._x == )

			int* constraints_list = CONTAINING_CONSTRAINTS[aux_info._pixNodeIndexInAuxNode];
			for(int l=0; l<NUM_CONTAINING_CONSTRAINTS; ++l)
			{
				int constraint_index = constraints_list[l];
				if(constraint_index == 0 || constraint_index == CODE_ALL_ONES_WORD)
					continue;
				ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];
				CODE_TYPE constraint_labeling = constraint_info._labeling;
				float constraint_slack = aux_node->_constraintSlacks[constraint_index];
				MY_ASSERT(constraint_slack >= 0);
				if(constraint_slack <= EPSILON)
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

void APGCDebug::ValidateSaturatedEdges(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_aux_nodes; ++i)
	{
		APGCAuxNode& aux_node = aux_nodes[i];
		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
		{
			APGCPixNodeInfo& node_info = aux_node._nodeInfo[j];
			APGCPixNode* pix_node = node_info._node;
			if(node_info._numTightConstraints > 0)
				MY_ASSERT((aux_node._saturatedEdgeCode & (1 << j)) != 0);
		}
		CODE_TYPE saturated_code = aux_node._saturatedEdgeCode;
		while(saturated_code != 0)
		{
			CODE_TYPE j = FIRST_NON_ZERO(saturated_code);
			saturated_code &= (1 << j) ^ CODE_ALL_ONES_WORD;
			APGCPixNodeInfo& node_info = aux_node._nodeInfo[j];
			MY_ASSERT(node_info._numTightConstraints > 0);
		}
	}
}

void APGCDebug::ValidateChildren(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_aux_nodes; ++i)
	{
		APGCAuxNode& aux_node = aux_nodes[i];
		CODE_TYPE children_code = aux_node._childrenCode;
		while(children_code != 0)
		{
			int j = FIRST_NON_ZERO(children_code);
			children_code &= (1 << j) ^ CODE_ALL_ONES_WORD;
			APGCPixNodeInfo& pix_node_info = aux_node._nodeInfo[j];
			APGCPixNode* pix_node = pix_node_info._node;
			MY_ASSERT(pix_node->_parentAuxNodeIndex == pix_node_info._cliqueIndexInNode);
			MY_ASSERT(pix_node->_parentPixNodeIndex == -1);
			MY_ASSERT(pix_node->_tree == aux_node._tree);
			if(aux_node._tree == TREE_SRC)
				MY_ASSERT((aux_node._saturatedEdgeCode & INDEX_CODE_ONE[j]) == 0);
		}
	}

	for(int i=0; i<num_pix_nodes; ++i)
	{
		APGCPixNode& node = pix_nodes[i];
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
			APGCAuxNodeInfo& aux_info = node._auxNodeInfo[node._parentAuxNodeIndex];
			APGCAuxNode* aux_node = aux_info._auxNode;
			if(node._parentPixNodeIndex != -1){
				APGCPixNodeInfo& parent_pix_info = aux_node->_nodeInfo[node._parentPixNodeIndex];
				APGCAuxNodeInfo& parent_aux_info = parent_pix_info._node->_auxNodeInfo[parent_pix_info._cliqueIndexInNode];
				MY_ASSERT((parent_aux_info._childrenCode & (1 << aux_info._pixNodeIndexInAuxNode)) != 0);
			}else
			{
				MY_ASSERT((aux_node->_childrenCode & (1 << aux_info._pixNodeIndexInAuxNode)) != 0);
			}
		}

		for(int j=0; j<node._numAuxNodes; ++j)
		{
			APGCAuxNodeInfo& aux_info = node._auxNodeInfo[j];
			APGCAuxNode* aux_node = aux_info._auxNode;
			int my_index = aux_info._pixNodeIndexInAuxNode;
			CODE_TYPE my_code = INDEX_CODE_ONE[my_index];
			if(node._tree == TREE_SRC)
				MY_ASSERT((aux_node->_srcLabeledNodesCode & my_code) != 0);
			else if(node._tree == TREE_SINK)
				MY_ASSERT((aux_node->_sinkLabeledNodesCode & my_code) != 0);
			else{
				MY_ASSERT((aux_node->_sinkLabeledNodesCode & my_code) == 0);
				MY_ASSERT((aux_node->_srcLabeledNodesCode & my_code) == 0);
				MY_ASSERT(aux_info._childrenCode == 0);
			}

			if(aux_info._childrenCode != 0)
				MY_ASSERT(node._parentAuxNodeIndex != j);

			if(aux_node->_parentPixNodeIndex == my_index)
				MY_ASSERT(node._parentAuxNodeIndex != j);

			CODE_TYPE temp_children_code = aux_info._childrenCode;
			while(temp_children_code != 0)
			{
				int k = FIRST_NON_ZERO(temp_children_code);
				temp_children_code &= (1 << k) ^ CODE_ALL_ONES_WORD;
				APGCPixNode* child_node = aux_info._auxNode->_nodeInfo[k]._node;
				MY_ASSERT(child_node->_parentPixNodeIndex == my_index);
				MY_ASSERT(child_node->_tree == node._tree);
			}
		}
	}
}

void APGCDebug::ValidatePaths(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		APGCPixNode& node = pix_nodes[i];
		if(fabs(node._excess) > EPSILON)
			continue;
		//DEBUG_COND(ITER_NUM == 13 && node._x == 0 && node._y == 0);
		if(node._parentAuxNodeIndex == -1)
		{
			for(int j=0; j<node._numAuxNodes; ++j)
			{
				APGCAuxNodeInfo& aux_info = node._auxNodeInfo[j];
				APGCAuxNode* aux_node = aux_info._auxNode;
				int my_index = aux_info._pixNodeIndexInAuxNode;
				APGCPixNodeInfo& node_info = aux_node->_nodeInfo[my_index];
				if(node_info._numTightConstraints == 0){
					MY_ASSERT(aux_node->_parentPixNodeIndex == -1 || aux_node->_active);
					continue;
				}
				if(aux_node->_tree == TREE_SINK){
					MY_ASSERT(aux_node->_active);
					continue;
				}
				ValidateNoP2PPath(node_info, my_index, aux_node);
			}
		}else
		{
			APGCAuxNodeInfo& aux_info = node._auxNodeInfo[node._parentAuxNodeIndex];
			APGCAuxNode* aux_node = aux_info._auxNode;
			int my_index = aux_info._pixNodeIndexInAuxNode;
			APGCPixNodeInfo& node_info = aux_node->_nodeInfo[my_index];
			MY_ASSERT(node._tree != TREE_NONE);
			if(node._parentPixNodeIndex == -1)
			{
				MY_ASSERT(node._distanceLabel == (aux_node->_distanceLabel+1));
				MY_ASSERT(aux_node->_parentPixNodeIndex != -1);
				MY_ASSERT(aux_node->_parentPixNodeIndex != my_index);
				MY_ASSERT((aux_node->_childrenCode & (1 << aux_info._pixNodeIndexInAuxNode)) != 0);
				MY_ASSERT(aux_node->_tree == node._tree);
				MY_ASSERT(aux_info._childrenCode == 0);
				if(aux_node->_tree == TREE_SRC)
					MY_ASSERT(node_info._numTightConstraints == 0);
			}
			else
			{
				int parent_index = node._parentPixNodeIndex;
				APGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[parent_index];
				APGCPixNode* parent_node = parent_node_info._node;
				MY_ASSERT(node._distanceLabel == (parent_node->_distanceLabel+2));
				if(node._tree == TREE_SINK)
				{
					MY_ASSERT(aux_node->_parentPixNodeIndex != my_index);
					MY_ASSERT(aux_node->_parentPixNodeIndex != parent_index);
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

void APGCDebug::ValidateNoP2PPath(APGCPixNodeInfo& node_info, int my_index, APGCAuxNode* aux_node)
{
	for(int k=0; k<GG_CLIQUE_SIZE; ++k)
	{
		APGCPixNodeInfo& possible_parent_info = aux_node->_nodeInfo[k];
		APGCPixNode* possible_parent = possible_parent_info._node;
		if(possible_parent->_active)
			continue;
		if(possible_parent->_parentAuxNodeIndex == -1)
			continue;
		if(!possible_parent->HasPathToRoot())
			continue;
		if(possible_parent->_parentAuxNodeIndex == aux_node->_nodeInfo[k]._cliqueIndexInNode)
			continue;
		if(possible_parent->_tree == TREE_SRC)
			MY_ASSERT((node_info._tightEdgeIntersectionCode & INDEX_CODE_ONE[k]) == 0);
		else 
			MY_ASSERT((possible_parent_info._tightEdgeIntersectionCode & INDEX_CODE_ONE[my_index]) == 0);
	}		
}

void APGCDebug::ValidateDistanceLabel( APGCPixNode* pix_nodes, int num_pix_nodes )
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		APGCPixNode* pix_node = pix_nodes + i;

		//DEBUG_COND(ITER_NUM >= 99 && pix_node->_x == 5 && pix_node->_y == 7);

		if(fabs(pix_node->_excess) > EPSILON){
			MY_ASSERT(pix_node->_distanceLabel == TERM_DIST);
			continue;
		}
		if(pix_node->_parentAuxNodeIndex == -1)
		{
			MY_ASSERT(pix_node->_distanceLabel == DIST_DEFAULT);
			MY_ASSERT(pix_node->_tree == TREE_NONE);
			continue;
		}
		if(pix_node->_parentPixNodeIndex == -1)
		{
			MY_ASSERT(pix_node->_tree != TREE_NONE);
			APGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[pix_node->_parentAuxNodeIndex];
			APGCAuxNode* aux_node = aux_info._auxNode;
			MY_ASSERT(aux_node->_distanceLabel != DIST_DEFAULT);
			MY_ASSERT(pix_node->_distanceLabel == (aux_node->_distanceLabel + 1));
			continue;
		}
		for(int j=0; j<pix_node->_numAuxNodes; ++j)
		{
			APGCAuxNodeInfo& my_aux_info = pix_node->_auxNodeInfo[j];
			int my_index = my_aux_info._pixNodeIndexInAuxNode;
			APGCAuxNode* aux_node = my_aux_info._auxNode;
			APGCPixNodeInfo& my_pix_info = aux_node->_nodeInfo[my_index];
			if(pix_node->_tree == TREE_SRC)
			{
				//MY_ASSERT((my_pix_info._numTightConstraints != 0) || (pix_node->_parentAuxNodeIndex != j));
				CODE_TYPE me_or_my_children = my_aux_info._childrenCode | (1 << my_index);
				CODE_TYPE nodes_i_can_seek_label = my_pix_info._tightEdgeIntersectionCode;
				nodes_i_can_seek_label &= (me_or_my_children ^ CODE_ALL_ONES_WORD);
				nodes_i_can_seek_label &= aux_node->_srcLabeledNodesCode;
				while(nodes_i_can_seek_label != 0)
				{
					int send_node_index = FIRST_NON_ZERO(nodes_i_can_seek_label);
					nodes_i_can_seek_label &= INDEX_CODE_ZERO[send_node_index];
					APGCPixNodeInfo& send_node_pix_info = aux_node->_nodeInfo[send_node_index]; 
					APGCPixNode* send_node = send_node_pix_info._node; 
					if(send_node->_parentAuxNodeIndex == send_node_pix_info._cliqueIndexInNode)
						continue;
					if(!send_node->HasPathToRoot())
						continue;
					MY_ASSERT(send_node->_distanceLabel >= (pix_node->_distanceLabel-2));
				}
			}else
			{
				CODE_TYPE nodes_i_can_seek_label = my_pix_info._tightEdgeUnionCode & INDEX_CODE_ZERO[my_index];
				while(nodes_i_can_seek_label != 0)
				{
					int rec_node_index = FIRST_NON_ZERO(nodes_i_can_seek_label);
					nodes_i_can_seek_label &= INDEX_CODE_ZERO[rec_node_index];
					APGCPixNodeInfo& rec_node_pix_info = aux_node->_nodeInfo[rec_node_index]; 
					MY_ASSERT(rec_node_pix_info._numTightConstraints > 0);
					if((rec_node_pix_info._tightEdgeIntersectionCode & INDEX_CODE_ONE[my_index]) == 0)
						continue;
					APGCPixNode* rec_node = rec_node_pix_info._node; 
					if(rec_node->_tree != TREE_SINK)
						continue;
					if(rec_node->_parentAuxNodeIndex == rec_node_pix_info._cliqueIndexInNode)
						continue;
					if(!rec_node->HasPathToRoot())
						continue;
					MY_ASSERT(rec_node->_distanceLabel >= (pix_node->_distanceLabel-2));
				}
			}
			
		}
	}
}

void APGCDebug::ValidateDistanceLabel(APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_aux_nodes; ++i)
	{
		APGCAuxNode* aux_node = aux_nodes + i;
		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
		{
			APGCPixNodeInfo& pix_info = aux_node->_nodeInfo[j];
			APGCPixNode* pix_node = pix_info._node;
			if(aux_node->_parentPixNodeIndex == -1)
			{
				MY_ASSERT(aux_node->_tree == TREE_NONE);
				MY_ASSERT(aux_node->_distanceLabel == DIST_DEFAULT);
				if(pix_node->_active || pix_node->_tree == TREE_NONE)
					continue;
				if(pix_node->_tree == TREE_SRC)
					MY_ASSERT(false);
				if(pix_node->_tree == TREE_SINK)
					MY_ASSERT((aux_node->_saturatedEdgeCode & INDEX_CODE_ONE[j]) != 0);
			}else
			{
				MY_ASSERT(aux_node->_tree != TREE_NONE);
				MY_ASSERT(aux_node->_distanceLabel != DIST_DEFAULT);
				if(pix_node->_active || pix_node->_tree == TREE_NONE)
					continue;
				if(pix_node->_distanceLabel >= (aux_node->_distanceLabel-1))
					continue;
				if(pix_node->_tree == TREE_SRC)
					MY_ASSERT(false);
				if(pix_node->_tree == TREE_SINK)
					MY_ASSERT((aux_node->_saturatedEdgeCode & INDEX_CODE_ONE[j]) != 0);
			}
		}
	}
}
#endif
