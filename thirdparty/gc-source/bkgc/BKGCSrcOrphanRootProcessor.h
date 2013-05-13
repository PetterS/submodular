#pragma once

#include "bkgc/BKGCDataStructures.h"
#include "bkgc/BKGCDebug.h"
#include "bkgc/BKGCDump.h"
#include "bkgc/BKGCOrphanUtil.h"

class BKGCSrcOrphanRootProcessor : public BKGCOrphanUtil
{
public:
	inline void FindPath(BKGCPixNode* root, std::list<BKGCNode*> *activeNodesList)
	{
		++TIME;
		root->_pathCheckTimeStamp = TIME;
		root->_pathCheckValidPath = false;
		int old_parent_aux_index, old_parent_pix_index;
		BKGCNode* path_found_node = FindPathInTheTree(root, old_parent_aux_index, old_parent_pix_index);
		if(path_found_node)
			InvertTree(activeNodesList, root, path_found_node, old_parent_aux_index, old_parent_pix_index);
		else
			DeclareTreeOrphan(TREE_SRC, root, activeNodesList);
	}

	inline void FindPath(BKGCAuxNode* aux_node, std::list<BKGCNode*> *activeNodesList)
	{
		++TIME;
		aux_node->_parentIndex = -1;
		int old_parent_aux_index, old_parent_pix_index;
		BKGCNode* path_found_node = FindPathInTheTree(aux_node, old_parent_aux_index, old_parent_pix_index);
		if(path_found_node)
			InvertTree(activeNodesList, aux_node, path_found_node, old_parent_aux_index, old_parent_pix_index);
		else
			DeclareTreeOrphan(TREE_SRC, aux_node, activeNodesList);
	}

private:
	BKGCNode* FindPathInTheTree(BKGCNode* root, int& old_parent_aux_index, int& old_parent_pix_index);
	void InvertTree(std::list<BKGCNode*> *activeNodesList, BKGCNode* root, BKGCNode* path_found_node, int old_parent_aux_index, int old_parent_pix_index);
	bool FindPathThroughPixNode(BKGCPixNode* node, int& old_parent_aux_index, int& old_parent_pix_index);
	bool FindPathThroughAuxNode(BKGCAuxNode* node, int& old_parent_pix_index);

public:
	static inline void ResetOutgoingEdges(BKGCAuxNodeInfo& send_node_info, BKGCAuxNode* aux_node)
	{
		//make my children, direct children of aux node
		DEBUG_STATEMENT(BKGCPixNode* send_node = aux_node->_nodeInfo[send_node_info._pixNodeIndexInAuxNode]._node);
		while(send_node_info._childrenCode != 0)
		{
			int child_index = FIRST_NON_ZERO_WORD[send_node_info._childrenCode];
			send_node_info._childrenCode &= (1 << child_index) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& child_node_info = aux_node->_nodeInfo[child_index];
			MY_ASSERT(child_node_info._node->_parentAuxNodeIndex == child_node_info._cliqueIndexInNode);
			MY_ASSERT(child_node_info._node->_parentPixNodeIndex == send_node_info._pixNodeIndexInAuxNode);
			DEBUG_STATEMENT(BKGCDump::UpdateEdge(send_node, child_node_info._node, aux_node, child_node_info._node));
			child_node_info._node->_parentPixNodeIndex = -1;
			aux_node->_childrenCode |= (1 << child_index);
		}
	}

	static inline void ResetOutgoingEdges(BKGCAuxNodeInfo& current_send_aux_info, BKGCAuxNodeInfo& new_send_aux_info, 
		int new_send_index, BKGCAuxNode* aux_node)
	{
		//make my children, direct children of my parent
		DEBUG_STATEMENT(BKGCPixNodeInfo& new_send_pix_info = aux_node->_nodeInfo[new_send_aux_info._pixNodeIndexInAuxNode]);
		DEBUG_STATEMENT(BKGCPixNode* current_send_pix_node = aux_node->_nodeInfo[current_send_aux_info._pixNodeIndexInAuxNode]._node);
		int new_send_index_code = 1 << new_send_index;
		int children_code = current_send_aux_info._childrenCode;
		while(children_code != 0)
		{
			int child_index = FIRST_NON_ZERO_WORD[children_code];
			int invert_code = (1 << child_index) ^ CODE_ALL_ONES_WORD;
			children_code &= invert_code;
			BKGCPixNodeInfo& child_node_info = aux_node->_nodeInfo[child_index];
			if((child_node_info._tightEdgeIntersectionCode & new_send_index_code) == 0)
				continue;
			current_send_aux_info._childrenCode &= invert_code;
			MY_ASSERT(child_node_info._node->_parentAuxNodeIndex == child_node_info._cliqueIndexInNode);
			MY_ASSERT(child_node_info._node->_parentPixNodeIndex == current_send_aux_info._pixNodeIndexInAuxNode);
			//MY_ASSERT((child_node_info._tightEdgeIntersectionCode & (1 << new_send_index)) != 0);
			DEBUG_STATEMENT(BKGCDump::UpdateEdge(current_send_pix_node, child_node_info._node,new_send_pix_info._node, child_node_info._node));
			child_node_info._node->_parentPixNodeIndex = new_send_index;
			new_send_aux_info._childrenCode |= (1 << child_index);
		}
	}

};
