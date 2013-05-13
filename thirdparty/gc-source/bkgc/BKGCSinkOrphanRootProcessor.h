#pragma once

#include "bkgc/BKGCDataStructures.h"
#include "bkgc/BKGCDebug.h"
#include "bkgc/BKGCDump.h"
#include "bkgc/BKGCOrphanUtil.h"

class BKGCSinkOrphanRootProcessor : public BKGCOrphanUtil
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
			DeclareTreeOrphan(TREE_SINK, root, activeNodesList);
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
			DeclareTreeOrphan(TREE_SINK, aux_node, activeNodesList);
	}

private:
	BKGCNode* FindPathInTheTree(BKGCNode* root, int& old_parent_aux_index, int& old_parent_pix_index);
	void InvertTree(std::list<BKGCNode*> *activeNodesList, BKGCNode* root, BKGCNode* path_found_node, int old_parent_aux_index, int old_parent_pix_index);
	bool FindPathThroughPixNode(BKGCPixNode* node, int& old_parent_aux_index, int& old_parent_pix_index);
	bool FindPathThroughAuxNode(BKGCAuxNode* node, int& old_parent_pix_index);

public:
	static inline void ResetIncomingEdges(BKGCAuxNodeInfo& rec_node_info, BKGCAuxNode* aux_node)
	{
		DEBUG_STATEMENT(BKGCPixNode* pix_node = aux_node->_nodeInfo[rec_node_info._pixNodeIndexInAuxNode]._node);
		while(rec_node_info._childrenCode != 0)
		{
			int child_index = FIRST_NON_ZERO_WORD[rec_node_info._childrenCode];
			rec_node_info._childrenCode &= (1 << child_index) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& child_node_info = aux_node->_nodeInfo[child_index];
			MY_ASSERT(child_node_info._node->_parentAuxNodeIndex == child_node_info._cliqueIndexInNode);
			MY_ASSERT(child_node_info._node->_parentPixNodeIndex == rec_node_info._pixNodeIndexInAuxNode);
			DEBUG_STATEMENT(PRINTF("Resetting edge "); BKGCDump::PrintNode(child_node_info._node)); 
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(pix_node));
			DEBUG_STATEMENT(PRINTF(" to "); BKGCDump::PrintNode(child_node_info._node)); 
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
			child_node_info._node->_parentPixNodeIndex = -1;
			aux_node->_childrenCode |= (1 << child_index);
		}
	}

	static inline void ResetIncomingEdges(BKGCAuxNodeInfo& current_rec_aux_info, BKGCAuxNodeInfo& new_rec_aux_info, 
		int new_rec_index, BKGCAuxNode* aux_node)
	{
		DEBUG_STATEMENT(BKGCPixNode* current_rec_pix_node = aux_node->_nodeInfo[current_rec_aux_info._pixNodeIndexInAuxNode]._node);
		BKGCPixNodeInfo& new_rec_pix_info = aux_node->_nodeInfo[new_rec_aux_info._pixNodeIndexInAuxNode];
		int children_code = current_rec_aux_info._childrenCode;
		while(children_code != 0)
		{
			int child_index = FIRST_NON_ZERO_WORD[children_code];
			int index_code = (1 << child_index);
			int index_code_inverted = index_code ^ CODE_ALL_ONES_WORD;
			children_code &= index_code_inverted;
			if((new_rec_pix_info._tightEdgeIntersectionCode & index_code) == 0)
				continue;
			current_rec_aux_info._childrenCode &= index_code_inverted;
			BKGCPixNodeInfo& child_node_info = aux_node->_nodeInfo[child_index];
			MY_ASSERT(child_node_info._node->_parentAuxNodeIndex == child_node_info._cliqueIndexInNode);
			MY_ASSERT(child_node_info._node->_parentPixNodeIndex == current_rec_aux_info._pixNodeIndexInAuxNode);
			//MY_ASSERT((new_rec_pix_info._tightEdgeIntersectionCode & (1 << child_index)) != 0);
			DEBUG_STATEMENT(PRINTF("Resetting edge "); BKGCDump::PrintNode(child_node_info._node)); 
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(current_rec_pix_node));
			DEBUG_STATEMENT(PRINTF(" to "); BKGCDump::PrintNode(child_node_info._node)); 
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(new_rec_pix_info._node); PRINTF("\n"));
			child_node_info._node->_parentPixNodeIndex = new_rec_index;
			new_rec_aux_info._childrenCode |= (1 << child_index);
		}
	}
};
