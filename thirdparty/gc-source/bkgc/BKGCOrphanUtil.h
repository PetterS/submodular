#pragma once

#include "bkgc/BKGCDataStructures.h"
#include "bkgc/BKGCDebug.h"

class BKGCOrphanUtil
{
public:
	BKGCOrphanUtil(void);
	virtual ~BKGCOrphanUtil(void);

protected:

	//TODO: I should keep a flag in the node which says the node was searched and keep it in some vector
	//the nodes in the subtree which has not searched for path and are being orphan should try to keep make
	//actives..rest should not do it. Upon finishing all the nodes that were searched (and kept in a vector) 
	//should reset their searched flags

	//TODO: while declaring tree orphan, there should be option that if I want to change the tree, I should change it and simply tell my children
	//that the tree has changed rather than declaring them orphan

	inline void DeclareSrcPixNodeOrphan(BKGCPixNode* pix_node)
	{
		//DEBUG_COND(ITER_NUM == 30 && pix_node->_x == 2 && pix_node->_y == 2);
		for(int i=0; i<pix_node->_numAuxNodes; ++i)
		{
			BKGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[i];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			aux_node->_srcLabeledNodesCode &= (1 << aux_info._pixNodeIndexInAuxNode) ^ CODE_ALL_ONES_WORD;
		}
		pix_node->_parentAuxNodeIndex = -1;
		pix_node->_parentPixNodeIndex = -1;
		pix_node->_tree = TREE_NONE;
	}

	inline void DeclareSinkPixNodeOrphan(BKGCPixNode* pix_node)
	{
		//DEBUG_COND(ITER_NUM == 30 && pix_node->_x == 2 && pix_node->_y == 2);
		for(int i=0; i<pix_node->_numAuxNodes; ++i)
		{
			BKGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[i];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			aux_node->_sinkLabeledNodesCode &= (1 << aux_info._pixNodeIndexInAuxNode) ^ CODE_ALL_ONES_WORD;
		}
		pix_node->_parentAuxNodeIndex = -1;
		pix_node->_parentPixNodeIndex = -1;
		pix_node->_tree = TREE_NONE;
	}

	void DeclareTreeOrphan(uchar tree_type, BKGCNode* root, std::list<BKGCNode*> *activeNodesList);

protected:
	std::list<BKGCNode*> *_currentNodeList, *_nextNodeList;
};
