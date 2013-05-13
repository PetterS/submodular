#include "StdAfx.h"
#include "apgc/APGC.h"
#include "apgc/APGCDump.h"

void APGC::RemoveLinkToParentAndPushInOraphanList(APGCAuxNode* aux_node)
{
	DEBUG_STATEMENT(PRINTF("Declaring orphan: "); APGCDump::PrintNode(aux_node); PRINTLN);
	//DEBUG_COND(ITER_NUM == 30 && pix_node->_x == 2 && pix_node->_y == 2);
	aux_node->_parentPixNodeIndex = -1;
	_orphanNodes[aux_node->_distanceLabel].push_back(aux_node);
}

void APGC::DeclareChildrenOrphan(APGCAuxNode* aux_node)
{
	CODE_TYPE temp_children = aux_node->_childrenCode;
	while (temp_children != 0)
	{
		int i = FIRST_NON_ZERO(temp_children);
		temp_children &= (1 << i) ^ CODE_ALL_ONES_WORD;
		APGCPixNodeInfo& node_info = aux_node->_nodeInfo[i];
		RemoveLinkToParentAndPushInOraphanList(node_info._node);
	}
}

void APGC::DeclareFree(APGCAuxNode* aux_node)
{
	MY_ASSERT(aux_node->_tree != TREE_NONE);
	DeclareChildrenOrphan(aux_node);
	aux_node->_tree = TREE_NONE;
	aux_node->_distanceLabel = DIST_DEFAULT;

	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		APGCPixNode* pix_node = aux_node->_nodeInfo[i]._node;
		if(pix_node->_tree == TREE_NONE)
			continue;
		DeclareActive(pix_node);
	}
}

void APGC::ProcessSrcOrphan(APGCAuxNode* aux_node)
{
	int parent_desired_label = aux_node->_distanceLabel-1;

	int src_labeled_nodes = aux_node->_srcLabeledNodesCode;
	while(src_labeled_nodes != 0)
	{
		int send_node_index = FIRST_NON_ZERO(src_labeled_nodes);
		src_labeled_nodes &= INDEX_CODE_ZERO[send_node_index];
		APGCPixNodeInfo& send_node_pix_info = aux_node->_nodeInfo[send_node_index];
		APGCPixNode* send_node = send_node_pix_info._node;
		if(!send_node->HasPathToRoot())
			continue;
		if(send_node->_distanceLabel == parent_desired_label)
		{
			DEBUG_STATEMENT(PRINTF("Creating PA edge "); APGCDump::PrintEdge(send_node, aux_node); PRINTLN);
			aux_node->_parentPixNodeIndex = send_node_index;
			DeclareActive(aux_node);
			return;
		}
		DeclareActive(send_node);
	}

	DeclareFree(aux_node);
}

void APGC::ProcessSinkOrphan(APGCAuxNode* aux_node)
{
	int parent_desired_label = aux_node->_distanceLabel-1;

	int possible_parents = aux_node->_sinkLabeledNodesCode & (aux_node->_saturatedEdgeCode ^ CODE_ALL_ONES_WORD);
	while(possible_parents != 0)
	{
		int rec_node_index = FIRST_NON_ZERO(possible_parents);
		possible_parents &= INDEX_CODE_ZERO[rec_node_index];
		APGCPixNode* rec_node = aux_node->_nodeInfo[rec_node_index]._node;
		if(!rec_node->HasPathToRoot())
			continue;
		if(rec_node->_distanceLabel == parent_desired_label)
		{
			DEBUG_STATEMENT(PRINTF("Creating AP edge "); APGCDump::PrintEdge(aux_node,rec_node); PRINTLN);
			aux_node->_parentPixNodeIndex = rec_node_index;
			DeclareActive(aux_node);
			return;
		}
		DeclareActive(rec_node);
	}

	DeclareFree(aux_node);
}