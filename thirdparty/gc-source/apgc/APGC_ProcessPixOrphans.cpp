#include "StdAfx.h"
#include "apgc/APGC.h"
#include "apgc/APGCDump.h"

void APGC::ProcessSinkOrphan(APGCPixNode* sending_node)
{
	//DEBUG_COND(ITER_NUM == 6);

	for(int aux_node_index=0; aux_node_index<sending_node->_numAuxNodes; ++aux_node_index)
	{
		APGCAuxNodeInfo& sending_node_aux_info = sending_node->_auxNodeInfo[aux_node_index];
		int sending_node_index = sending_node_aux_info._pixNodeIndexInAuxNode;
		APGCAuxNode* aux_node = sending_node_aux_info._auxNode;
		APGCPixNodeInfo& sending_node_pix_info = aux_node->_nodeInfo[sending_node_index];

		if(sending_node_pix_info._numTightConstraints == 0)
		{
			if(aux_node->_tree == TREE_SINK)
			{
				if(!aux_node->HasPathToRoot())
					continue;
				if(aux_node->_distanceLabel == (sending_node->_distanceLabel-1))
				{
					sending_node->_parentAuxNodeIndex = aux_node_index;
					sending_node->_parentPixNodeIndex = -1;
					DEBUG_STATEMENT(PRINTF("Creating PA edge "); APGCDump::PrintEdge(sending_node,aux_node); PRINTLN); 
					aux_node->_childrenCode |= (1 << sending_node_index);
					DeclareActive(sending_node);
					return;
				}
			}
			continue;
		}

		CODE_TYPE union_constraint_edges = sending_node_pix_info._tightEdgeUnionCode & INDEX_CODE_ZERO[sending_node_index];
		CODE_TYPE node_seeking_label_from = union_constraint_edges & aux_node->_sinkLabeledNodesCode;
		while(node_seeking_label_from != 0)
		{
			int rec_node_index = FIRST_NON_ZERO(node_seeking_label_from);
			node_seeking_label_from &= INDEX_CODE_ZERO[rec_node_index];
			APGCPixNodeInfo& rec_node_pix_info = aux_node->_nodeInfo[rec_node_index]; 
			MY_ASSERT(rec_node_pix_info._numTightConstraints > 0);
			if((rec_node_pix_info._tightEdgeIntersectionCode & (1 << sending_node_index)) == 0)
				continue;
			APGCPixNode* rec_node = rec_node_pix_info._node; 

			//DEBUG_COND(sending_node->_x == 4 && sending_node->_y == 0 && sending_node->_nodeIndex == 2 && rec_node->_x == 4 && rec_node->_y == 1 && rec_node->_nodeIndex == 2 && ITER_NUM == 56);

			if(rec_node->_parentAuxNodeIndex == rec_node_pix_info._cliqueIndexInNode)
				continue;
			if(!rec_node->HasPathToRoot())
				continue;
			MY_ASSERT(rec_node->_distanceLabel >= (sending_node->_distanceLabel-2));
			if(rec_node->_distanceLabel == (sending_node->_distanceLabel-2))
			{
				sending_node->_parentAuxNodeIndex = aux_node_index;
				sending_node->_parentPixNodeIndex = rec_node_index;
				APGCAuxNodeInfo& rec_node_aux_info = rec_node->_auxNodeInfo[rec_node_pix_info._cliqueIndexInNode];
				DEBUG_STATEMENT(PRINTF("Creating PP edge "); APGCDump::PrintEdge(aux_node, sending_node,rec_node); PRINTLN); 
				rec_node_aux_info._childrenCode |= (1 << sending_node_index);				
				DeclareActive(sending_node);
				return;
			}
			DeclareActive(rec_node);
		}
	}

	DeclareFree(sending_node);
}

void APGC::ProcessSrcOrphan(APGCPixNode* rec_node)
{
	for(int aux_node_index=0; aux_node_index<rec_node->_numAuxNodes; ++aux_node_index)
	{
		APGCAuxNodeInfo& rec_node_aux_info = rec_node->_auxNodeInfo[aux_node_index];
		int rec_node_index = rec_node_aux_info._pixNodeIndexInAuxNode;
		APGCAuxNode* aux_node = rec_node_aux_info._auxNode;
		APGCPixNodeInfo& rec_node_pix_info = aux_node->_nodeInfo[rec_node_index];
		if(rec_node_pix_info._numTightConstraints == 0)
		{
			if(aux_node->_tree == TREE_SRC)
			{
				if(!aux_node->HasPathToRoot())
					continue;
				if(aux_node->_distanceLabel == (rec_node->_distanceLabel-1))
				{
					rec_node->_parentAuxNodeIndex = aux_node_index;
					rec_node->_parentPixNodeIndex = -1;
					DEBUG_STATEMENT(PRINTF("Creating AP edge "); APGCDump::PrintEdge(aux_node,rec_node); PRINTLN); 
					aux_node->_childrenCode |= (1 << rec_node_index);
					DeclareActive(rec_node);
					return;
				}
			}
			continue;
		}

		CODE_TYPE me_or_my_children = rec_node_aux_info._childrenCode | INDEX_CODE_ONE[rec_node_index];
		CODE_TYPE nodes_i_can_seek_label = rec_node_pix_info._tightEdgeIntersectionCode & (me_or_my_children ^ CODE_ALL_ONES_WORD);
		nodes_i_can_seek_label &= aux_node->_srcLabeledNodesCode;
		while(nodes_i_can_seek_label != 0)
		{
			int send_node_index = FIRST_NON_ZERO(nodes_i_can_seek_label);
			nodes_i_can_seek_label &= INDEX_CODE_ZERO[send_node_index];
			APGCPixNodeInfo& send_node_pix_info = aux_node->_nodeInfo[send_node_index]; 
			APGCPixNode* send_node = send_node_pix_info._node; 

			//DEBUG_COND(send_node->_x == 2 && send_node->_y == 2 && send_node->_nodeIndex == 1 && rec_node->_x == 1 && rec_node->_y == 2 && rec_node->_nodeIndex == 1 && ITER_NUM == 19);

			if(send_node->_parentAuxNodeIndex == send_node_pix_info._cliqueIndexInNode)
				continue;

			if(!send_node->HasPathToRoot())
				continue;

			MY_ASSERT(send_node->_distanceLabel >= (rec_node->_distanceLabel-2));
			if(send_node->_distanceLabel == (rec_node->_distanceLabel-2))
			{
				rec_node->_parentAuxNodeIndex = aux_node_index;
				rec_node->_parentPixNodeIndex = send_node_index;
				APGCAuxNodeInfo& send_node_aux_info = send_node->_auxNodeInfo[send_node_pix_info._cliqueIndexInNode];
				DEBUG_STATEMENT(PRINTF("Creating PP edge "); APGCDump::PrintEdge(aux_node, send_node, rec_node); PRINTLN);
				send_node_aux_info._childrenCode |= (1 << rec_node_index);
				DeclareActive(rec_node);
				return;
			}
			DeclareActive(send_node);
		}
	}

	DeclareFree(rec_node);
}

void APGC::DeclareFree(APGCPixNode* node)
{
	MY_ASSERT(node->_tree != TREE_NONE);
	for(int i=0; i<node->_numAuxNodes; ++i)
	{
		APGCAuxNodeInfo& aux_info = node->_auxNodeInfo[i];
		APGCAuxNode* aux_node = aux_info._auxNode;
		if(node->_tree == TREE_SINK)
			aux_node->_sinkLabeledNodesCode &= INDEX_CODE_ZERO[aux_info._pixNodeIndexInAuxNode];
		else
			aux_node->_srcLabeledNodesCode &= INDEX_CODE_ZERO[aux_info._pixNodeIndexInAuxNode];
	}
	DeclareChildrenOrphan(node);
	node->_tree = TREE_NONE;
	node->_distanceLabel = DIST_DEFAULT;

	for(int i=0; i<node->_numAuxNodes; ++i)
	{
		APGCAuxNodeInfo& my_aux_info = node->_auxNodeInfo[i];
		APGCAuxNode* aux_node = my_aux_info._auxNode;
		if(aux_node->_tree != TREE_NONE)
			DeclareActive(aux_node);
		for(int i=0; i<GG_CLIQUE_SIZE; ++i)
		{
			APGCPixNode* pix_node = aux_node->_nodeInfo[i]._node;
			if(pix_node->_tree == TREE_NONE)
				continue;
			DeclareActive(pix_node);
		}
	}
}

void APGC::RemoveLinkToParentAndPushInOraphanList(APGCPixNode* pix_node)
{
	DEBUG_STATEMENT(PRINTF("Declaring orphan: "); APGCDump::PrintNode(pix_node); PRINTLN);
	//DEBUG_COND(ITER_NUM == 30 && pix_node->_x == 2 && pix_node->_y == 2);
	if(pix_node->_parentAuxNodeIndex != -1)
	{
		if(pix_node->_parentPixNodeIndex == -1)
		{
			APGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[pix_node->_parentAuxNodeIndex];
			APGCAuxNode* aux_node = aux_info._auxNode;
			aux_node->_childrenCode &= INDEX_CODE_ZERO[aux_info._pixNodeIndexInAuxNode];
		}else
		{
			APGCAuxNodeInfo& my_aux_info = pix_node->_auxNodeInfo[pix_node->_parentAuxNodeIndex];
			APGCPixNodeInfo& parent_node_info = my_aux_info._auxNode->_nodeInfo[pix_node->_parentPixNodeIndex];
			APGCAuxNodeInfo& parent_aux_info = parent_node_info._node->_auxNodeInfo[parent_node_info._cliqueIndexInNode];
			parent_aux_info._childrenCode &= INDEX_CODE_ZERO[my_aux_info._pixNodeIndexInAuxNode];
		}
		pix_node->_parentAuxNodeIndex = -1;
		pix_node->_parentPixNodeIndex = -1;
	}else{
		MY_ASSERT(pix_node->_distanceLabel == TERM_DIST);
	}

	_orphanNodes[pix_node->_distanceLabel].push_back(pix_node);
}

void APGC::DeclareChildrenOrphan(APGCPixNode* rec_node)
{
	for(int i=0; i<rec_node->_numAuxNodes; ++i)
	{
		APGCAuxNodeInfo& rec_node_aux_info = rec_node->_auxNodeInfo[i];
		int rec_node_index = rec_node_aux_info._pixNodeIndexInAuxNode;
		APGCAuxNode* aux_node = rec_node_aux_info._auxNode;
		if(aux_node->_parentPixNodeIndex == rec_node_index)
		{
			MY_ASSERT(aux_node->_tree == rec_node->_tree);
			RemoveLinkToParentAndPushInOraphanList(aux_node);
		}
		CODE_TYPE temp_children = rec_node_aux_info._childrenCode;
		while (temp_children != 0)
		{
			int i = FIRST_NON_ZERO(temp_children);
			temp_children &= (1 << i) ^ CODE_ALL_ONES_WORD;
			APGCPixNodeInfo& node_info = aux_node->_nodeInfo[i];
			RemoveLinkToParentAndPushInOraphanList(node_info._node);
		}
	}
}