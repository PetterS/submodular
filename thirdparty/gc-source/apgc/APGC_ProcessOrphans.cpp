#include "StdAfx.h"
#include "apgc/APGC.h"
#include "apgc/APGCDump.h"

void APGC::ProcessTightConstraints()
{
	//DEBUG_COND(ITER_NUM == 16);

	APGCPixNode* join_edge_rec_node = NULL;
	if(_joinEdge->_receivingNodeIndex != -1)
		join_edge_rec_node = _joinEdge->_auxNode->_nodeInfo[_joinEdge->_receivingNodeIndex]._node;

	while(!TIGHT_CONSTRAINTS.empty())
	{
		std::pair<APGCAuxNode*, int> top = TIGHT_CONSTRAINTS.front();
		TIGHT_CONSTRAINTS.pop_front();
		APGCAuxNode* aux_node = top.first;
		int constraint_index = top.second;
		ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];

		//DEBUG_ITER(10);

		DEBUG_STATEMENT(PRINTF("processing tight constraint %d-%s, of ", constraint_index, BINARY(constraint_info._labeling)); APGCDump::PrintNode(aux_node); PRINTLN);
		CODE_TYPE temp_labeling = constraint_info._labeling;
		while(temp_labeling != 0)
		{
			uchar pix_node_index = FIRST_NON_ZERO(temp_labeling);
			temp_labeling &= (1 << pix_node_index) ^ CODE_ALL_ONES_WORD;
			APGCPixNodeInfo& pix_node_info = aux_node->_nodeInfo[pix_node_index];
			APGCPixNode* pix_node = pix_node_info._node;
			if(pix_node->_tree == TREE_SINK)
			{
				if(pix_node == join_edge_rec_node){
					if(_joinEdge->_sendingNodeIndex == -1 || (pix_node_info._tightEdgeIntersectionCode & (1 << _joinEdge->_sendingNodeIndex))== 0)
						_joinEdge->_auxNode = NULL;
				}

				if(aux_node->_parentPixNodeIndex == pix_node_index){
					MY_ASSERT(aux_node->_tree == TREE_SINK);
					RemoveLinkToParentAndPushInOraphanList(aux_node);
				}

				APGCAuxNodeInfo& rec_node_aux_info = pix_node->_auxNodeInfo[pix_node_info._cliqueIndexInNode];
				CODE_TYPE children = rec_node_aux_info._childrenCode;
				while(children != 0)
				{
					int send_node_index = FIRST_NON_ZERO(children);
					children &= (1 << send_node_index) ^ CODE_ALL_ONES_WORD;
					if((pix_node_info._tightEdgeIntersectionCode & (1 << send_node_index)) != 0)
						continue;
					rec_node_aux_info._childrenCode &= (1 << send_node_index) ^ CODE_ALL_ONES_WORD;
					APGCPixNode* sending_node = aux_node->_nodeInfo[send_node_index]._node;
					RemoveLinkToParentAndPushInOraphanList(sending_node);
				}
			}else
			{
				if(pix_node->_parentAuxNodeIndex != pix_node_info._cliqueIndexInNode)
					continue;				
				if(pix_node->_parentPixNodeIndex != -1 && 
					(pix_node_info._tightEdgeIntersectionCode & (1 << pix_node->_parentPixNodeIndex)) != 0)
					continue;
				if(pix_node->_parentAuxNodeIndex != -1){
					if(pix_node->_parentPixNodeIndex != -1){
						APGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[pix_node->_parentPixNodeIndex];
						APGCAuxNodeInfo& parent_aux_info = parent_node_info._node->_auxNodeInfo[parent_node_info._cliqueIndexInNode];
						parent_aux_info._childrenCode &= (1 << pix_node_index) ^ CODE_ALL_ONES_WORD;
					}else{
						aux_node->_childrenCode &= (1 << pix_node_index) ^ CODE_ALL_ONES_WORD;
					}
				}
				RemoveLinkToParentAndPushInOraphanList(pix_node);
			}
		}
	}
}

void APGC::ProcessOrphans()
{
	for(int i=0; i<=_maxAssignedDistLabel; ++i)
	{
		APGCNodeList& ol = _orphanNodes[i];
		while(!ol.empty())
		{
			APGCNode* node = ol.front();
			DEBUG_STATEMENT(PRINTF("\nProcessing orphan node"); APGCDump::PrintNode(node); PRINTLN);
			ol.pop_front();
			++TIME;
			if(node->_nodeType == NODE_TYPE_AUX){
				if(node->_tree == TREE_SRC)
					ProcessSrcOrphan((APGCAuxNode*)node);
				else
					ProcessSinkOrphan((APGCAuxNode*)node);
			}else{
				if(node->_tree == TREE_SRC)
					ProcessSrcOrphan((APGCPixNode*)node);
				else
					ProcessSinkOrphan((APGCPixNode*)node);
			}
		}
	}
}
