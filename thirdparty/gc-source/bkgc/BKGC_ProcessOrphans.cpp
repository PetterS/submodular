#include "StdAfx.h"
#include "bkgc/BKGC.h"
#include "bkgc/BKGCSinkOrphanRootProcessor.h"
#include "bkgc/BKGCSrcOrphanRootProcessor.h"

void BKGC::ProcessTightConstraints()
{
//	DEBUG_COND(ITER_NUM == 119);
	if(EQUAL(_currentPathSinkNode->_excess,0)){
		DEBUG_STATEMENT(BKGCDump::PrintNode(_currentPathSinkNode); PRINTF(" became saturated\n"));
		_sinkTreeOrphanRootProcessor->FindPath(_currentPathSinkNode, _activeNodes);
	}

	if(_currentPathSrcNode->_excess <= EPSILON){
		DEBUG_STATEMENT(BKGCDump::PrintNode(_currentPathSrcNode); PRINTF(" became saturated\n"));
		_srcTreeOrphanRootProcessor->FindPath(_currentPathSrcNode, _activeNodes);
	}

	//DEBUG_COND(ITER_NUM == 16);

	BKGCPixNode* join_edge_rec_node = NULL;
	if(_joinEdge->_receivingNodeIndex != -1)
		join_edge_rec_node = _joinEdge->_auxNode->_nodeInfo[_joinEdge->_receivingNodeIndex]._node;

	while(!_tightConstraints.empty())
	{
		std::pair<BKGCAuxNode*, int> top = _tightConstraints.front();
		_tightConstraints.pop_front();
		BKGCAuxNode* aux_node = top.first;
		int index = top.second;
		DEBUG_STATEMENT(PRINTF("processing tight constraint %d of ",index); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
		int temp_labeling = CONSTRAINT_LABELING(index);
		while(temp_labeling != 0)
		{
			uchar pix_node_index = FIRST_NON_ZERO_WORD[temp_labeling];
			temp_labeling &= (1 << pix_node_index) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& pix_node_info = aux_node->_nodeInfo[pix_node_index];
			BKGCPixNode* pix_node = pix_node_info._node;
			if(pix_node->_tree == TREE_SINK)
			{
				if(pix_node == join_edge_rec_node){
					if(_joinEdge->_sendingNodeIndex == -1 || (pix_node_info._tightEdgeIntersectionCode & (1 << _joinEdge->_sendingNodeIndex))== 0)
						_joinEdge->_auxNode = NULL;
				}

				if(aux_node->_parentIndex == pix_node_index){
					MY_ASSERT(aux_node->_tree == TREE_SINK);
					_sinkTreeOrphanRootProcessor->FindPath(aux_node, _activeNodes);
					//continue;
				}

				BKGCAuxNodeInfo& rec_node_aux_info = pix_node->_auxNodeInfo[pix_node_info._cliqueIndexInNode];
				int children = rec_node_aux_info._childrenCode;
				while(children != 0)
				{
					int send_node_index = FIRST_NON_ZERO_WORD[children];
					children &= (1 << send_node_index) ^ CODE_ALL_ONES_WORD;
					if((pix_node_info._tightEdgeIntersectionCode & (1 << send_node_index)) != 0)
						continue;
					rec_node_aux_info._childrenCode &= (1 << send_node_index) ^ CODE_ALL_ONES_WORD;
					BKGCPixNode* sending_node = aux_node->_nodeInfo[send_node_index]._node;
					_sinkTreeOrphanRootProcessor->FindPath(sending_node, _activeNodes);
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
						BKGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[pix_node->_parentPixNodeIndex];
						BKGCAuxNodeInfo& parent_aux_info = parent_node_info._node->_auxNodeInfo[parent_node_info._cliqueIndexInNode];
						parent_aux_info._childrenCode &= (1 << pix_node_index) ^ CODE_ALL_ONES_WORD;
					}else{
						aux_node->_childrenCode &= (1 << pix_node_index) ^ CODE_ALL_ONES_WORD;
					}
				}
				_srcTreeOrphanRootProcessor->FindPath(pix_node, _activeNodes);
			}
		}
	}
}
