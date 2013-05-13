#include "StdAfx.h"
#include "bkgc/BKGC.h"
#include "bkgc/BKGCDataStructures.h"
#include "bkgc/BKGCDump.h"

void BKGC::ProcessActives()
{
	while(!_activeNodes->empty())
	{
		BKGCNode* graph_node = _activeNodes->front();
		//DEBUG_COND(ITER_NUM == 4);
		if(graph_node->_tree != TREE_NONE)
		{
			if(graph_node->_tree == TREE_SINK)
			{
				if(graph_node->_nodeType == NODE_TYPE_AUX)
					ProcessSinkAuxActive((BKGCAuxNode*)graph_node);
				else
					ProcessSinkPixActive((BKGCPixNode*)graph_node);
			}else 
			{
				if(graph_node->_nodeType == NODE_TYPE_AUX)
					ProcessSrcAuxActive((BKGCAuxNode*)graph_node);
				else
					ProcessSrcPixActive((BKGCPixNode*)graph_node);
			}
			if(_joinEdge->_auxNode != NULL){
				DEBUG_STATEMENT(PRINTF("found joining edge "); BKGCDump::PrintEdge(*_joinEdge); PRINTF("\n"));
				return;
			}
		}
		graph_node->_active = false;
		_activeNodes->pop_front();
	}
}

void BKGC::ProcessSinkPixActive(BKGCPixNode* node)
{
	MY_ASSERT(node->_active);
	DEBUG_STATEMENT(PRINTF("Processing sink active node "); BKGCDump::PrintNode(node); PRINTF("\n"));
	MY_ASSERT(node->_excess < -EPSILON || node->_parentAuxNodeIndex != -1);
	//DEBUG_COND(ITER_NUM == 34);
	for(int i=0; i<node->_numAuxNodes; ++i)
	{
		if(i == node->_parentAuxNodeIndex)
			continue;
		BKGCAuxNodeInfo& receiving_node_aux_info = node->_auxNodeInfo[i];
		BKGCAuxNode* aux_node = receiving_node_aux_info._auxNode;
		int receiving_node_index = receiving_node_aux_info._pixNodeIndexInAuxNode;
		BKGCPixNodeInfo& receiving_node_info = aux_node->_nodeInfo[receiving_node_index];
		if(receiving_node_info._numTightConstraints == 0)
		{
			if(aux_node->_parentIndex == -1)
			{
				aux_node->_parentIndex = receiving_node_aux_info._pixNodeIndexInAuxNode;
				aux_node->_tree = TREE_SINK;
				DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintNode(aux_node)); 
				DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(receiving_node_info._node); PRINTF("\n"));
				if(!aux_node->_active){
					aux_node->_active = true;
					_activeNodes->push_back(aux_node);
				}
			}else
			{
				if(aux_node->_tree == TREE_SRC)	
				{
					_joinEdge->_auxNode = aux_node;
					_joinEdge->_sendingNodeIndex = -1;
					_joinEdge->_receivingNodeIndex = receiving_node_index;
					return;
				}
			}
		}else
		{
			int possible_join_nodes = aux_node->_srcLabeledNodesCode & receiving_node_info._tightEdgeIntersectionCode;
			if(possible_join_nodes != 0)
			{
				_joinEdge->_auxNode = aux_node;
				_joinEdge->_sendingNodeIndex = FIRST_NON_ZERO_WORD[possible_join_nodes];
				_joinEdge->_receivingNodeIndex = receiving_node_index;
				return;
			}

			if(aux_node->_parentIndex != -1 && aux_node->_tree == TREE_SINK)
				continue;

			int unlabeled_nodes = aux_node->_sinkLabeledNodesCode ^ CODE_ALL_ONES_WORD;
			int nodes_i_can_give_label = unlabeled_nodes & receiving_node_info._tightEdgeIntersectionCode;
			while(nodes_i_can_give_label != 0)
			{
				uchar j = FIRST_NON_ZERO_WORD[nodes_i_can_give_label];
				nodes_i_can_give_label &= (1 << j) ^ CODE_ALL_ONES_WORD;
				BKGCPixNodeInfo& sending_node_info = aux_node->_nodeInfo[j];
				BKGCPixNode* sending_node = sending_node_info._node;
				MY_ASSERT(sending_node->_parentAuxNodeIndex == -1);
				MY_ASSERT(sending_node->_tree == TREE_NONE);
				sending_node->_parentAuxNodeIndex = sending_node_info._cliqueIndexInNode;
				sending_node->_parentPixNodeIndex = receiving_node_index;
				receiving_node_aux_info._childrenCode |= (1 << j);
				DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintNode(sending_node)); 
				DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(receiving_node_info._node); PRINTF("\n"));
				for(int i=0; i<sending_node->_numAuxNodes; ++i){
					BKGCAuxNodeInfo& aux_info = sending_node->_auxNodeInfo[i];
					aux_info._auxNode->_sinkLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
				}
				sending_node->_tree = TREE_SINK;
				if(!sending_node->_active){
					sending_node->_active = true;
					_activeNodes->push_back(sending_node);
				}
			}
		}
	}
}

void BKGC::ProcessSinkAuxActive(BKGCAuxNode* aux_node)
{
	MY_ASSERT(aux_node->_active);
	DEBUG_STATEMENT(PRINTF("Processing sink active node "); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
	//DEBUG_COND(ITER_NUM == 9);
	MY_ASSERT(aux_node->_parentIndex != -1);
	if(aux_node->_srcLabeledNodesCode != 0)
	{
		_joinEdge->_auxNode = aux_node;
		_joinEdge->_receivingNodeIndex = -1;
		_joinEdge->_sendingNodeIndex = FIRST_NON_ZERO_WORD[aux_node->_srcLabeledNodesCode];
		return;
	}

	int unlabeled_nodes = aux_node->_sinkLabeledNodesCode ^ CODE_ALL_ONES_WORD;
	while(unlabeled_nodes != 0)
	{
		uchar j = FIRST_NON_ZERO_WORD[unlabeled_nodes];
		unlabeled_nodes &= (1 << j) ^ CODE_ALL_ONES_WORD;
		BKGCPixNodeInfo& sending_node_info = aux_node->_nodeInfo[j];
		BKGCPixNode* sending_node = sending_node_info._node;
		MY_ASSERT(sending_node->_parentAuxNodeIndex == -1);
		sending_node->_parentAuxNodeIndex = sending_node_info._cliqueIndexInNode;
		aux_node->_childrenCode |= (1 << j);
		DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintNode(sending_node)); 
		DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
		for(int i=0; i<sending_node->_numAuxNodes; ++i){
			BKGCAuxNodeInfo& aux_info = sending_node->_auxNodeInfo[i];
			aux_info._auxNode->_sinkLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
		}
		sending_node->_tree = TREE_SINK;
		if(!sending_node->_active){
			sending_node->_active = true;
			_activeNodes->push_back(sending_node);
		}
	}
}

void BKGC::ProcessSrcPixActive(BKGCPixNode* node)
{
	MY_ASSERT(node->_active);
	DEBUG_STATEMENT(PRINTF("Processing src active node "); BKGCDump::PrintNode(node); PRINTF("\n"));
	MY_ASSERT(node->_excess > 0 || node->_parentAuxNodeIndex != -1);
	for(int i=0; i<node->_numAuxNodes; ++i)
	{
		//DEBUG_COND(ITER_NUM == 6);// && node->_x == 1 && node->_y == 0);

		if(i == node->_parentAuxNodeIndex)
			continue;
		BKGCAuxNodeInfo& sending_node_aux_info = node->_auxNodeInfo[i];
		BKGCAuxNode* aux_node = sending_node_aux_info._auxNode;
		int sending_node_index = sending_node_aux_info._pixNodeIndexInAuxNode;
		BKGCPixNodeInfo& sending_node_info = aux_node->_nodeInfo[sending_node_index];

		if(aux_node->_parentIndex == -1)
		{
			aux_node->_parentIndex = sending_node_index;
			aux_node->_tree = TREE_SRC;
			DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintNode(sending_node_info._node)); 
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
			if(!aux_node->_active){
				aux_node->_active = true;
				_activeNodes->push_back(aux_node);
			}
		}else
		{
			if(aux_node->_tree == TREE_SINK)	
			{
				_joinEdge->_auxNode = aux_node;
				_joinEdge->_sendingNodeIndex = sending_node_index;
				_joinEdge->_receivingNodeIndex = -1;
				return;
			}
		}

		if(sending_node_info._numTightConstraints > 0)
		{
			int possible_join_nodes = aux_node->_sinkLabeledNodesCode & sending_node_info._tightEdgeUnionCode;
			while(possible_join_nodes != 0)
			{
				int j = FIRST_NON_ZERO_WORD[possible_join_nodes];
				possible_join_nodes &= (1 << j) ^ CODE_ALL_ONES_WORD;
				BKGCPixNodeInfo& rec_node_info = aux_node->_nodeInfo[j];
				if((rec_node_info._tightEdgeIntersectionCode & (1 << sending_node_index)) != 0)
				{
					_joinEdge->_auxNode = aux_node;
					_joinEdge->_sendingNodeIndex = sending_node_index;
					_joinEdge->_receivingNodeIndex = j;
					return;
				}
			}

			int unlabeled_not_me_nodes = (aux_node->_sinkLabeledNodesCode | aux_node->_srcLabeledNodesCode 
									| (1 << sending_node_index)) ^ CODE_ALL_ONES_WORD;
			int nodes_i_can_give_label = unlabeled_not_me_nodes & sending_node_info._tightEdgeUnionCode;
			while(nodes_i_can_give_label != 0)
			{
				int j = FIRST_NON_ZERO_WORD[nodes_i_can_give_label];
				nodes_i_can_give_label &= (1 << j) ^ CODE_ALL_ONES_WORD;
				BKGCPixNodeInfo& rec_node_info = aux_node->_nodeInfo[j];
				if((rec_node_info._tightEdgeIntersectionCode & (1 << sending_node_index)) == 0)
					continue;
				BKGCPixNode* rec_node = rec_node_info._node;
				MY_ASSERT(rec_node->_tree == TREE_NONE);
				rec_node->_parentAuxNodeIndex = rec_node_info._cliqueIndexInNode;
				rec_node->_parentPixNodeIndex = sending_node_index;
				sending_node_aux_info._childrenCode |= (1 << j);
				DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintNode(sending_node_info._node)); 
				DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(rec_node); PRINTF("\n"));
				for(int k=0; k<rec_node->_numAuxNodes; ++k){
					BKGCAuxNodeInfo& aux_info = rec_node->_auxNodeInfo[k];
					aux_info._auxNode->_srcLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
				}
				rec_node->_tree = TREE_SRC;
				if(!rec_node->_active){
					rec_node->_active = true;
					_activeNodes->push_back(rec_node);
				}
			}
		}
	}
}

void BKGC::ProcessSrcAuxActive(BKGCAuxNode* aux_node)
{
	MY_ASSERT(aux_node->_active);
	DEBUG_STATEMENT(PRINTF("Processing src active node "); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
	//DEBUG_COND(ITER_NUM == 9);
	MY_ASSERT(aux_node->_parentIndex != -1);
	int sink_labeled_unsaturated_nodes = aux_node->_sinkLabeledNodesCode & (aux_node->_saturatedEdgeCode ^ CODE_ALL_ONES_WORD);
	if(sink_labeled_unsaturated_nodes != 0)
	{
		_joinEdge->_auxNode = aux_node;
		_joinEdge->_sendingNodeIndex = -1;
		_joinEdge->_receivingNodeIndex = FIRST_NON_ZERO_WORD[sink_labeled_unsaturated_nodes];
		return;
	}

	int unlabeled_unsaturated_nodes = (aux_node->_sinkLabeledNodesCode | aux_node->_srcLabeledNodesCode | aux_node->_saturatedEdgeCode) ^ CODE_ALL_ONES_WORD;
	while(unlabeled_unsaturated_nodes != 0)
	{
		int j = FIRST_NON_ZERO_WORD[unlabeled_unsaturated_nodes];
		unlabeled_unsaturated_nodes &= (1 << j) ^ CODE_ALL_ONES_WORD;
		BKGCPixNodeInfo& rec_node_info = aux_node->_nodeInfo[j];
		BKGCPixNode* rec_node = rec_node_info._node;
		rec_node->_parentAuxNodeIndex = rec_node_info._cliqueIndexInNode;
		aux_node->_childrenCode |= (1 << j);
		DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintNode(aux_node)); 
		DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(rec_node); PRINTF("\n"));
		for(int i=0; i<rec_node->_numAuxNodes; ++i){
			BKGCAuxNodeInfo& aux_info = rec_node->_auxNodeInfo[i];
			aux_info._auxNode->_srcLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
		}
		rec_node->_tree = TREE_SRC;
		if(!rec_node->_active){
			rec_node->_active = true;
			_activeNodes->push_back(rec_node);
		}
	}
}
