#include "StdAfx.h"
#include "bkgc/BKGCSrcOrphanRootProcessor.h"

void BKGCSrcOrphanRootProcessor::InvertTree(std::list<BKGCNode*> *activeNodesList, BKGCNode* root, BKGCNode* path_found_node, 
										int old_parent_aux_index, int old_parent_pix_index)
{
	DEBUG_STATEMENT(PRINTF("Found path - Inverting subtree rooted at "); BKGCDump::PrintNode(root); PRINTF("\n"));

	//DEBUG_COND(ITER_NUM == 17);

	if(!path_found_node->_active){
		activeNodesList->push_back(path_found_node);
		path_found_node->_active = true;
	}

	while(path_found_node != root)
	{
		if(path_found_node->_nodeType == NODE_TYPE_AUX)		
		{
			BKGCAuxNode* aux_node = (BKGCAuxNode*)path_found_node;
			aux_node->_childrenCode |= (1 << old_parent_pix_index);
			BKGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[old_parent_pix_index];
			BKGCPixNode* parent_node = parent_node_info._node;
			BKGCAuxNodeInfo& parent_aux_info = parent_node->_auxNodeInfo[parent_node_info._cliqueIndexInNode];
			MY_ASSERT(parent_node_info._numTightConstraints == 0);
			if(!parent_node->_active){
				activeNodesList->push_back(parent_node);
				parent_node->_active = true;
			}
			ResetOutgoingEdges(parent_aux_info, aux_node);
			DEBUG_STATEMENT(PRINTF("inverting edge "); BKGCDump::PrintEdge(aux_node,parent_node));
			DEBUG_STATEMENT(PRINTF(" to "); BKGCDump::PrintEdge(parent_node, aux_node); PRINTF("\n"));
			old_parent_aux_index = parent_node->_parentAuxNodeIndex;
			old_parent_pix_index = parent_node->_parentPixNodeIndex;
			parent_node->_parentAuxNodeIndex = parent_node_info._cliqueIndexInNode;
			parent_node->_parentPixNodeIndex = -1;
			path_found_node = parent_node;
		}else
		{
			BKGCPixNode* pix_node = (BKGCPixNode*)path_found_node;
			BKGCAuxNodeInfo& my_aux_node_info = pix_node->_auxNodeInfo[old_parent_aux_index];
			int my_index = my_aux_node_info._pixNodeIndexInAuxNode;
			BKGCAuxNode* aux_node = my_aux_node_info._auxNode;
			if(old_parent_pix_index == -1)
			{
				DEBUG_STATEMENT(PRINTF("inverting edge "); BKGCDump::PrintEdge(aux_node,pix_node));
				DEBUG_STATEMENT(PRINTF(" to "); BKGCDump::PrintEdge(pix_node,aux_node); PRINTF("\n"));
				if(!aux_node->_active){
					activeNodesList->push_back(aux_node);
					aux_node->_active = true;
				}
				aux_node->_childrenCode &= (1 << my_index) ^ CODE_ALL_ONES_WORD;
				//MY_ASSERT(pix_node->_parentAuxNodeIndex != old_parent_aux_index);
				old_parent_aux_index = -1;
				old_parent_pix_index = aux_node->_parentIndex;
				aux_node->_parentIndex = my_index;
				path_found_node = aux_node;
			}else
			{
				BKGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[old_parent_pix_index];
				BKGCPixNode* parent_node = parent_node_info._node;
				BKGCAuxNodeInfo& parent_aux_info = parent_node->_auxNodeInfo[parent_node_info._cliqueIndexInNode];
				DEBUG_STATEMENT(PRINTF("inverting edge "); BKGCDump::PrintEdge(parent_node,pix_node));
				DEBUG_STATEMENT(PRINTF(" to "); BKGCDump::PrintEdge(pix_node,parent_node); PRINTF("\n"));
				if(!parent_node->_active){
					activeNodesList->push_back(parent_node);
					parent_node->_active = true;
				}
				parent_aux_info._childrenCode &= (1 << my_index) ^ CODE_ALL_ONES_WORD;
				my_aux_node_info._childrenCode |= (1 << old_parent_pix_index);
				ResetOutgoingEdges(parent_aux_info, my_aux_node_info, my_index, aux_node);			
				old_parent_aux_index = parent_node->_parentAuxNodeIndex;
				old_parent_pix_index = parent_node->_parentPixNodeIndex;
				parent_node->_parentAuxNodeIndex = parent_node_info._cliqueIndexInNode;
				parent_node->_parentPixNodeIndex = my_index;
				path_found_node = parent_node;
			}
		}
	}
}

BKGCNode* BKGCSrcOrphanRootProcessor::FindPathInTheTree(BKGCNode* root, int& old_parent_aux_index, int& old_parent_pix_index)
{
	DEBUG_STATEMENT(PRINTF("Finding path for subtree rooted at "); BKGCDump::PrintNode(root); PRINTF("\n"));
	bool path_found = false;
	_currentNodeList->clear();
	_nextNodeList->clear();
	_currentNodeList->push_back(root);
	old_parent_aux_index = old_parent_pix_index = -1;
	while(true)
	{
		while(!_currentNodeList->empty())
		{
			BKGCNode* graph_node = _currentNodeList->front(); _currentNodeList->pop_front();
			if(graph_node->_nodeType == NODE_TYPE_AUX)
				path_found = FindPathThroughAuxNode((BKGCAuxNode*)graph_node, old_parent_pix_index);
			else
				path_found = FindPathThroughPixNode((BKGCPixNode*)graph_node, old_parent_aux_index, old_parent_pix_index);
			if(path_found)
				return graph_node;
		}
		if(_nextNodeList->empty())
			break;
		std::list<BKGCNode*> *temp = _currentNodeList;
		_currentNodeList = _nextNodeList;
		_nextNodeList = temp;
	}
	return NULL;
}

bool BKGCSrcOrphanRootProcessor::FindPathThroughPixNode(BKGCPixNode* rec_node, int& old_parent_aux_index, 
													int& old_parent_pix_index)
{
	DEBUG_STATEMENT(PRINTF("finding path through "); BKGCDump::PrintNode(rec_node); PRINTF("\n"));
	for(int aux_node_index=0; aux_node_index<rec_node->_numAuxNodes; ++aux_node_index)
	{
		BKGCAuxNodeInfo& rec_node_aux_info = rec_node->_auxNodeInfo[aux_node_index];
		int rec_node_index = rec_node_aux_info._pixNodeIndexInAuxNode;
		BKGCAuxNode* aux_node = rec_node_aux_info._auxNode;
		BKGCPixNodeInfo& rec_node_pix_info = aux_node->_nodeInfo[rec_node_index];
		if(rec_node_pix_info._numTightConstraints == 0 && aux_node->_tree == TREE_SRC
			&& aux_node->_parentIndex != rec_node_index && aux_node->_parentIndex != -1 
			&& aux_node->_nodeInfo[aux_node->_parentIndex]._node->HasPathToRoot())
		{
			old_parent_aux_index = rec_node->_parentAuxNodeIndex;
			old_parent_pix_index = rec_node->_parentPixNodeIndex;
			rec_node->_parentAuxNodeIndex = aux_node_index;
			rec_node->_parentPixNodeIndex = -1;
			aux_node->_childrenCode |= (1 << rec_node_index);
			DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintEdge(aux_node, rec_node); PRINTF("\n"));
			ResetOutgoingEdges(rec_node_aux_info, aux_node);
			return true;
		}

		if(rec_node_pix_info._numTightConstraints != 0)
		{
			int me_or_my_children = rec_node_aux_info._childrenCode | (1 << rec_node_index);
			int nodes_i_can_seek_label = rec_node_pix_info._tightEdgeIntersectionCode & 
										(me_or_my_children ^ CODE_ALL_ONES_WORD) &
										aux_node->_srcLabeledNodesCode;
			while(nodes_i_can_seek_label != 0)
			{
				int send_node_index = FIRST_NON_ZERO_WORD[nodes_i_can_seek_label];
				nodes_i_can_seek_label &= (1 << send_node_index) ^ CODE_ALL_ONES_WORD;
				BKGCPixNodeInfo& send_node_pix_info = aux_node->_nodeInfo[send_node_index]; 
				BKGCPixNode* send_node = send_node_pix_info._node; 
				if(send_node->_parentAuxNodeIndex == send_node_pix_info._cliqueIndexInNode)
					continue;
				if(!send_node->HasPathToRoot())
					continue;
				old_parent_aux_index = rec_node->_parentAuxNodeIndex;
				old_parent_pix_index = rec_node->_parentPixNodeIndex;
				rec_node->_parentAuxNodeIndex = aux_node_index;
				rec_node->_parentPixNodeIndex = send_node_index;
				BKGCAuxNodeInfo& send_node_aux_info = send_node->_auxNodeInfo[send_node_pix_info._cliqueIndexInNode];
				DEBUG_STATEMENT(PRINTF("Creating PP edge "); BKGCDump::PrintEdge(send_node, rec_node));
				DEBUG_STATEMENT(PRINTF(" through "); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
				send_node_aux_info._childrenCode |= (1 << rec_node_index);
				ResetOutgoingEdges(rec_node_aux_info, send_node_aux_info, send_node_index, aux_node);
				return true;
			}
		}

		//prepare next wavefront
		if(aux_node->_parentIndex == rec_node_aux_info._pixNodeIndexInAuxNode && rec_node_pix_info._numTightConstraints == 0)
			_nextNodeList->push_front(aux_node);
		
		int temp_children = rec_node_aux_info._childrenCode;
		while (temp_children != 0)
		{
			int i = FIRST_NON_ZERO_WORD[temp_children];
			temp_children &= (1 << i) ^ CODE_ALL_ONES_WORD;
			BKGCPixNodeInfo& node_info = aux_node->_nodeInfo[i];
			if((rec_node_pix_info._tightEdgeIntersectionCode & (1 << i)) == 0){
				DEBUG_STATEMENT(BKGCDump::PrintNode(rec_node); PRINTF(" is not putting child "); BKGCDump::PrintNode(node_info._node));
				DEBUG_STATEMENT(PRINTF(" in search list since the edge is saturated and I will not be able to take path from it\n"));
				continue;
			}
			node_info._node->_pathCheckValidPath = false;
			node_info._node->_pathCheckTimeStamp = TIME;
			_nextNodeList->push_back(node_info._node);
		}
	}
	return false;
}

bool BKGCSrcOrphanRootProcessor::FindPathThroughAuxNode(BKGCAuxNode* aux_node, int& old_parent_pix_index)
{
	DEBUG_STATEMENT(PRINTF("finding path through "); BKGCDump::PrintNode(aux_node); PRINTF("\n"));
	int nodes_i_can_seek_path = aux_node->_srcLabeledNodesCode & (aux_node->_childrenCode ^ CODE_ALL_ONES_WORD);
	int i;
	while(nodes_i_can_seek_path != 0)
	{
		i = FIRST_NON_ZERO_WORD[nodes_i_can_seek_path];
		nodes_i_can_seek_path &= (1 << i) ^ CODE_ALL_ONES_WORD;
		BKGCPixNodeInfo& node_info = aux_node->_nodeInfo[i];
		BKGCPixNode* pix_node = node_info._node; 
		if(pix_node->_parentAuxNodeIndex == node_info._cliqueIndexInNode)
			continue;
		if(!pix_node->HasPathToRoot())
			continue;
		old_parent_pix_index = aux_node->_parentIndex;
		aux_node->_parentIndex = i;
		DEBUG_STATEMENT(PRINTF("Creating edge "); BKGCDump::PrintEdge(pix_node, aux_node); PRINTF("\n"));
		return true;
	}
	int temp_labeling = aux_node->_childrenCode;
	while(temp_labeling != 0)
	{
		i = FIRST_NON_ZERO_WORD[temp_labeling];
		temp_labeling &= (1 << i) ^ CODE_ALL_ONES_WORD;
		BKGCPixNodeInfo& node_info = aux_node->_nodeInfo[i];
		node_info._node->_pathCheckValidPath = false;
		node_info._node->_pathCheckTimeStamp = TIME;
		_nextNodeList->push_back(node_info._node);
	}
	return false;
}
