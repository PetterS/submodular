#include "StdAfx.h"
#include "apgc/APGC.h"
#include "apgc/APGCDump.h"

void APGC::ExpandGraph()
{
	PRINTF("Max node dist=%d\n", _maxAssignedDistLabel);

	for(int level = 0; level <= _maxAssignedDistLabel; ++level)
	{
		MY_ASSERT((int)_activeNodes.size() > level);
		if(_activeNodes[level].empty()){
			if(level == _maxAssignedDistLabel)
				break;
			continue;
		}

		while((int)_activeNodes.size() <= (level+2)){
			_activeNodes.push_back(APGCNodeList());
			_orphanNodes.push_back(APGCNodeList());
		}

		PRINTF("\nProcessing active nodes at level %d\n", level);

		std::multimap<int, APGCEdge> found_paths;

		APGCNodeList& current_actives = _activeNodes[level];
		std::vector<APGCNode*> backup_nodes;

		while(!current_actives.empty())
		{
			_foundPathLength = INT_MAX;

			APGCNode* gnode = current_actives.front();
			if(gnode->_distanceLabel == DIST_DEFAULT){
				gnode->_active = false;
			}else if(gnode->_distanceLabel == level)
			{
				if(gnode->_nodeType == NODE_TYPE_AUX){			
					if(gnode->_tree == TREE_SRC)
						ProcessSrcActive((APGCAuxNode*)gnode);
					else
						ProcessSinkActive((APGCAuxNode*)gnode);
				}else{
					if(gnode->_tree == TREE_SRC)
						ProcessSrcActive((APGCPixNode*)gnode);
					else
						ProcessSinkActive((APGCPixNode*)gnode);
				}
				if(_foundPathLength != INT_MAX)
				{
					found_paths.insert(std::pair<int, APGCEdge>(_foundPathLength, *_joinEdge));
					_joinEdge->_auxNode = NULL;
					backup_nodes.push_back(gnode);
				}
				gnode->_active = false;
			}else{
				DEBUG_STATEMENT(if(gnode->_distanceLabel < level){ APGCDump::PrintNode(gnode); MY_ASSERT(false);});
				DEBUG_STATEMENT(PRINTF("Node "); APGCDump::PrintNode(gnode); PRINTF(" at active level %d, moving to %d\n", level, gnode->_distanceLabel));
				_activeNodes[gnode->_distanceLabel].push_back(gnode);
			}
			current_actives.pop_front();
		}

		if(!found_paths.empty())
		{
			MY_ASSERT(!backup_nodes.empty());
			int num_backup_nodes = (int)backup_nodes.size();
			for(int i=0; i<num_backup_nodes; ++i)
			{
				APGCNode* node = backup_nodes[i];
				MY_ASSERT(!node->_active);
				node->_active = true;
				current_actives.push_back(node);
			}
			std::multimap<int, APGCEdge>::iterator iter = found_paths.begin();
			*_joinEdge = iter->second;
			DEBUG_STATEMENT(PRINTF("found joining edge "); APGCDump::PrintEdge(*_joinEdge); PRINTLN);
			return;
		}else
		{
			MY_ASSERT(backup_nodes.empty());
		}
	}
}

void APGC::ProcessSinkActive(APGCPixNode* node)
{
	MY_ASSERT(node->_active);
	DEBUG_STATEMENT(PRINTF("\nProcessing sink active node "); APGCDump::PrintNode(node); PRINTLN);
	MY_ASSERT(node->_tree == TREE_SINK);

	//DEBUG_COND(ITER_NUM == 100 && node->_x == 4 && node->_y == 8);

	for(int i=0; i<node->_numAuxNodes; ++i)
	{
		if(i == node->_parentAuxNodeIndex)
			continue;
		APGCAuxNodeInfo& receiving_node_aux_info = node->_auxNodeInfo[i];
		APGCAuxNode* aux_node = receiving_node_aux_info._auxNode;
		int receiving_node_index = receiving_node_aux_info._pixNodeIndexInAuxNode;
		APGCPixNodeInfo& receiving_node_info = aux_node->_nodeInfo[receiving_node_index];
		if(receiving_node_info._numTightConstraints == 0)
		{
			if(aux_node->_parentPixNodeIndex == -1)
			{
				DEBUG_STATEMENT(PRINTF("Creating edge "); APGCDump::PrintEdge(aux_node, node); PRINTLN); 
				aux_node->_parentPixNodeIndex = receiving_node_index;
				aux_node->_tree = TREE_SINK;
				AssignDistLabel(aux_node, node->_distanceLabel+1);
			}else
			{
				if(aux_node->_tree == TREE_SRC)	
				{
					int expected_path_length = aux_node->_distanceLabel + node->_distanceLabel + 1;
					if(expected_path_length < _foundPathLength)
					{
						_foundPathLength = expected_path_length;
						_joinEdge->_auxNode = aux_node;
						_joinEdge->_sendingNodeIndex = -1;
						_joinEdge->_receivingNodeIndex = receiving_node_index;
					}
				}
			}
		}else
		{
			CODE_TYPE possible_join_nodes = aux_node->_srcLabeledNodesCode & receiving_node_info._tightEdgeIntersectionCode;
			while(possible_join_nodes != 0)
			{
				int index = FIRST_NON_ZERO(possible_join_nodes);
				possible_join_nodes &= INDEX_CODE_ZERO[index];
				APGCPixNode* send_node = aux_node->_nodeInfo[index]._node;
				int expected_path_length = send_node->_distanceLabel + node->_distanceLabel + 2;
				if(expected_path_length < _foundPathLength)
				{
					_foundPathLength = expected_path_length;
					_joinEdge->_auxNode = aux_node;
					_joinEdge->_sendingNodeIndex = index;
					_joinEdge->_receivingNodeIndex = receiving_node_index;
				}
			}

			CODE_TYPE unlabeled_nodes = (aux_node->_sinkLabeledNodesCode | aux_node->_srcLabeledNodesCode) ^ CODE_ALL_ONES_WORD;
			CODE_TYPE nodes_i_can_give_label = unlabeled_nodes & receiving_node_info._tightEdgeIntersectionCode;
			while(nodes_i_can_give_label != 0)
			{
				uchar j = FIRST_NON_ZERO(nodes_i_can_give_label);
				nodes_i_can_give_label &= (1 << j) ^ CODE_ALL_ONES_WORD;
				APGCPixNodeInfo& sending_node_info = aux_node->_nodeInfo[j];
				APGCPixNode* sending_node = sending_node_info._node;
				MY_ASSERT(sending_node->_parentAuxNodeIndex == -1);
				MY_ASSERT(sending_node->_tree == TREE_NONE);
				sending_node->_parentAuxNodeIndex = sending_node_info._cliqueIndexInNode;
				sending_node->_parentPixNodeIndex = receiving_node_index;
				receiving_node_aux_info._childrenCode |= (1 << j);
				DEBUG_STATEMENT(PRINTF("Creating edge "); APGCDump::PrintEdge(aux_node, sending_node, node); PRINTLN);
				
				sending_node->_tree = TREE_SINK;
				AssignDistLabel(sending_node, node->_distanceLabel+2);
				for(int i=0; i<sending_node->_numAuxNodes; ++i){
					APGCAuxNodeInfo& aux_info = sending_node->_auxNodeInfo[i];
					aux_info._auxNode->_sinkLabeledNodesCode |= INDEX_CODE_ONE[aux_info._pixNodeIndexInAuxNode];
				}
			}
		}
	}
}

void APGC::ProcessSrcActive(APGCPixNode* node)
{
	MY_ASSERT(node->_active);
	DEBUG_STATEMENT(PRINTF("\nProcessing src active node "); APGCDump::PrintNode(node); PRINTLN);
	MY_ASSERT(node->_tree == TREE_SRC);
	for(int i=0; i<node->_numAuxNodes; ++i)
	{
		//DEBUG_COND(ITER_NUM == 9 && node->_x == 1 && node->_y == 1);

		if(i == node->_parentAuxNodeIndex)
			continue;
		APGCAuxNodeInfo& sending_node_aux_info = node->_auxNodeInfo[i];
		APGCAuxNode* aux_node = sending_node_aux_info._auxNode;
		int sending_node_index = sending_node_aux_info._pixNodeIndexInAuxNode;
		APGCPixNodeInfo& sending_node_info = aux_node->_nodeInfo[sending_node_index];

		if(aux_node->_parentPixNodeIndex == -1)
		{
			DEBUG_STATEMENT(PRINTF("Creating edge "); APGCDump::PrintEdge(node,aux_node); PRINTLN);
			aux_node->_parentPixNodeIndex = sending_node_index;
			aux_node->_tree = TREE_SRC;
			AssignDistLabel(aux_node, node->_distanceLabel+1);
		}else
		{
			if(aux_node->_tree == TREE_SINK)	
			{
				int expected_path_length = aux_node->_distanceLabel + node->_distanceLabel + 1;
				if(expected_path_length < _foundPathLength)
				{
					_foundPathLength = expected_path_length;
					_joinEdge->_auxNode = aux_node;
					_joinEdge->_sendingNodeIndex = sending_node_index;
					_joinEdge->_receivingNodeIndex = -1;
				}
			}
		}

		if(sending_node_info._numTightConstraints > 0)
		{
			CODE_TYPE possible_join_nodes = aux_node->_sinkLabeledNodesCode & sending_node_info._tightEdgeUnionCode;
			while(possible_join_nodes != 0)
			{
				int j = FIRST_NON_ZERO(possible_join_nodes);
				possible_join_nodes &= (1 << j) ^ CODE_ALL_ONES_WORD;
				APGCPixNodeInfo& rec_node_info = aux_node->_nodeInfo[j];
				if((rec_node_info._tightEdgeIntersectionCode & INDEX_CODE_ONE[sending_node_index]) != 0)
				{
					int expected_path_length = rec_node_info._node->_distanceLabel + node->_distanceLabel + 2;
					if(expected_path_length < _foundPathLength)
					{
						_foundPathLength = expected_path_length;
						_joinEdge->_auxNode = aux_node;
						_joinEdge->_sendingNodeIndex = sending_node_index;
						_joinEdge->_receivingNodeIndex = j;
					}
				}
			}

			CODE_TYPE unlabeled_not_me_nodes = (aux_node->_sinkLabeledNodesCode | aux_node->_srcLabeledNodesCode 
									| (1 << sending_node_index)) ^ CODE_ALL_ONES_WORD;
			CODE_TYPE nodes_i_can_give_label = unlabeled_not_me_nodes & sending_node_info._tightEdgeUnionCode;
			while(nodes_i_can_give_label != 0)
			{
				int j = FIRST_NON_ZERO(nodes_i_can_give_label);
				nodes_i_can_give_label &= (1 << j) ^ CODE_ALL_ONES_WORD;
				APGCPixNodeInfo& rec_node_info = aux_node->_nodeInfo[j];
				if((rec_node_info._tightEdgeIntersectionCode & (1 << sending_node_index)) == 0)
					continue;
				APGCPixNode* rec_node = rec_node_info._node;
				MY_ASSERT(rec_node->_tree == TREE_NONE);
				rec_node->_parentAuxNodeIndex = rec_node_info._cliqueIndexInNode;
				rec_node->_parentPixNodeIndex = sending_node_index;
				sending_node_aux_info._childrenCode |= (1 << j);
				DEBUG_STATEMENT(PRINTF("Creating edge "); APGCDump::PrintEdge(aux_node, sending_node_info._node,rec_node); PRINTLN);
				rec_node->_tree = TREE_SRC;
				AssignDistLabel(rec_node, node->_distanceLabel+2);	
				for(int i=0; i<rec_node->_numAuxNodes; ++i){
					APGCAuxNodeInfo& aux_info = rec_node->_auxNodeInfo[i];
					aux_info._auxNode->_srcLabeledNodesCode |= INDEX_CODE_ONE[aux_info._pixNodeIndexInAuxNode];
				}
			}
		}
	}
}

void APGC::ProcessSrcActive(APGCAuxNode* aux_node)
{
	MY_ASSERT(aux_node->_active);
	DEBUG_STATEMENT(PRINTF("\nProcessing src active node "); APGCDump::PrintNode(aux_node); PRINTLN);
	//DEBUG_COND(ITER_NUM == 8);
	MY_ASSERT(aux_node->_parentPixNodeIndex != -1 && aux_node->_distanceLabel != DIST_DEFAULT);
	CODE_TYPE sink_labeled_unsaturated_nodes = aux_node->_sinkLabeledNodesCode & (aux_node->_saturatedEdgeCode ^ CODE_ALL_ONES_WORD);
	while(sink_labeled_unsaturated_nodes != 0)
	{
		int rec_node_index = FIRST_NON_ZERO(sink_labeled_unsaturated_nodes);
		sink_labeled_unsaturated_nodes &= INDEX_CODE_ZERO[rec_node_index];
		APGCPixNode* rec_node = aux_node->_nodeInfo[rec_node_index]._node;
		MY_ASSERT(rec_node->_tree == TREE_SINK);
		int expected_path_length = aux_node->_distanceLabel + rec_node->_distanceLabel + 1;
		if(expected_path_length < _foundPathLength)
		{
			_foundPathLength = expected_path_length;
			_joinEdge->_auxNode = aux_node;
			_joinEdge->_sendingNodeIndex = -1;
			_joinEdge->_receivingNodeIndex = rec_node_index;
		}
	}

	CODE_TYPE unlabeled_unsaturated_nodes = (aux_node->_sinkLabeledNodesCode | aux_node->_srcLabeledNodesCode | aux_node->_saturatedEdgeCode) ^ CODE_ALL_ONES_WORD;
	while(unlabeled_unsaturated_nodes != 0)
	{
		int j = FIRST_NON_ZERO(unlabeled_unsaturated_nodes);
		unlabeled_unsaturated_nodes &= (1 << j) ^ CODE_ALL_ONES_WORD;
		APGCPixNodeInfo& rec_node_info = aux_node->_nodeInfo[j];
		APGCPixNode* rec_node = rec_node_info._node; MY_ASSERT(rec_node->_tree == TREE_NONE);
		rec_node->_parentAuxNodeIndex = rec_node_info._cliqueIndexInNode;
		aux_node->_childrenCode |= (1 << j);
		DEBUG_STATEMENT(PRINTF("Creating edge "); APGCDump::PrintEdge(aux_node,rec_node); PRINTLN);
		for(int i=0; i<rec_node->_numAuxNodes; ++i){
			APGCAuxNodeInfo& aux_info = rec_node->_auxNodeInfo[i];
			aux_info._auxNode->_srcLabeledNodesCode |= (1 << aux_info._pixNodeIndexInAuxNode);
		}
		rec_node->_tree = TREE_SRC;
		AssignDistLabel(rec_node, aux_node->_distanceLabel + 1);	
		for(int i=0; i<rec_node->_numAuxNodes; ++i){
			APGCAuxNodeInfo& aux_info = rec_node->_auxNodeInfo[i];
			aux_info._auxNode->_srcLabeledNodesCode |= INDEX_CODE_ONE[aux_info._pixNodeIndexInAuxNode];
		}
	}
}

void APGC::ProcessSinkActive(APGCAuxNode* aux_node)
{
	MY_ASSERT(aux_node->_active);
	DEBUG_STATEMENT(PRINTF("\nProcessing sink active node "); APGCDump::PrintNode(aux_node); PRINTLN);
	//DEBUG_COND(ITER_NUM == 9);
	MY_ASSERT(aux_node->_parentPixNodeIndex != -1 && aux_node->_distanceLabel != DIST_DEFAULT); 
	CODE_TYPE src_labeled_nodes = aux_node->_srcLabeledNodesCode;
	while(src_labeled_nodes != 0)
	{
		int send_node_index = FIRST_NON_ZERO(src_labeled_nodes);
		src_labeled_nodes &= INDEX_CODE_ZERO[send_node_index];
		APGCPixNode* src_node = aux_node->_nodeInfo[send_node_index]._node;
		MY_ASSERT(src_node->_tree == TREE_SRC);
		int expected_path_length = aux_node->_distanceLabel + src_node->_distanceLabel + 1;
		if(expected_path_length < _foundPathLength)
		{
			_foundPathLength = expected_path_length;
			_joinEdge->_auxNode = aux_node;
			_joinEdge->_receivingNodeIndex = -1;
			_joinEdge->_sendingNodeIndex = send_node_index;
		}
	}

	int unlabeled_nodes = (aux_node->_sinkLabeledNodesCode | aux_node->_srcLabeledNodesCode) ^ CODE_ALL_ONES_WORD;
	while(unlabeled_nodes != 0)
	{
		uchar j = FIRST_NON_ZERO_WORD[unlabeled_nodes];
		unlabeled_nodes &= (1 << j) ^ CODE_ALL_ONES_WORD;
		APGCPixNodeInfo& sending_node_info = aux_node->_nodeInfo[j];
		APGCPixNode* sending_node = sending_node_info._node;
		MY_ASSERT(sending_node->_parentAuxNodeIndex == -1);
		DEBUG_STATEMENT(PRINTF("Creating edge "); APGCDump::PrintEdge(sending_node,aux_node); PRINTLN);
		sending_node->_parentAuxNodeIndex = sending_node_info._cliqueIndexInNode;
		aux_node->_childrenCode |= (1 << j);
		sending_node->_tree = TREE_SINK;
		AssignDistLabel(sending_node, aux_node->_distanceLabel+1);
		for(int i=0; i<sending_node->_numAuxNodes; ++i){
			APGCAuxNodeInfo& aux_info = sending_node->_auxNodeInfo[i];
			aux_info._auxNode->_sinkLabeledNodesCode |= INDEX_CODE_ONE[aux_info._pixNodeIndexInAuxNode];
		}
	}
}

void APGC::DeclareActive(APGCNode* node)
{
	MY_ASSERT(node->_tree != TREE_NONE);
	if(node->_active)
		return;
	DEBUG_STATEMENT(PRINTF("Pushing node"); APGCDump::PrintNode(node); PRINTF(" in active list\n"));
	node->_active = true;
	_activeNodes[node->_distanceLabel].push_back(node);
}

void APGC::AssignDistLabel(APGCNode* node, int dist)
{
	MY_ASSERT(node->_distanceLabel < dist);
	node->_distanceLabel = dist;
	if(_maxAssignedDistLabel < node->_distanceLabel)
		_maxAssignedDistLabel = node->_distanceLabel;
	DeclareActive(node);
}