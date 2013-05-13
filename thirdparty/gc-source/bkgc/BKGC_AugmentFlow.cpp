#include "StdAfx.h"
#include "bkgc/BKGC.h"
#include "bkgc/BKGCDump.h"
#include "bkgc/BKGCDataStructures.h"

float BKGC::AugmentFlow()
{
	//DEBUG_COND(ITER_NUM == 261);

	std::list<BKGCEdge> path;
	float flow_to_augment = FindMinResidual(path);

	MY_ASSERT(flow_to_augment > EPSILON/GG_CLIQUE_SIZE);
	DEBUG_STATEMENT(PRINTF("augmenting flow = %f\n",flow_to_augment));
	DEBUG_STATEMENT(PRINTF("path = "); BKGCDump::PrintPath(path); PRINTF("\n"));

	MY_ASSERT(_currentPathSinkNode != NULL);

	//DEBUG_COND(ITER_NUM == 21);

	for(std::list<BKGCEdge>::reverse_iterator iter = path.rbegin(); iter != path.rend(); ++iter)
	{
		BKGCEdge& edge = *iter;
		DEBUG_STATEMENT(PRINTF("Increasing flow in edge "); BKGCDump::PrintEdge(edge); PRINTF("\n"));
		edge._auxNode->AugmentFlow(edge._sendingNodeIndex, edge._receivingNodeIndex, flow_to_augment);
	}

	_currentPathSinkNode->_excess += flow_to_augment;

	_currentPathSrcNode->_excess -= flow_to_augment;
	if(_currentPathSrcNode->_excess <= EPSILON && _currentPathSrcNode->_parentAuxNodeIndex != -1){
		DEBUG_STATEMENT(BKGCDump::PrintNode(_currentPathSrcNode); PRINTF(" got saturated, making it active\n"));
		if(!_currentPathSrcNode->_active){
			_currentPathSrcNode->_active = true;
			_activeNodes->push_back(_currentPathSrcNode);
		}
	}

	return flow_to_augment;
}

void BKGC::CountConstraints(BKGCEdge& edge, std::map<FlowConstraint, int>& counted_constraints)
{
	MY_ASSERT(edge._sendingNodeIndex != -1 && edge._receivingNodeIndex != -1 && edge._auxNode != NULL);

	BKGCPixNodeInfo& rec_node_info = edge._auxNode->_nodeInfo[edge._receivingNodeIndex];
	CODE_TYPE sending_node_code = 1 << edge._sendingNodeIndex;
	int* constraints_list = CONTAINING_CONSTRAINTS[edge._receivingNodeIndex];
	for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
	{
		int constraint_index = constraints_list[i];
		if(constraint_index == 0 || constraint_index == CODE_ALL_ONES_WORD)
			continue;
		CODE_TYPE constraint_labeling = CONSTRAINT_LABELING(constraint_index);
		if((constraint_labeling & sending_node_code) == 0)
		{
			MY_ASSERT(edge._auxNode->_constraintSlacks[constraint_index] > EPSILON);
			FlowConstraint constraint_to_add(edge._auxNode->_index, constraint_index);
			std::map<FlowConstraint, int>::iterator iter = counted_constraints.find(constraint_to_add);
			if(iter == counted_constraints.end())
				counted_constraints[constraint_to_add] = 1;
			else
				++(iter->second);
		}
	}
}

float BKGC::FindMinResidual(std::list<BKGCEdge>& path)
{
	//DEBUG_COND(ITER_NUM == 28);

	std::map<FlowConstraint, int> counted_constraints;

	BKGCEdge edge = *_joinEdge;
	edge._auxNode = _joinEdge->_auxNode;
	if(_joinEdge->_sendingNodeIndex == -1){
		MY_ASSERT(_joinEdge->_auxNode->_tree == TREE_SRC);
		MY_ASSERT(_joinEdge->_auxNode->_parentIndex != -1);
		MY_ASSERT(_joinEdge->_receivingNodeIndex != -1);
		edge._sendingNodeIndex = _joinEdge->_auxNode->_parentIndex;
		DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(edge._auxNode));
		DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node); PRINTF("\n"));
	}else if(_joinEdge->_receivingNodeIndex == -1){
		MY_ASSERT(_joinEdge->_auxNode->_tree == TREE_SINK);
		MY_ASSERT(_joinEdge->_auxNode->_parentIndex != -1);
		edge._receivingNodeIndex = _joinEdge->_auxNode->_parentIndex;
		DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node));
		DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode); PRINTF("\n"));
	}else{
		DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node));
		DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node); PRINTF("\n"));
	}
	path.push_back(edge);

	CountConstraints(edge, counted_constraints);

	BKGCPixNode* sending_node = edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node;
	while(true)
	{
		if(sending_node->_excess < -EPSILON){
			_currentPathSinkNode = sending_node;
			break;
		}
		MY_ASSERT(sending_node->_parentAuxNodeIndex != -1);
		BKGCAuxNodeInfo& aux_info = sending_node->_auxNodeInfo[sending_node->_parentAuxNodeIndex];
		edge._sendingNodeIndex = aux_info._pixNodeIndexInAuxNode;
		edge._auxNode = aux_info._auxNode;
		if(sending_node->_parentPixNodeIndex == -1){
			edge._receivingNodeIndex = aux_info._auxNode->_parentIndex;
			DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(sending_node));
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode));
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node); PRINTF("\n"));
		}else{
			edge._receivingNodeIndex = sending_node->_parentPixNodeIndex;
			DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(sending_node));
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node); PRINTF("\n"));
		}
		path.push_back(edge);
		CountConstraints(edge, counted_constraints);
		sending_node = edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node;
	}

	edge = path.front();
	BKGCPixNode* receiving_node = edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node;
	while(true)
	{
		if(receiving_node->_excess > EPSILON){
			_currentPathSrcNode = receiving_node;
			break;
		}
		BKGCAuxNodeInfo& aux_info = receiving_node->_auxNodeInfo[receiving_node->_parentAuxNodeIndex];
		edge._receivingNodeIndex = aux_info._pixNodeIndexInAuxNode;
		edge._auxNode = aux_info._auxNode;
		if(receiving_node->_parentPixNodeIndex == -1){
			edge._sendingNodeIndex = aux_info._auxNode->_parentIndex;
			DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node));
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(edge._auxNode));
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(receiving_node); PRINTF("\n"));
		}else{
			edge._sendingNodeIndex = receiving_node->_parentPixNodeIndex;
			DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node));
			DEBUG_STATEMENT(PRINTF("->"); BKGCDump::PrintNode(receiving_node); PRINTF("\n"));
		}
		path.push_front(edge);
		CountConstraints(edge, counted_constraints);
		receiving_node = aux_info._auxNode->_nodeInfo[edge._sendingNodeIndex]._node;
	}
	
	MY_ASSERT(_currentPathSrcNode->_excess > EPSILON);
	MY_ASSERT(_currentPathSinkNode->_excess < -EPSILON);
	float min_residual = MIN(_currentPathSrcNode->_excess, -_currentPathSinkNode->_excess);
	for(std::map<FlowConstraint, int>::iterator iter=counted_constraints.begin(); iter!=counted_constraints.end(); ++iter)
	{
		FlowConstraint fc = iter->first;
		float slack = _auxNodes[fc._cliqueIndex]._constraintSlacks[fc._constraintIndex];
		MY_ASSERT(slack > EPSILON);
		min_residual = MIN(min_residual, slack/iter->second);
	}

	return min_residual;
}

