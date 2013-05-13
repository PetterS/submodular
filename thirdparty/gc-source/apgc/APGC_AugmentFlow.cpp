#include "StdAfx.h"
#include "apgc/APGC.h"
#include "apgc/APGCDump.h"

float APGC::AugmentFlow()
{
	//DEBUG_COND(ITER_NUM == 261);

	std::list<APGCEdge> src_side_path, sink_side_path;
	float flow_to_augment = FindMinResidual(src_side_path, sink_side_path);

	MY_ASSERT(flow_to_augment > EPSILON/GG_CLIQUE_SIZE);
	DEBUG_STATEMENT(PRINTF("\nAugmenting flow = %f\n",flow_to_augment));

	DEBUG_STATEMENT(std::list<APGCEdge> complete_path = src_side_path);
	DEBUG_STATEMENT(complete_path.insert(complete_path.end(), sink_side_path.begin(), sink_side_path.end()));
	DEBUG_STATEMENT(PRINTF("path = "); APGCDump::PrintPath(complete_path); PRINTLN);

	int current_path_length = (int)(src_side_path.size()+sink_side_path.size());
	MY_ASSERT(current_path_length >= _lastPathLength);
	_lastPathLength = current_path_length;

	MY_ASSERT(_currentPathSinkNode != NULL);

	//DEBUG_COND(ITER_NUM == 21);

	_currentPathSinkNode->_excess += flow_to_augment;
	if(EQUAL(_currentPathSinkNode->_excess,0)){
		DEBUG_STATEMENT(PRINTF("\n"); APGCDump::PrintNode(_currentPathSinkNode); PRINTF(" became saturated\n"));
		RemoveLinkToParentAndPushInOraphanList(_currentPathSinkNode);
	}

	_currentPathSrcNode->_excess -= flow_to_augment;
	if(EQUAL(_currentPathSrcNode->_excess,0)){
		DEBUG_STATEMENT(PRINTF("\n"); APGCDump::PrintNode(_currentPathSrcNode); PRINTF(" became saturated\n"));
		RemoveLinkToParentAndPushInOraphanList(_currentPathSrcNode);
	}

	for(std::list<APGCEdge>::reverse_iterator iter = sink_side_path.rbegin(); 
		iter != sink_side_path.rend(); ++iter)
	{
		APGCEdge& edge = *iter;
		DEBUG_STATEMENT(PRINTF("Increasing flow in edge "); APGCDump::PrintEdge(edge); PRINTLN);
		edge._auxNode->AugmentFlow(edge._sendingNodeIndex, edge._receivingNodeIndex, flow_to_augment);
	}

	for(std::list<APGCEdge>::iterator iter = src_side_path.begin(); 
		iter != src_side_path.end(); ++iter)
	{
		APGCEdge& edge = *iter;
		DEBUG_STATEMENT(PRINTF("Increasing flow in edge "); APGCDump::PrintEdge(edge); PRINTLN);
		edge._auxNode->AugmentFlow(edge._sendingNodeIndex, edge._receivingNodeIndex, flow_to_augment);
	}

	return flow_to_augment;
}

void APGC::CountConstraints(APGCEdge& edge, std::map<FlowConstraint, int>& counted_constraints)
{
	MY_ASSERT(edge._sendingNodeIndex != -1 && edge._receivingNodeIndex != -1 && edge._auxNode != NULL);

	APGCPixNodeInfo& rec_node_info = edge._auxNode->_nodeInfo[edge._receivingNodeIndex];
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

float APGC::FindMinResidual(std::list<APGCEdge>& src_side_path,
								std::list<APGCEdge>& sink_side_path)
{
	std::map<FlowConstraint, int> counted_constraints;

	//DEBUG_COND(ITER_NUM == 28);

	APGCEdge edge = *_joinEdge;
	if(_joinEdge->_sendingNodeIndex == -1){
		MY_ASSERT(_joinEdge->_auxNode->_tree == TREE_SRC);
		MY_ASSERT(_joinEdge->_auxNode->_parentPixNodeIndex != -1);
		MY_ASSERT(_joinEdge->_receivingNodeIndex != -1);
		edge._sendingNodeIndex = _joinEdge->_auxNode->_parentPixNodeIndex;
	}else if(_joinEdge->_receivingNodeIndex == -1){
		MY_ASSERT(_joinEdge->_auxNode->_tree == TREE_SINK);
		MY_ASSERT(_joinEdge->_auxNode->_parentPixNodeIndex != -1);
		edge._receivingNodeIndex = _joinEdge->_auxNode->_parentPixNodeIndex;
	}

	DEBUG_STATEMENT(PRINTF("\nPopulating path - pushed edge "); APGCDump::PrintEdge(edge); PRINTLN);

	src_side_path.push_back(edge);

	CountConstraints(edge, counted_constraints);

	APGCPixNode* sending_node = edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node;
	while(true)
	{
		if(sending_node->_excess < -EPSILON){
			_currentPathSinkNode = sending_node;
			break;
		}
		MY_ASSERT(sending_node->_parentAuxNodeIndex != -1);
		APGCAuxNodeInfo& aux_info = sending_node->_auxNodeInfo[sending_node->_parentAuxNodeIndex];
		edge._sendingNodeIndex = aux_info._pixNodeIndexInAuxNode;
		edge._auxNode = aux_info._auxNode;
		if(sending_node->_parentPixNodeIndex == -1)
			edge._receivingNodeIndex = aux_info._auxNode->_parentPixNodeIndex;
		else
			edge._receivingNodeIndex = sending_node->_parentPixNodeIndex;
		DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); APGCDump::PrintEdge(edge); PRINTLN);
		sink_side_path.push_back(edge);
		CountConstraints(edge, counted_constraints);
		sending_node = edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node;
	}

	edge = src_side_path.front();
	APGCPixNode* receiving_node = edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node;
	while(true)
	{
		if(receiving_node->_excess > EPSILON){
			_currentPathSrcNode = receiving_node;
			break;
		}

		APGCAuxNodeInfo& aux_info = receiving_node->_auxNodeInfo[receiving_node->_parentAuxNodeIndex];
		edge._receivingNodeIndex = aux_info._pixNodeIndexInAuxNode;
		edge._auxNode = aux_info._auxNode;
		if(receiving_node->_parentPixNodeIndex == -1)
			edge._sendingNodeIndex = aux_info._auxNode->_parentPixNodeIndex;
		else
			edge._sendingNodeIndex = receiving_node->_parentPixNodeIndex;
		DEBUG_STATEMENT(PRINTF("populating path - pushed edge "); APGCDump::PrintEdge(edge); PRINTLN);
		src_side_path.push_front(edge);
		CountConstraints(edge, counted_constraints);
		receiving_node = edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node;
	}
	
	MY_ASSERT(_currentPathSinkNode->_excess < -EPSILON && _currentPathSrcNode->_excess > EPSILON);
	float min_residual = MIN(-_currentPathSinkNode->_excess, _currentPathSrcNode->_excess);
	for(std::map<FlowConstraint, int>::iterator iter=counted_constraints.begin(); iter!=counted_constraints.end(); ++iter)
	{
		FlowConstraint fc = iter->first;
		float slack = _auxNodes[fc._cliqueIndex]._constraintSlacks[fc._constraintIndex];
		MY_ASSERT(slack > EPSILON);
		min_residual = MIN(min_residual, slack/iter->second);
	}

	return min_residual;
}

