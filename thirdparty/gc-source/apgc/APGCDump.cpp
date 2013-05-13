#include "StdAfx.h"
#include "apgc/APGCDump.h"

#ifdef CHECK_FOR_ERRORS

void APGCDump::PrintNode(APGCAuxNode* node)
{
	PRINTF("(%d,d=%d,%c,A)", node->_index, node->_distanceLabel, node->_tree);
}

void APGCDump::PrintNode(APGCPixNode* node)
{
	PRINTF("(%d,%d,%c,d=%d,%.1f)", node->_x, node->_y, node->_tree, node->_distanceLabel, node->_excess) ;
}

void APGCDump::PrintEdge(const APGCEdge& edge)
{
	if(edge._receivingNodeIndex == -1){
		APGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node);
		PRINTF("->");
		APGCDump::PrintNode(edge._auxNode);
	}else if(edge._sendingNodeIndex == -1){
		APGCDump::PrintNode(edge._auxNode);
		PRINTF("->");
		APGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node);
	}else{
		PRINTF("[%d,%d]",edge._auxNode->_x, edge._auxNode->_y);
		APGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node);
		PRINTF("->");
		APGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node);
	}
}

void APGCDump::PrintPath(std::list<APGCEdge>& path)
{
	int size = (int)path.size();
	MY_ASSERT(size > 0);
	PRINTF("[");
	for(std::list<APGCEdge>::iterator iter = path.begin(); iter!=path.end(); ++iter)
	{
		APGCEdge& edge = *iter;
		if(edge._sendingNodeIndex != -1)
			APGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node);
		else
			APGCDump::PrintNode(edge._auxNode);
		PRINTF("->");
	}
	APGCEdge& edge = path.back();
	APGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node);
	PRINTF("]");
}

void APGCDump::DumpState(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes)
{
	PRINTF("=======================\n");
	for(int i=0; i<num_pix_nodes; ++i)
	{
		APGCPixNode& node = pix_nodes[i];
		APGCDump::PrintNode(&node);
		if(node._parentAuxNodeIndex != -1){
			PRINTF(" parent clique (%d)",node._auxNodeInfo[node._parentAuxNodeIndex]._auxNode->_index);
			if(node._parentPixNodeIndex != -1){
				APGCPixNode* parent_node = node._auxNodeInfo[node._parentAuxNodeIndex]._auxNode->_nodeInfo[node._parentPixNodeIndex]._node;
				PRINTF(" parent pix node = "); APGCDump::PrintNode(parent_node);
			}
		}
		PRINTLN;
	}
	for(int i=0;i<num_aux_nodes; ++i)
	{
		APGCAuxNode& aux_node = aux_nodes[i];
		PRINTLN;
		APGCDump::PrintNode(&aux_node); PRINTLN;
		for(int k=0; k<NUM_CONSTRAINTS; ++k)
		{
			PRINTF("\tconstraint %03d-%s, slack = %f\n", k, BINARY(CONSTRAINT_INFO[k]._labeling), aux_node._constraintSlacks[k]);
		}
	}
	PRINTLN;
}

#endif