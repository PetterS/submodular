#include "StdAfx.h"
#include "bkgc/BKGCDump.h"

#ifdef CHECK_FOR_ERRORS

void BKGCDump::PrintNode(BKGCAuxNode* node)
{
	PRINTF("(%d,%d,A)",node->_x,node->_y);
}

void BKGCDump::PrintNode(BKGCPixNode* node)
{
	PRINTF("(%d,%d)",node->_x,node->_y);
}

void BKGCDump::PrintEdge(const BKGCEdge& edge)
{
	if(edge._receivingNodeIndex == -1){
		BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node);
		PRINTF("->");
		BKGCDump::PrintNode(edge._auxNode);
	}else if(edge._sendingNodeIndex == -1){
		BKGCDump::PrintNode(edge._auxNode);
		PRINTF("->");
		BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node);
	}else{
		PRINTF("[%d,%d]",edge._auxNode->_x,edge._auxNode->_y);
		BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node);
		PRINTF("->");
		BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node);
	}
}

void BKGCDump::PrintPath(std::list<BKGCEdge>& path)
{
	int size = (int)path.size();
	MY_ASSERT(size > 0);
	PRINTF("[");
	for(std::list<BKGCEdge>::iterator iter = path.begin(); iter!=path.end(); ++iter)
	{
		BKGCEdge& edge = *iter;
		BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._sendingNodeIndex]._node);
		PRINTF("->");
	}
	BKGCEdge& edge = path.back();
	BKGCDump::PrintNode(edge._auxNode->_nodeInfo[edge._receivingNodeIndex]._node);
	PRINTF("]");
}

void BKGCDump::DumpState(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes)
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		BKGCPixNode& node = pix_nodes[i];
		BKGCDump::PrintNode(&node);
		PRINTF(" excess = %f",node._excess);
		PRINTF(" tree = %c",node._tree);
		if(node._parentAuxNodeIndex != -1){
			PRINTF(" parent clique (%d,%d)",node._auxNodeInfo[node._parentAuxNodeIndex]._auxNode->_x, node._auxNodeInfo[node._parentAuxNodeIndex]._auxNode->_y);
			if(node._parentPixNodeIndex != -1){
				BKGCPixNode* parent_node = node._auxNodeInfo[node._parentAuxNodeIndex]._auxNode->_nodeInfo[node._parentPixNodeIndex]._node;
				PRINTF(" parent pix node = "); BKGCDump::PrintNode(parent_node);
			}
		}
		PRINTF("\n");
	}
	for(int i=0;i<num_aux_nodes; ++i)
	{
		BKGCAuxNode& aux_node = aux_nodes[i];
		PRINTLN;
		BKGCDump::PrintNode(&aux_node);
		PRINTF(" parent index %d, tree = %c\n",(int)aux_node._parentIndex, aux_node._tree);
		for(int k=1; k<CODE_ALL_ONES_WORD; ++k)
		{
			float slack = aux_node._constraintSlacks[k];
			if(EQUAL(slack,0)){
				PRINTF("\tconstraint %d is tight \n", k);
			}
		}
	}
}

#endif