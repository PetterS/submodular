#pragma once

#include "bkgc/BKGCDataStructures.h"

#ifdef CHECK_FOR_ERRORS

class BKGCDump
{
public:
	static void PrintNode(BKGCPixNode* node);
	static void PrintNode(BKGCAuxNode* node);
	static void PrintNode(BKGCNode* node){
		if(node->_nodeType == NODE_TYPE_AUX)
			PrintNode((BKGCAuxNode*)node);
		else
			PrintNode((BKGCPixNode*)node);
	};
	static void PrintEdge(const BKGCEdge& edge);
	static void PrintEdge(BKGCNode* send_node, BKGCNode* rec_node){
		PrintNode(send_node); 
		PRINTF("->");
		PrintNode(rec_node); 
	};
	static void UpdateEdge(BKGCNode* send_node, BKGCNode* rec_node, BKGCNode* new_send_node, BKGCNode* new_rec_node)
	{
		PRINTF("Resetting edge "); BKGCDump::PrintEdge(send_node, rec_node); 
		PRINTF(" to "); BKGCDump::PrintEdge(new_send_node, new_rec_node); PRINTF("\n");
	}
	static void PrintPath(std::list<BKGCEdge>& path);
	static void DumpState(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes);

};

#endif
