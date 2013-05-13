#pragma once

#ifdef CHECK_FOR_ERRORS

class APGCDump
{
public:
	static void PrintNode(APGCPixNode* node);
	static void PrintNode(APGCAuxNode* node);
	static void PrintNode(APGCNode* node){
		if(node->_nodeType == NODE_TYPE_AUX)
			PrintNode((APGCAuxNode*)node);
		else
			PrintNode((APGCPixNode*)node);
	};
	static void PrintEdge(const APGCEdge& edge);
	static void PrintEdge(APGCAuxNode* aux_node, APGCPixNode* send_node, APGCPixNode* rec_node)
	{
		PRINTF("[%d,%d]",aux_node->_x, aux_node->_y);
		PrintNode(send_node); 
		PRINTF("->");
		PrintNode(rec_node); 
	};
	static void PrintEdge(APGCNode* send_node, APGCNode* rec_node)
	{
		PrintNode(send_node); 
		PRINTF("->");
		PrintNode(rec_node); 
	};
	static void UpdateEdge(APGCAuxNode* aux_node, APGCPixNode* send_node, APGCPixNode* rec_node, APGCPixNode* new_send_node, APGCPixNode* new_rec_node)
	{
		PRINTF("Resetting edge "); APGCDump::PrintEdge(aux_node, send_node, rec_node); 
		PRINTF(" to "); APGCDump::PrintEdge(aux_node, new_send_node, new_rec_node); 
	}
	static void UpdateEdge(APGCAuxNode* aux_node, APGCPixNode* send_node, APGCPixNode* rec_node, APGCAuxNode* new_send_node, APGCPixNode* new_rec_node)
	{
		PRINTF("Resetting edge "); APGCDump::PrintEdge(aux_node, send_node, rec_node); 
		PRINTF(" to "); APGCDump::PrintEdge(new_send_node, new_rec_node); 
	}

	static void PrintPath(std::list<APGCEdge>& path);
	static void DumpState(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes);

};

#endif
