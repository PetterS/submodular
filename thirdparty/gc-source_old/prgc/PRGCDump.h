#pragma once

#ifdef CHECK_FOR_ERRORS

class PRGCDump
{
public:
	static void PrintNode(PRGCNode* node);
	static void PrintClique(PRGCClique* node);
	static void DumpState(PRGCNode* nodes, int num_pix_nodes, PRGCClique* cliques, int num_cliques);

};

#endif
