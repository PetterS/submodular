#include "StdAfx.h"
#include "prgc/PRGCDump.h"

#ifdef CHECK_FOR_ERRORS

void PRGCDump::PrintClique(PRGCClique* clique)
{
	PRINTF("[%d,%d,%d]", clique->_x, clique->_y, clique->_index);
}

void PRGCDump::PrintNode(PRGCNode* node)
{
	PRINTF("(%d,%d,%d,%.1f)", node->_x, node->_y, node->_distanceLabel, node->_excess) ;
}


void PRGCDump::DumpState(PRGCNode* nodes, int num_pix_nodes, 
						 PRGCClique* cliques, int num_cliques)
{
	PRINTF("=======================\n");
	for(int i=0; i<num_pix_nodes; ++i)
	{
		PRGCDump::PrintNode(&nodes[i]); PRINTLN;
	}
	for(int i=0;i<num_cliques; ++i)
	{
		PRGCClique& clique = cliques[i];
		PRINTLN; PRGCDump::PrintClique(&clique); PRINTLN;
		for(int k=0; k<NUM_CONSTRAINTS; ++k)
			PRINTF("\tconstraint %04d-%s, slack = %f\n", k, BINARY(CONSTRAINT_INFO[k]._labeling), clique._constraintSlacks[k]);
	}
	PRINTLN;
}

#endif