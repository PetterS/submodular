#include "StdAfx.h"
#include "prgc/PRGCDebug.h"
#include "prgc/PRGC.h"

#ifdef CHECK_FOR_ERRORS

void PRGCDebug::ValidateAll(PRGCNode* nodes, int num_pix_nodes, PRGCClique* cliques, int num_cliques)
{
	ValidateCodes(cliques, num_cliques);
	ValidateCodes(nodes, num_pix_nodes);
	ValidateDistanceLabel(nodes, num_pix_nodes);
}

void PRGCDebug::ValidateCodes(PRGCClique* cliques, int num_cliques)
{
	for(int i=0; i<num_cliques; ++i)
	{
		PRGCClique* clique = cliques + i;
		for(int j=0; j<GG_CLIQUE_SIZE; ++j)
		{
			CODE_TYPE tightEdgeIntersectionCode = CODE_ALL_ONES_WORD;
			CODE_TYPE tightEdgeUnionCode = 0;

			int *constraints = CONTAINING_CONSTRAINTS[j];
			for(int i=0; i<NUM_CONTAINING_CONSTRAINTS; ++i)
			{
				int constraint_index = constraints[i];
				if(clique->_constraintSlacks[constraint_index] > EPSILON)
					continue;
				ConstraintInfo& constraint_info = CONSTRAINT_INFO[constraint_index];
				CODE_TYPE constraint_labeling = constraint_info._labeling;
				tightEdgeIntersectionCode &= constraint_labeling;
				tightEdgeUnionCode |= constraint_labeling;
			}

			PRGCEdgeInfo& edge_info = clique->_edgeInfo[j];
			MY_ASSERT(edge_info._tightEdgeIntersectionCode == tightEdgeIntersectionCode);
			MY_ASSERT(edge_info._tightEdgeUnionCode == tightEdgeUnionCode);
		}
	}
}

void PRGCDebug::ValidateCodes(PRGCNode* pix_nodes, int num_pix_nodes)
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		PRGCNode* pix_node = pix_nodes + i;
		for(int j=0; j<pix_node->_numCliques; ++j)
		{
			PRGCEdgeInfo& my_edge_info = *(pix_node->_edgeInfo[j]);
			int my_index = my_edge_info._indexInClique;
			PRGCClique* clique = my_edge_info._clique;

			if(my_edge_info._tightEdgeUnionCode != 0)
				MY_ASSERT((clique->_saturatedEdgeCode & INDEX_CODE_ONE[my_index]) != 0);
			else
				MY_ASSERT((clique->_saturatedEdgeCode & INDEX_CODE_ONE[my_index]) == 0);

			if(pix_node->_distanceLabel != DIST_DEFAULT)
				MY_ASSERT((clique->_labeledNodesCode & INDEX_CODE_ONE[my_index]) != 0);
			else
				MY_ASSERT((clique->_labeledNodesCode & INDEX_CODE_ONE[my_index]) == 0);
		}
	}
}

void PRGCDebug::ValidateDistanceLabel(PRGCNode* pix_nodes, int num_pix_nodes )
{
	for(int i=0; i<num_pix_nodes; ++i)
	{
		PRGCNode* pix_node = pix_nodes + i;

		//DEBUG_COND(ITER_NUM == 359 && pix_node->_x == 8 && pix_node->_y == 8 && pix_node->_type == 'b');

		if(pix_node->_excess < -EPSILON){
			MY_ASSERT(pix_node->_distanceLabel == TERM_DIST);
			continue;
		}
		if(pix_node->_excess > EPSILON){
			if(pix_node->_distanceLabel != DIST_DEFAULT)
				MY_ASSERT(pix_node->_inExcessList == true);
		}
		if(pix_node->_distanceLabel != DIST_DEFAULT)
			continue;

		for(int j=0; j<pix_node->_numCliques; ++j)
		{
			//DEBUG_COND(ITER_NUM >= 205 && pix_node->_x == 24 && pix_node->_y == 4 && pix_node->_type == 'a' && j == 8);

			PRGCEdgeInfo& my_edge_info = *(pix_node->_edgeInfo[j]);
			int my_index = my_edge_info._indexInClique;
			PRGCClique* clique = my_edge_info._clique;

			int not_me = (1 << my_index) ^ CODE_ALL_ONES_WORD;
			CODE_TYPE possible_rec_nodes = clique->_labeledNodesCode & not_me;

			while(possible_rec_nodes != 0)
			{
				int rec_node_index = FIRST_NON_ZERO(possible_rec_nodes);
				possible_rec_nodes &= INDEX_CODE_ZERO[rec_node_index];
				PRGCEdgeInfo& rec_node_edge_info = clique->_edgeInfo[rec_node_index]; 
				if((rec_node_edge_info._tightEdgeIntersectionCode & INDEX_CODE_ONE[my_index]) == 0)
					continue;
				PRGCNode* rec_node = rec_node_edge_info._node; 
				if(rec_node->_distanceLabel >= PRGC_MAX_DIST_LABEL)
					continue;
				if(rec_node->_active)
					continue;
				MY_ASSERT(false);
			}
			
		}
	}
}

#endif
