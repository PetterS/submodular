#pragma once

#ifdef CHECK_FOR_ERRORS

#include "bkgc/BKGCDataStructures.h"

class BKGCDebug
{
public:
	static void ValidateAll(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes);

private:
	static void ValidatePaths(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateConstraints(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateSaturatedEdges(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateChildren(BKGCPixNode* pix_nodes, int num_pix_nodes, BKGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidatePixNodePath(BKGCPixNode* pix_node);
	static void ValidateAuxNodePath(BKGCAuxNode* aux_node);
	static void ValidateNoP2PPath(BKGCPixNodeInfo& node_info, int my_index, BKGCAuxNode* aux_node);
};

#endif
