#pragma once

#ifdef CHECK_FOR_ERRORS

class APGCDebug
{
public:
	static void ValidateAll(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes);

private:
	static void ValidatePaths(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateDistanceLabel(APGCPixNode* pix_nodes, int num_pix_nodes);
	static void ValidateDistanceLabel(APGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateConstraints(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateSaturatedEdges(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidateChildren(APGCPixNode* pix_nodes, int num_pix_nodes, APGCAuxNode* aux_nodes, int num_aux_nodes);
	static void ValidatePixNodePath(APGCPixNode* pix_node);
	static void ValidateAuxNodePath(APGCAuxNode* aux_node);
	static void ValidateNoP2PPath(APGCPixNodeInfo& node_info, int my_index, APGCAuxNode* aux_node);
};

#endif
