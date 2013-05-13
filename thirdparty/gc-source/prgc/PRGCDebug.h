#pragma once

#ifdef CHECK_FOR_ERRORS

class PRGCDebug
{
public:
	static void ValidateAll(PRGCNode* nodes, int num_pix_nodes, PRGCClique* cliques, int num_cliques);

private:
	static void ValidateDistanceLabel(PRGCNode* pix_nodes, int num_pix_nodes);
	static void ValidateCodes(PRGCNode* pix_nodes, int num_pix_nodes);
	static void ValidateCodes(PRGCClique* cliques, int num_cliques);
};

#endif
