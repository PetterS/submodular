#pragma once

#include "util/GadgetGraph.h"

class PRGC : public GadgetGraph
{
public:
	PRGC(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node);
	~PRGC(void);

	double FindMaxFlow();

	inline int GetLabel(int node_id){
		if(_nodes[node_id]._distanceLabel == DIST_DEFAULT)
			return LABEL_SOURCE_0;
		else 
			return LABEL_SINK_1;
	}

	void AddHigherTerm(int* node_indexes, float* labeling_cost);

	void AddUnaryTerm(int index, float e_0_b, float e_1_a);

	void DumpState();

private:
	void InitializeConfigVariables();

	void RelabelAll();
	void ExpandGraph();
	void ExpandNode(PRGCNode* pix_node);
	void PushFlow(PRGCNode* pix_node);
	void PushFlow(PRGCClique* clique);
	float FindResidual(PRGCClique* clique, int send_node_index, int rec_node_index);
	void AdjustConstraintSlacksAfterSendFlow(PRGCClique* clique, int send_node_index, int rec_node_index, float flow);

private:
	int _numCliques;
	double _flow;

	PRGCClique* _cliques;
	PRGCNode* _nodes;

	PRGCPNodeList _excessNodes;
	PRGCPNodeList _activeNodes;
};
