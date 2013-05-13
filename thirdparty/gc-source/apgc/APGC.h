#pragma once

#include "util/GadgetGraph.h"

class APGC : public GadgetGraph
{
public:
	APGC(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node);
	~APGC(void);

	double FindMaxFlow();

	inline int GetLabel(int node_id){
		if(_nodes[node_id]._tree == TREE_SRC)
			return LABEL_SOURCE_0;
		else 
			return LABEL_SINK_1;
	}

	void AddHigherTerm(int* node_indexes, float* labeling_cost);

	void AddUnaryTerm(int index, float e_0_b, float e_1_a);

	void DumpState();

private:
	void InitializeConfigVariables();
	void InitializeStucturesBeforeFlow();

	void ExpandGraph();
	void ProcessSrcActive(APGCPixNode* pix_node);
	void ProcessSinkActive(APGCPixNode* pix_node);
	void ProcessSrcActive(APGCAuxNode* aux_node);
	void ProcessSinkActive(APGCAuxNode* aux_node);

	float AugmentFlow();
	float FindMinResidual(std::list<APGCEdge>& src_side_path, std::list<APGCEdge>& sink_side_path);
	void CountConstraints(APGCEdge& edge, std::map<FlowConstraint, int>& counted_constraints);

	void ProcessTightConstraints();
	void ProcessOrphans();
	void DeclareActive(APGCNode* node);
	void AssignDistLabel(APGCNode* node, int dist);
	//void DeclareActive(APGCNode* node, APGCNodeList& next_actives);

	void ProcessSrcOrphan(APGCPixNode* node);
	void ProcessSinkOrphan(APGCPixNode* node);
	void RemoveLinkToParentAndPushInOraphanList(APGCPixNode* pix_node);
	void DeclareChildrenOrphan(APGCPixNode* node);
	void DeclareFree(APGCPixNode* node);

	void ProcessSrcOrphan(APGCAuxNode* aux_node);
	void ProcessSinkOrphan(APGCAuxNode* aux_node);
	void RemoveLinkToParentAndPushInOraphanList(APGCAuxNode* aux_node);
	void DeclareChildrenOrphan(APGCAuxNode* aux_node);
	void DeclareFree(APGCAuxNode* aux_node);
	
private:
	double _flow;
	APGCAuxNode* _auxNodes; int _numAuxNodes;
	APGCPixNode* _nodes;

	APGCEdge* _joinEdge;
	APGCPixNode* _currentPathSinkNode, *_currentPathSrcNode;

	std::vector<APGCNodeList> _activeNodes, _orphanNodes;
	int _maxAssignedDistLabel;

	int _lastPathLength;

	int _foundPathLength;
};
