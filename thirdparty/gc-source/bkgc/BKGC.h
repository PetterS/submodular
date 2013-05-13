#pragma once

class BKGCSinkOrphanRootProcessor;
class BKGCSrcOrphanRootProcessor;
#include "bkgc/BKGCDataStructures.h"
#include "util/GadgetGraph.h"

class BKGC : public GadgetGraph
{
public:
	enum EDGE_DIRECTION {REV_DIRECTION=0, FWD_DIRECTION=1};
	typedef std::string Labeling;
	typedef int Cost;
	typedef std::pair<Labeling, Cost> Potential;
	typedef std::pair<BKGCEdge*, EDGE_DIRECTION> D_EDGE;

	BKGC(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node);
	~BKGC(void);

	double FindMaxFlow();

	int GetLabel(int node_id);

	void AddHigherTerm(int* node_indexes, float* labeling_cost);

	void AddUnaryTerm(int index, float e_0_b, float e_1_a);

	double FindCutValue();

	DEBUG_STATEMENT(void DumpState());

private:
	void InitializeConfigVariables();
	void InitializeStucturesBeforeFlow();

	void ProcessActives();
	void ProcessSinkPixActive(BKGCPixNode* pix_node);
	void ProcessSinkAuxActive(BKGCAuxNode* aux_node);
	void ProcessSrcPixActive(BKGCPixNode* pix_node);
	void ProcessSrcAuxActive(BKGCAuxNode* aux_node);

	float AugmentFlow();
	float FindMinResidual(std::list<BKGCEdge>& path);
	void CountConstraints(BKGCEdge& edge, std::map<FlowConstraint, int>& counted_constraints);

	void ProcessTightConstraints();

private:
	double _flow;
	BKGCAuxNode* _auxNodes; int _numAuxNodes;
	BKGCPixNode* _nodes;

	BKGCEdge* _joinEdge;
	BKGCPixNode* _currentPathSinkNode, *_currentPathSrcNode;

	std::list<std::pair<BKGCAuxNode*, int> > _tightConstraints;
	std::list<BKGCNode*> *_activeNodes;

	BKGCSrcOrphanRootProcessor* _srcTreeOrphanRootProcessor;
	BKGCSinkOrphanRootProcessor* _sinkTreeOrphanRootProcessor;
};
