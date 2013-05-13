#pragma once

class GadgetGraph
{
public:
	//Constructs a gadget graph with num_ndoes nodes, num_cliques cliques. 
	//Nodes are numbered 0-(num_nodes-1)
	//Size of the cliques is specified with clique_size
	//A node may participate in num_cliques_per_node or lesser number of cliques.
	GadgetGraph(int num_nodes, int num_cliques, int clique_size, int num_cliques_per_node);

	virtual ~GadgetGraph();

	// Adds unary term with cost values e_0_b for labeling 0 and e_1_a for labeling 1.
	virtual void AddUnaryTerm(int index, float e_0_b, float e_1_a) = 0;

	//Adds a new clique 
	//node_indexes specifies nodes (0 based index) participating in the clique.
	//It is expected that the length of node_indexes will be exctly equal to clique_size
	//labeling cost specifies clique potential such that labeling_cost[7] represents costs for 
	//labeling nodes of the clique '0111'.
	//Leftmost digit represents label of first node. 
	//e.g. '0111' represents first node (in node_indexes) is labeled 0 and rest 1.
	virtual void AddHigherTerm(int* node_indexes, float* labeling_cost) = 0;

	//finds and returns the value of maximum flow. If the clique potential specified is submodular for all cliques
	//then the returned value should be equal to value of primal.
	virtual double FindMaxFlow() = 0;

	//returns the label of node node_id (0/1).
	virtual int GetLabel(int node_id) = 0;

	//debug function to dump the state of the flowgraph
	virtual void DumpState() {};
};
