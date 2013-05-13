#pragma once

class PRGCClique;
class PRGCNode;
typedef std::list<PRGCNode*> PRGCPNodeList;

class PRGCEdgeInfo
{
public:
	PRGCEdgeInfo() : _tightEdgeIntersectionCode(CODE_ALL_ONES_WORD), 
		_tightEdgeUnionCode(0){}

	PRGCClique* _clique;
	PRGCNode* _node;
	uchar _indexInClique, _indexInNode;
	CODE_TYPE _tightEdgeIntersectionCode, _tightEdgeUnionCode;
};

class PRGCNode
{
public:
	PRGCNode(void)
		:_numCliques(0), _excess(0), _distanceLabel(DIST_DEFAULT), _inExcessList(false)
	{
		_edgeInfo = new PRGCEdgeInfo*[NUM_CLIQUES_PER_NODE];
		_active = false;
	}
	~PRGCNode() { 
		delete[] _edgeInfo; 
	}
	void SetAsSink(PRGCPNodeList& active_list);
	void DeclareOrphan();
	void FindNewDistanceLabel();
	void ReduceExcess(float flow);
	void IncreaseExcess(float flow, PRGCPNodeList& excess_list, PRGCPNodeList& active_list);
public:
	DEBUG_STATEMENT(int _x);
	DEBUG_STATEMENT(int _y);

	PRGCEdgeInfo** _edgeInfo;
	uchar _numCliques;
	bool _inExcessList;
	float _excess;
	DIST_TYPE _distanceLabel;
	bool _active;
};

class PRGCClique
{
public:
	PRGCClique(void)
		:_labeledNodesCode(0), _saturatedEdgeCode(0)
	{
		_edgeInfo = new PRGCEdgeInfo[GG_CLIQUE_SIZE];
		_constraintSlacks = new float[NUM_CONSTRAINTS];
		memset(_constraintSlacks, 0, NUM_CONSTRAINTS*sizeof(float));
	};
	~PRGCClique(){
		delete[] _edgeInfo;
		delete[] _constraintSlacks;
	}
	inline void InsertNode(PRGCNode* node, int node_index)
	{
		PRGCEdgeInfo* edge_info = &_edgeInfo[node_index];
		node->_edgeInfo[node->_numCliques] = edge_info;

		edge_info->_node = node;
		edge_info->_clique = this;
		edge_info->_indexInClique = node_index;
		edge_info->_indexInNode = node->_numCliques++;
	}
	void RecalculateTightEdgeCodes(int index);
	void RecalculateSaturatedEdgeCode();
	void Reparametrize(PRGCPNodeList& excess_list, PRGCPNodeList& active_list, double& flow);
	void DeclareConstraintNotTight(int constraint_index, CODE_TYPE constraint_labeling, PRGCPNodeList& active_list);
	void DeclareConstraintTight(int constraint_index, CODE_TYPE constraint_labeling);

private:
	float FindMinContainingSlack(int node_index);
	float FindMinNonContainingSlack(int node_index);

public:
	DEBUG_STATEMENT(int _x);
	DEBUG_STATEMENT(int _y);
	DEBUG_STATEMENT(int _index);

	PRGCEdgeInfo* _edgeInfo;
	float* _constraintSlacks;

	CODE_TYPE _labeledNodesCode;
	CODE_TYPE _saturatedEdgeCode;
};
