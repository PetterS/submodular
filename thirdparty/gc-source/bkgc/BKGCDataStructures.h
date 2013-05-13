#pragma once

class BKGCAuxNode;

class BKGCNode
{
public:
	char _nodeType;
	char _tree;
	bool _active;	
};

class BKGCEdge
{
public:
	BKGCAuxNode* _auxNode;
	char _sendingNodeIndex;
	char _receivingNodeIndex;
};

class BKGCAuxNodeInfo
{
public:
	BKGCAuxNode* _auxNode;
	uchar _pixNodeIndexInAuxNode;
	CODE_TYPE _childrenCode;
};

class BKGCPixNode : public BKGCNode
{
public:
	BKGCPixNode(void)
		:_numAuxNodes(0), _parentAuxNodeIndex(-1), _parentPixNodeIndex(-1), 
		_pathCheckTimeStamp(-1), _pathCheckValidPath(false), _excess(0)
	{
		_auxNodeInfo = new BKGCAuxNodeInfo[NUM_CLIQUES_PER_NODE];
		_nodeType = NODE_TYPE_PIX;
		_active = false;
	};

	~BKGCPixNode(void)
	{
		delete[] _auxNodeInfo;
	}
	void DeclareActives(std::list<BKGCNode*> *activeNodesList, int aux_index);
	bool HasPathToRoot();

public:
	DEBUG_STATEMENT(int _x);
	DEBUG_STATEMENT(int _y);

	BKGCAuxNodeInfo* _auxNodeInfo;

	float _excess;
	int _pathCheckTimeStamp;
	bool _pathCheckValidPath;

	uchar _numAuxNodes;
	char _parentAuxNodeIndex;
	char _parentPixNodeIndex;	
};

class BKGCPixNodeInfo
{
public:
	BKGCPixNodeInfo() : _numTightConstraints(0),_tightEdgeIntersectionCode(CODE_ALL_ONES_WORD),_tightEdgeUnionCode(0){}
	BKGCPixNode* _node;
	uchar _cliqueIndexInNode;
	CODE_TYPE _numTightConstraints;
	CODE_TYPE _tightEdgeIntersectionCode;
	CODE_TYPE _tightEdgeUnionCode;
};

class BKGCAuxNode : public BKGCNode
{
public:
	BKGCAuxNode(void)
		:_sinkLabeledNodesCode(0),_srcLabeledNodesCode(0), _saturatedEdgeCode(0),
		_parentIndex(-1), _childrenCode(0)
	{
		_nodeInfo = new BKGCPixNodeInfo[GG_CLIQUE_SIZE];
		_constraintSlacks = new float[NUM_CONSTRAINTS];
		memset(_constraintSlacks, 0, NUM_CONSTRAINTS*sizeof(float));
		_nodeType = NODE_TYPE_AUX;
		_active = false;
	};

	~BKGCAuxNode(void)
	{
		delete[] _nodeInfo;
		delete[] _constraintSlacks;
	}

	void DeclareActives(std::list<BKGCNode*> *activeNodesList);
	inline void InsertNode(BKGCPixNode* node, int node_index)
	{
		BKGCAuxNodeInfo& aux_node_info = node->_auxNodeInfo[node->_numAuxNodes++];
		aux_node_info._auxNode = this;
		aux_node_info._pixNodeIndexInAuxNode = node_index;
		aux_node_info._childrenCode = 0;

		BKGCPixNodeInfo& pix_node_info = _nodeInfo[node_index];
		pix_node_info._node = node;
		pix_node_info._cliqueIndexInNode = node->_numAuxNodes-1;
	}
	void AugmentFlow(int sending_node, int receiving_node, float flow);
	void RecalculateTightEdgeCodes(int index);
	void RecalculateSaturatedEdgeCode();
	void Reparametrize(double& flow);
	float FindMinContainingSlack(int node_index);
	float FindMinNonContainingSlack(int node_index);

public:
	DEBUG_STATEMENT(int _x);
	DEBUG_STATEMENT(int _y);

	int _index;

	BKGCPixNodeInfo* _nodeInfo;
	float* _constraintSlacks;

	CODE_TYPE _sinkLabeledNodesCode;
	CODE_TYPE _srcLabeledNodesCode;
	CODE_TYPE _saturatedEdgeCode;
	CODE_TYPE _childrenCode;
	char _parentIndex;
};
