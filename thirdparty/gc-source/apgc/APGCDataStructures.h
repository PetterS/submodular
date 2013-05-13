#pragma once

class APGCAuxNode;

class APGCNode
{
public:
	class APGCNode() :_parentPixNodeIndex(-1), 
		_pathCheckTimeStamp(-1), _pathCheckValidPath(false),_tree(TREE_NONE),
		_distanceLabel(DIST_DEFAULT),_active(false){}

	char _nodeType;
	bool _active;	
	int _pathCheckTimeStamp;
	bool _pathCheckValidPath;
	char _tree;
	char _parentPixNodeIndex;	
	DIST_TYPE _distanceLabel;
};

class APGCEdge
{
public:
	APGCAuxNode* _auxNode;
	char _sendingNodeIndex;
	char _receivingNodeIndex;
};

class APGCAuxNodeInfo
{
public:
	APGCAuxNode* _auxNode;
	uchar _pixNodeIndexInAuxNode;
	CODE_TYPE _childrenCode;
};

class APGCPixNode : public APGCNode
{
public:
	APGCPixNode(void) :_parentAuxNodeIndex(-1), _numAuxNodes(0), _excess(0) {
		_auxNodeInfo = new APGCAuxNodeInfo[NUM_CLIQUES_PER_NODE];
		_nodeType = NODE_TYPE_PIX;
	}

	bool HasPathToRoot();
	
public:
	DEBUG_STATEMENT(int _x);
	DEBUG_STATEMENT(int _y);
	DEBUG_STATEMENT(int _nodeIndex);

	APGCAuxNodeInfo* _auxNodeInfo;
	uchar _numAuxNodes;

	float _excess;
	char _parentAuxNodeIndex;
};

class APGCPixNodeInfo
{
public:
	APGCPixNodeInfo() : _numTightConstraints(0),_tightEdgeIntersectionCode(CODE_ALL_ONES_WORD),
		_tightEdgeUnionCode(0){}
	APGCPixNode* _node;
	uchar _cliqueIndexInNode;
	CODE_TYPE _numTightConstraints;
	CODE_TYPE _tightEdgeIntersectionCode;
	CODE_TYPE _tightEdgeUnionCode;
};

class APGCAuxNode : public APGCNode
{
public:
	APGCAuxNode(void)
		:_sinkLabeledNodesCode(0),_srcLabeledNodesCode(0), _saturatedEdgeCode(0),_childrenCode(0)
	{
		_nodeInfo = new APGCPixNodeInfo[GG_CLIQUE_SIZE];
		_constraintSlacks = new float[NUM_CONSTRAINTS];
		memset(_constraintSlacks, 0, NUM_CONSTRAINTS*sizeof(float));
		_nodeType = NODE_TYPE_AUX;
		_active = false;
	};
	~APGCAuxNode(){
		delete[] _nodeInfo;
		delete[] _constraintSlacks;
	}
	inline void InsertNode(APGCPixNode* node, int node_index)
	{
		APGCAuxNodeInfo& aux_node_info = node->_auxNodeInfo[node->_numAuxNodes++];
		aux_node_info._auxNode = this;
		aux_node_info._pixNodeIndexInAuxNode = node_index;
		aux_node_info._childrenCode = 0;

		APGCPixNodeInfo& pix_node_info = _nodeInfo[node_index];
		pix_node_info._node = node;
		pix_node_info._cliqueIndexInNode = node->_numAuxNodes-1;
	}
	bool HasPathToRoot();
	void AugmentFlow(int sending_node, int receiving_node, float flow);
	void DeclareConstraintNotTight(int constraint_index, CODE_TYPE constraint_labeling);
	void DeclareConstraintTight(int constraint_index, CODE_TYPE constraint_labeling);
	void RecalculateTightEdgeCodes(int index);
	void RecalculateSaturatedEdgeCode();
	void Reparametrize(double& flow);
	float FindMinContainingSlack(int node_index);
	float FindMinNonContainingSlack(int node_index);

public:
	DEBUG_STATEMENT(int _x);
	DEBUG_STATEMENT(int _y);

	int _index;

	APGCPixNodeInfo* _nodeInfo;
	float* _constraintSlacks;

	CODE_TYPE _sinkLabeledNodesCode;
	CODE_TYPE _srcLabeledNodesCode;
	CODE_TYPE _saturatedEdgeCode;
	CODE_TYPE _childrenCode;
};
