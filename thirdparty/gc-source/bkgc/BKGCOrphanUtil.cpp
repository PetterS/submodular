#include "StdAfx.h"
#include "bkgc/BKGCOrphanUtil.h"
#include "bkgc/BKGCDump.h"

BKGCOrphanUtil::BKGCOrphanUtil(void)
{
	_currentNodeList = new std::list<BKGCNode*>();
	_nextNodeList = new std::list<BKGCNode*>();
}

BKGCOrphanUtil::~BKGCOrphanUtil(void)
{
	delete _currentNodeList;
	delete _nextNodeList;
}

void BKGCOrphanUtil::DeclareTreeOrphan(uchar tree_type, BKGCNode* root, std::list<BKGCNode*> *activeNodesList)
{
	//DEBUG_COND(ITER_NUM == 13);

	_currentNodeList->clear();
	_nextNodeList->clear();

	if(root->_nodeType == NODE_TYPE_AUX){
		((BKGCAuxNode*)root)->_parentIndex = -1;
		root->_tree = TREE_NONE;
	}else{
		BKGCPixNode* pix_node = (BKGCPixNode*)root;
		if(pix_node->_parentAuxNodeIndex != -1){
			BKGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[pix_node->_parentAuxNodeIndex];
			BKGCAuxNode* aux_node = aux_info._auxNode;
			if(pix_node->_parentPixNodeIndex == -1)
				aux_node->_childrenCode &= (1 << aux_info._pixNodeIndexInAuxNode) ^ CODE_ALL_ONES_WORD;
			else{
				BKGCPixNodeInfo& parent_node_info = aux_node->_nodeInfo[pix_node->_parentPixNodeIndex];
				BKGCAuxNodeInfo& parent_aux_info = parent_node_info._node->_auxNodeInfo[parent_node_info._cliqueIndexInNode];
				parent_aux_info._childrenCode &= (1<<aux_info._pixNodeIndexInAuxNode) ^ CODE_ALL_ONES_WORD;
			}		
		}
		if(tree_type == TREE_SRC)
			DeclareSrcPixNodeOrphan(pix_node);
		else
			DeclareSinkPixNodeOrphan(pix_node);
	}
	DEBUG_STATEMENT(PRINTF("Could not find path - Declaring orphan subtree rooted at "); BKGCDump::PrintNode(root); PRINTF("\n"));
	_currentNodeList->push_back(root);
	while(true)
	{
		while(!_currentNodeList->empty())
		{
			BKGCNode* graph_node = _currentNodeList->front(); _currentNodeList->pop_front();
			if(graph_node->_nodeType == NODE_TYPE_AUX)
			{
				BKGCAuxNode* aux_node = (BKGCAuxNode*)graph_node;
				DEBUG_STATEMENT(BKGCDump::PrintNode(aux_node); PRINTF(" declared orphan\n"));
				aux_node->DeclareActives(activeNodesList);

				//some of the nodes may be waiting for taking path from me
				//since I am now becoming orphan, there may be a possibility to take p2p path
				if(tree_type == TREE_SINK && aux_node->_active)
				{
					int nodes_to_declare_active = aux_node->_saturatedEdgeCode & 
												  aux_node->_sinkLabeledNodesCode & 
												  (aux_node->_childrenCode ^ CODE_ALL_ONES_WORD);
					while(nodes_to_declare_active != 0)
					{
						int i = FIRST_NON_ZERO_WORD[nodes_to_declare_active];
						nodes_to_declare_active &= (1 << i) ^ CODE_ALL_ONES_WORD;
						BKGCPixNode* pix_node = aux_node->_nodeInfo[i]._node;
						if(!pix_node->_active){
							pix_node->_active = true;
							activeNodesList->push_back(pix_node);
						}
					}
				}

				while(aux_node->_childrenCode != 0)
				{
					int i = FIRST_NON_ZERO_WORD[aux_node->_childrenCode];
					int i_code_inverse = (1 << i) ^ CODE_ALL_ONES_WORD;
					aux_node->_childrenCode &= i_code_inverse;
					BKGCPixNode* pix_node = aux_node->_nodeInfo[i]._node;
					if(tree_type == TREE_SRC)
						DeclareSrcPixNodeOrphan(pix_node);
					else
						DeclareSinkPixNodeOrphan(pix_node);
					_nextNodeList->push_back(pix_node);
				}
			}else
			{
				BKGCPixNode* pix_node = (BKGCPixNode*)graph_node;
				DEBUG_STATEMENT(BKGCDump::PrintNode(pix_node); PRINTF(" declared orphan\n"));
				//DEBUG_COND(ITER_NUM == 28 && pix_node->_x == 3 && pix_node->_y == 4);
				for(int j=0; j<pix_node->_numAuxNodes; ++j)
				{
					BKGCAuxNodeInfo& aux_info = pix_node->_auxNodeInfo[j];
					BKGCAuxNode* aux_node = aux_info._auxNode;
					int my_index = aux_info._pixNodeIndexInAuxNode;
					BKGCPixNodeInfo& my_info = aux_node->_nodeInfo[my_index];

					if(aux_node->_parentIndex == my_index){
						aux_node->_tree = TREE_NONE;
						aux_node->_parentIndex = -1;
						_nextNodeList->push_front(aux_node);
					}					
					pix_node->DeclareActives(activeNodesList, j);
					while(aux_info._childrenCode != 0)
					{
						int k = FIRST_NON_ZERO_WORD[aux_info._childrenCode];
						aux_info._childrenCode &= (1 << k) ^ CODE_ALL_ONES_WORD;
						BKGCPixNode* pix_node = aux_node->_nodeInfo[k]._node;
						MY_ASSERT(pix_node->_parentPixNodeIndex == aux_info._pixNodeIndexInAuxNode);
						if(tree_type == TREE_SRC)
							DeclareSrcPixNodeOrphan(pix_node);
						else
							DeclareSinkPixNodeOrphan(pix_node);
						_nextNodeList->push_back(pix_node);
					}
				}
			}
		}
		if(_nextNodeList->empty())
			break;
		std::list<BKGCNode*> *temp = _currentNodeList;
		_currentNodeList = _nextNodeList;
		_nextNodeList = temp;
	}
}