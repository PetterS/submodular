#pragma once

typedef int DIST_TYPE;
static const DIST_TYPE DIST_DEFAULT = -1;
static const int TERM_DIST = 0;

enum LABEL {LABEL_SOURCE_0 = 0, LABEL_SINK_1 = 1, LABEL_NONE=2};

#define TREE_NONE 'n'
#define TREE_SRC 's'
#define TREE_SINK 't'

#define NODE_TYPE_A_AUX 'a'
#define NODE_TYPE_B_AUX 'b'
#define NODE_TYPE_PIX   'p'
#define NODE_TYPE_AUX   'a'

extern float EPSILON;

#define EQUAL(f1, f2) (fabs(f1-f2) <= EPSILON)
#define LESS_THAN(f1, f2) ((f2-f1) > EPSILON)
#define GREATER_THAN(f1, f2) ((f1-f2) > EPSILON)

extern bool GADGET_GRAPH_DEBUG_VERBOSE;

extern int NUM_NODES;

extern int GG_CLIQUE_SIZE;
extern int D_GG_CLIQUE_SIZE;

extern int NUM_CLIQUES_PER_NODE;

typedef int CODE_TYPE;
extern CODE_TYPE CODE_ALL_ONES_WORD;
extern CODE_TYPE CODE_ALL_ONES_DWORD;
extern CODE_TYPE* INDEX_CODE_ONE;
extern CODE_TYPE* INDEX_CODE_ZERO;

#define IS_BIT_1(X,INDEX) (((X) & INDEX_CODE_ONE[(INDEX)]) != 0)
#define IS_BIT_0(X,INDEX) (((X) & INDEX_CODE_ONE[(INDEX)]) == 0)
#define SET_BIT_0(X,INDEX) ((X) &= INDEX_CODE_ZERO[(INDEX)])

class ConstraintInfo;
extern int NUM_CONSTRAINTS;
extern ConstraintInfo* CONSTRAINT_INFO;

#define CONSTRAINT_LABELING(X) CONSTRAINT_INFO[(X)]._labeling

extern int NUM_CONTAINING_CONSTRAINTS;
extern int** CONTAINING_CONSTRAINTS;
extern int** NON_CONTAINING_CONSTRAINTS;

extern __int64 ITER_NUM;
extern int TIME;

#include "util/FirstNonZero.h"
static inline uchar FIRST_NON_ZERO(CODE_TYPE code)
{
	int index = (code > USHRT_MAX) ? (FIRST_NON_ZERO_WORD[(code >> 16)] + 16) : FIRST_NON_ZERO_WORD[code];
	return index;
}

static char* D_BINARY(CODE_TYPE val)
{
	static const int size = sizeof(CODE_TYPE)*8;
	static char out[size+1];
	out[D_GG_CLIQUE_SIZE] = '\0';
	for(int i=0; i<D_GG_CLIQUE_SIZE; ++i)
	{
		if(IS_BIT_1(val, i))
			out[D_GG_CLIQUE_SIZE-1-i] = '1';
		else
			out[D_GG_CLIQUE_SIZE-1-i] = '0';
	}
	return out;
}

static char* BINARY(CODE_TYPE val)
{
	static const int size = sizeof(CODE_TYPE)*8;
	static char out[size+1];
	out[GG_CLIQUE_SIZE] = '\0';
	for(int i=0; i<GG_CLIQUE_SIZE; ++i)
	{
		if(IS_BIT_1(val, i))
			out[GG_CLIQUE_SIZE-1-i] = '1';
		else
			out[GG_CLIQUE_SIZE-1-i] = '0';
	}
	return out;
}

#define DEBUG_ITER(X)	DEBUG_COND(ITER_NUM == (X))

class ConstraintInfo
{
public:
	ConstraintInfo(){ _containedNodes = new int[GG_CLIQUE_SIZE]; }
	~ConstraintInfo() { delete[] _containedNodes; }

	int* _containedNodes;
	CODE_TYPE _labeling;
};

class FlowConstraint
{
public:
	FlowConstraint(int clique_index=-1, int constraint_index=-1) :
	  _cliqueIndex(clique_index), _constraintIndex(constraint_index){}
	  int _cliqueIndex;
	  int _constraintIndex;
	  bool operator < (const FlowConstraint& rhs) const
	  {
		  if(_cliqueIndex < rhs._cliqueIndex)
			  return true;
		  if(_cliqueIndex > rhs._cliqueIndex)
			  return false;
		  if(_constraintIndex < rhs._constraintIndex)
			  return true;
		  if(_constraintIndex > rhs._constraintIndex)
			  return false;
		  return false;
	  }
	  inline bool operator == (const FlowConstraint& rhs) const
	  {
		  return ((_cliqueIndex == rhs._cliqueIndex) && (_constraintIndex == rhs._constraintIndex));
	  }
};


