#pragma once

extern int APGC_NUM_CONSTRAINTS_PER_CLIQUE;

class APGCNode;
typedef std::list<APGCNode*> APGCNodeList;

class APGCPixNode;
typedef std::list<APGCPixNode*> APGCPNodeList;

class APGCAuxNode;
extern std::list<std::pair<APGCAuxNode*, int>> TIGHT_CONSTRAINTS;

#include "apgc/APGCDataStructures.h"
