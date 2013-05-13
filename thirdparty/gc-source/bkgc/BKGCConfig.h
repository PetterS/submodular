#pragma once

class BKGCNode;
typedef std::list<BKGCNode*> BKGCNodeList;

class BKGCPixNode;
typedef std::list<BKGCPixNode*> BKGCPNodeList;

class BKGCAuxNode;

extern std::list<std::pair<BKGCAuxNode*, int>>* BKGC_TIGHT_CONSTRAINTS;
extern std::list<BKGCNode*>* BKGC_ACTIVE_NODES;

#include "bkgc/BKGCDataStructures.h"
