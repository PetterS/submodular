// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>
#include "opencv2/opencv.hpp"
#include <vector>
#include <list>
#include <iostream>
#pragma warning (disable : 4996)

#include "util/DebugUtils.h"
#include "util/GadgetGraphConfig.h"
#include "apgc/APGCConfig.h"
#include "prgc/PRGCConfig.h"

static const int BINARY_0_COLOR = 0;
static const int BINARY_1_COLOR = 255;
static const int BINARY_UNLABEL_COLOR = 150;

extern int CLIQUE_WIDTH;
extern int CLIQUE_HEIGHT;
extern int CLIQUE_SIZE;
extern int IMG_WIDTH;
extern int IMG_HEIGHT;
