#include "stdafx.h"
#include "DebugUtils.h"

void DEBUG_COND_FUNC(bool cond)
{
	if(cond)
		printf("debug condition hit\n");
}

void MY_ASSERT_FUNC(bool cond)
{
	if(!cond)
	{
		printf("\n\n assertion failure\n");
		if(STATE_FILE){
			fflush(STATE_FILE);
			fclose(STATE_FILE);
		}
		exit(1);
	}
}

FILE* STATE_FILE = NULL;
bool PRINT_ENABLED = true;

