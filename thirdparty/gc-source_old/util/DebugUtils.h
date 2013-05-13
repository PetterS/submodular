#pragma once

//#define CHECK_FOR_ERRORS
#ifdef _DEBUG
#define CHECK_FOR_ERRORS
#endif

#define XPRINTF(p,...)			{ \
	if((p) == 0x01) { \
	printf(__VA_ARGS__); \
	} else if((p) == 0x10) { \
	if(STATE_FILE) fprintf(STATE_FILE, __VA_ARGS__); \
	} else { \
	printf(__VA_ARGS__); \
	if(STATE_FILE) fprintf(STATE_FILE, __VA_ARGS__); \
								  } \
								}
#define PRINTF(...)				{ if(PRINT_ENABLED) { XPRINTF(3,__VA_ARGS__); fflush(STATE_FILE); } }
#define PRINTLN					PRINTF("\n")

#ifdef CHECK_FOR_ERRORS
#define DEBUG_COND(cond)		DEBUG_COND_FUNC((cond))
#define DEBUG_STATEMENT(X)		X
#define MY_ASSERT(X)			MY_ASSERT_FUNC(X)
#else
#define DEBUG_STATEMENT(X)
#define DEBUG_COND(cond)
#define MY_ASSERT(X)			
#endif

#pragma warning (disable : 4996)
extern FILE* STATE_FILE;
extern bool PRINT_ENABLED;
extern void MY_ASSERT_FUNC(bool cond);
extern void DEBUG_COND_FUNC(bool cond);
