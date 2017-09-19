/* Include file for GAPS basics */



/* Initialization functions */

int RNInitBase();
void RNStopBase();



/* Class definition */

class RNBase {};



/* TRUE/FALSE constant definitions */

typedef int RNBoolean;
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif



/* Programming convenience macros */

#define brcase break; case



/* Standard types */

typedef int RNMark;



/* Global variables */

extern RNMark RNmark;



/* Useful functions */

extern void RNBreakDebug(void);
