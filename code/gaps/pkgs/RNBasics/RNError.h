/* Include file for GAPS error utility */



/* Initialization functions */

int RNInitError();
void RNStopError();



/* Error file functions */

void RNSetErrorFile(FILE *fp);
void RNSetErrorLevel(int level);



/* Error reporting functions */

void RNAbort(const char *fmt, ...);
void RNFail(const char *fmt, ...);
void RNWarning(const char *fmt, ...);



/* Define my own assert function */

#ifdef assert
#  undef assert
#endif

#ifdef NDEBUG
#define assert(__x__)
#else
#define assert(__x__) \
  if (!(__x__)) RNAbort("Assertion error %s at line %d in file %s", #__x__, __LINE__, __FILE__)
#endif



