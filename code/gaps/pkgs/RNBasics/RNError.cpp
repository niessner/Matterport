/* Source file for GAPS error utility  */



/* Include files */

#include "RNBasics.h"



/* Public variable */

static FILE *RNerror_file = stderr;
static int RNerror_level = 0;



int RNInitError() 
{
    /* Return OK status */
    return TRUE;
}



void RNStopError()
{
}



void RNSetErrorFile(FILE *fp)
{
    // Set error file (or NULL if want no errors printed)
    RNerror_file = fp;
}



void RNSetErrorLevel(int level)
{
    // Set error level 
    RNerror_level = level;
}



void RNAbort(const char *fmt, ...)
{
#if (0) && (RN_OS == RN_WINDOWS)
    /* Display error message */
    char outbuf[256];
    outbuf[0] = '\0';
    if (fmt) {
		char buf[256];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);
		strcat(outbuf, buf);
        va_end(args);
        strcat(outbuf, "\n");
    }
    MessageBox(NULL, outbuf, "GAPS FATAL ERROR -- ABORTING", MB_ICONSTOP | MB_OK);
#endif

    /* Print error message */
    if (RNerror_file) {
        perror("GAPS FATAL ERROR.  ABORTING:");
	if (fmt) {
	    va_list args;
	    va_start(args, fmt);
	    vfprintf(RNerror_file, fmt, args);
	    va_end(args);
	    fprintf(RNerror_file, "\n");
	}
	fflush(RNerror_file);
    }

    /* Abort */
    fflush(stdout);
    fflush(stderr);
    abort();
}
    


void RNFail(const char *fmt, ...)
{
#if (0) && (RN_OS == RN_WINDOWS)
    /* Display error message */
    char outbuf[256];
    outbuf[0] = '\0';
    if (fmt) {
		char buf[256];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);
		strcat(outbuf, buf);
        va_end(args);
        strcat(outbuf, "\n");
    }
    MessageBox(NULL, outbuf, "GAPS ERROR", MB_ICONSTOP | MB_OK);
#endif

    /* Print error message */
    if (RNerror_file) {
        perror("GAPS ERROR");
        if (fmt) {
	    va_list args;
	    va_start(args, fmt);
	    vfprintf(RNerror_file, fmt, args);
	    va_end(args);
	    fprintf(RNerror_file, "\n");
	}
    }
}



void RNWarning(const char *fmt, ...)
{
    /* Check error level */
    if (RNerror_level > 0) return;

#if (0) && (RN_OS == RN_WINDOWS)
    /* Display error message */
    char outbuf[256];
    outbuf[0] = '\0';
    if (fmt) {
		char buf[256];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);
		strcat(outbuf, buf);
        va_end(args);
        strcat(outbuf, "\n");
    }
    MessageBox(NULL, outbuf, "GAPS WARNING", MB_OK);
#endif
    
    /* Print error message */
    if (RNerror_file) {
        perror("GAPS WARNING");
        if (fmt) {
	    va_list args;
	    va_start(args, fmt);
	    vfprintf(RNerror_file, fmt, args);
	    va_end(args);
	    fprintf(RNerror_file, "\n");
	}
    }
}
    




