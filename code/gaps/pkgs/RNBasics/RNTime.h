/* Include file for GAPS wallclock time class */



/* Initialization functions */

int RNInitTime();
void RNStopTime();



/* Class definition */

class RNTime /* : public RNBase */ {
    public:
        // Constructor functions
        RNTime(void);
	RNTime(const RNTime& tm);

        // Relationship functions/operators
        RNScalar Elapsed(const RNTime& tm) const;
        RNScalar Elapsed(void) const;

        // Arithmetic operators
	RNScalar operator-(const RNTime& tm) const;

        // Manipulation functions/operators
	void Read(void);

    private:
#       if (RN_OS == RN_WINDOWS)
            LARGE_INTEGER timefreq;
            LARGE_INTEGER timevalue;
#       elif (RN_OS == OLD_RN_WINDOWS)
	    DWORD timevalue;
#       else
	    struct timeval timevalue;
#       endif 
};



/* Public functions */

RNTime RNCurrentTime(void);
void RNSleep(RNScalar seconds);



/* Inline functions */

inline void RNTime::
Read (void) 
{
    /* Read the current time */
#   if (RN_OS == RN_WINDOWS)
        QueryPerformanceCounter(&timevalue);
#   elif (RN_OS == OLD_RN_WINDOWS)
	timevalue = GetTickCount();
#   else
	gettimeofday(&timevalue, NULL);
#   endif 
}



