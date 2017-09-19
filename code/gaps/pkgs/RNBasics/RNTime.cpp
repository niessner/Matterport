/* Source file for GAPS real wallclock time class  */



/* Include files */

#include "RNBasics.h"



int RNInitTime() 
{
    /* Return OK status */
    return TRUE;
}



void RNStopTime()
{
}



RNTime::
RNTime(void)
{
#   if (RN_OS == RN_WINDOWS)
        QueryPerformanceFrequency(&timefreq);
#   endif
}



RNTime::
RNTime(const RNTime& tm)
    : timevalue(tm.timevalue)
{
#   if (RN_OS == RN_WINDOWS)
        timefreq = timefreq;
#   endif
}



RNScalar RNTime::
operator- (const RNTime& tm) const
{
    /* Return the difference between this and tm times (in seconds) */
#   if (RN_OS == RN_WINDOWS)
        return ((RNScalar) this->timevalue.QuadPart - (RNScalar) tm.timevalue.QuadPart) / ((RNScalar) this->timefreq.QuadPart);
#   elif (RN_OS == OLD_RN_WINDOWS)
	return (RNScalar) (this->timevalue - tm.timevalue) / 1.0E3;
#   else
	return (RNScalar) (this->timevalue.tv_sec - tm.timevalue.tv_sec + 
		1.0E-6F * (this->timevalue.tv_usec - tm.timevalue.tv_usec));
#   endif 
}



RNScalar RNTime::
Elapsed (const RNTime& tm) const
{
    // Return number of seconds elapsed between times
    RNScalar delta = *this - tm;
    assert(delta >= 0.0);
    return delta;
}



RNScalar RNTime::
Elapsed (void) const
{
    // Return number of seconds elapsed since time
    RNTime tm; tm.Read();
    return tm.Elapsed(*this);
}



RNTime 
RNCurrentTime(void)
{
    RNTime tm;
    tm.Read();
    return tm;
}




#if (RN_OS != RN_WINDOWS)
#   include <unistd.h>
#endif

void 
RNSleep(RNScalar seconds)
{
#if (RN_OS == RN_IRIX)
    sginap((long) (seconds * CLK_TCK));
#elif (RN_OS == RN_WINDOWS)
    Sleep((unsigned long) (1000 * seconds));
#elif (RN_OS == RN_LINUX)
    usleep((unsigned long) (1000000 * seconds));
#else
    RNAbort("Not implemented");
#endif
}



