/* Source file for GAPS basics module */



/* Include files */

#include "RNBasics.h"



/* Private variables */

static int RNbasics_active_count = 0;



int RNInitBasics(void)
{
    // Check whether are already initialized 
    if ((RNbasics_active_count++) > 0) return TRUE;

    // Initialize submodules 
    RNSeedRandomScalar();

    // Return OK status 
    return TRUE;
}



void RNStopBasics(void)
{
    // Check whether have been initialized 
    if ((--RNbasics_active_count) > 0) return;

    // Stop submodules 
    // ???
}




