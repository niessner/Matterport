/* Source file for GAPS shapes module */



/* Include files */

#include "R2Shapes.h"



/* Private variables */

static int R2shapes_active_count = 0;



int R2InitShapes(void)
{
    // Check whether are already initialized 
    if ((R2shapes_active_count++) > 0) return TRUE;

    // Initialize dependencies
    if (!RNInitBasics()) return FALSE;

    // Initialize submodules 
    if (!R2InitCircle()) return FALSE;

    // return OK status 
    return TRUE;
}



void R2StopShapes(void)
{
    // Check whether have been initialized 
    if ((--R2shapes_active_count) > 0) return;

    // Stop submodules 
    R2StopCircle();

    // Stop dependencies
    RNStopBasics();
}






