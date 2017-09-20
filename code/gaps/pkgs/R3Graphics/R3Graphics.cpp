/* Source file for GAPS graphics module */



/* Include files */

#include "R3Graphics.h"



/* Private variables */

static int R3graphics_active_count = 0;



int R3InitGraphics(void)
{
    // Check whether are already initialized 
    if ((R3graphics_active_count++) > 0) return TRUE;

    // Initialize dependencies
    if (!R3InitShapes()) return FALSE;

    // Initialize submodules 
    if (!R3InitMaterial()) return FALSE;

    // return OK status 
    return TRUE;
}



void R3StopGraphics(void)
{
    // Check whether have been initialized 
    if ((--R3graphics_active_count) > 0) return;

    // Stop submodules 
    R3StopMaterial();

    // Stop dependencies
    R3StopShapes();
}






