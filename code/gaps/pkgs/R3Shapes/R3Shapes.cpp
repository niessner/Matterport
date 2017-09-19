/* Source file for GAPS shapes module */



/* Include files */

#include "R3Shapes.h"



/* Public 3DR variables -- must be set by the application ??? */

#if (RN_3D_GRFX == RN_3DR)
    R3dHandle_t R3dr_rc = NULL;
    G3dHandle_t R3dr_gc = NULL;
#endif



/* Private variables */

static int R3shapes_active_count = 0;



int R3InitShapes(void)
{
    // Check whether are already initialized 
    if ((R3shapes_active_count++) > 0) return TRUE;

    // Initialize dependencies
    if (!R2InitShapes()) return FALSE;

    // Initialize submodules 
    if (!R3InitCircle()) return FALSE;

    // return OK status 
    return TRUE;
}



void R3StopShapes(void)
{
    // Check whether have been initialized 
    if ((--R3shapes_active_count) > 0) return;

    // Stop submodules 
    R3StopCircle();

    // Stop dependencies
    R2StopShapes();
}




