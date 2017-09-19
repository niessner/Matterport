/* Source file for the R3 surface class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public functions */

int 
R3InitSurface()
{
    /* Return success */
    return TRUE;
}



void 
R3StopSurface()
{
}



R3Surface::
R3Surface(void)
{
}



R3Surface::
~R3Surface(void)
{
}



const RNBoolean R3Surface::
IsSurface(void) const
{
    // All surface shapes are surfaces
    return TRUE;
}




