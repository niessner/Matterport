/* Source file for the R3 solid class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public functions */

int 
R3InitSolid()
{
    /* Return success */
    return TRUE;
}



void 
R3StopSolid()
{
}



R3Solid::
R3Solid(void)
{
}



R3Solid::
~R3Solid(void)
{
}



const RNBoolean R3Solid::
IsSolid(void) const
{
    // All solid shapes are solids
    return TRUE;
}




