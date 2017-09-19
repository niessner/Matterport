/* Source file for the R2 solid class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int 
R2InitSolid()
{
    /* Return success */
    return TRUE;
}



void 
R2StopSolid()
{
}



R2Solid::
R2Solid(void)
{
}



R2Solid::
~R2Solid(void)
{
}



const RNBoolean R2Solid::
IsSolid(void) const
{
    // All solid shapes are solids
    return TRUE;
}




