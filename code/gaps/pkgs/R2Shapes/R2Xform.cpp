/* Source file for the R2 transformation class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int 
R2InitTransformation()
{
    /* Return success */
    return TRUE;
}



void 
R2StopTransformation()
{
}



R2Transformation::
~R2Transformation(void)
{
}



const RNBoolean R2Transformation::
IsMirrored(void) const
{
    // Return whether transformation is mirrored (default is FALSE)
    return FALSE;
}



const RNBoolean R2Transformation::
IsLinear(void) const
{
    // Return whether transformation is linear (default is FALSE)
    return FALSE;
}



const RNBoolean R2Transformation::
IsAffine(void) const
{
    // Return whether transformation is affine (default is FALSE)
    return FALSE;
}



const RNBoolean R2Transformation::
IsIsotropic(void) const
{
    // Return whether transformation is isotropic (default is FALSE)
    return FALSE;
}




