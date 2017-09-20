/* Source file for the R3 transformation class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public functions */

int 
R3InitTransformation()
{
    /* Return success */
    return TRUE;
}



void 
R3StopTransformation()
{
}



R3Transformation::
~R3Transformation(void)
{
}



const RNBoolean R3Transformation::
IsIdentity(void) const
{
    // Return whether transformation is identity (default is FALSE)
    return FALSE;
}



const RNBoolean R3Transformation::
IsMirrored(void) const
{
    // Return whether transformation is mirrored (default is FALSE)
    return FALSE;
}



const RNBoolean R3Transformation::
IsLinear(void) const
{
    // Return whether transformation is linear (default is FALSE)
    return FALSE;
}



const RNBoolean R3Transformation::
IsAffine(void) const
{
    // Return whether transformation is affine (default is FALSE)
    return FALSE;
}



const RNBoolean R3Transformation::
IsIsotropic(void) const
{
    // Return whether transformation is isotropic (default is FALSE)
    return FALSE;
}



const RNScalar R3Transformation::
ScaleFactor(void) const
{
    // Return identity scaling
    return 1.0;
}





