/* Source file for the R2 curve class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int 
R2InitCurve()
{
    /* Return success */
    return TRUE;
}



void 
R2StopCurve()
{
}



R2Curve::
R2Curve(void)
{
}



R2Curve::
~R2Curve(void)
{
}



const RNBoolean R2Curve::
IsCurve(void) const
{
    // All curve shapes are curves
    return TRUE;
}




