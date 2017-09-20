/* Source file for the R3 curve class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public functions */

int 
R3InitCurve()
{
    /* Return success */
    return TRUE;
}



void 
R3StopCurve()
{
}



R3Curve::
R3Curve(void)
{
}



R3Curve::
~R3Curve(void)
{
}



const R3Point R3Curve::
StartPosition(void) const
{
    // Return position at start of curve
    return PointPosition(StartParameter());
}



const R3Point R3Curve::
EndPosition(void) const
{
    // Return position at end of curve
    return PointPosition(EndParameter());
}



R3Vector R3Curve::
PointDirection(RNScalar u) const
{
    // Return tangent direction of curve at parametric value u
    R3Vector direction = PointDerivative(u);
    direction.Normalize();
    return direction;
}



const RNBoolean R3Curve::
IsCurve(void) const
{
    // All curve shapes are curves
    return TRUE;
}




