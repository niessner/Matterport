/* Source file for the GAPS halfspace class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Halfspace R2null_halfspace(0.0, 0.0, 0.0);
const R2Halfspace R2posx_halfspace(1.0, 0.0, 0.0);
const R2Halfspace R2posy_halfspace(0.0, 1.0, 0.0);
const R2Halfspace R2negx_halfspace(-1.0, 0.0, 0.0);
const R2Halfspace R2negy_halfspace(0.0, -1.0, 0.0);



/* Public functions */

int 
R2InitHalfspace()
{
    /* Return success */
    return TRUE;
}



void 
R2StopHalfspace()
{
}



R2Halfspace::
R2Halfspace (void)
{
}



R2Halfspace::
R2Halfspace (const R2Halfspace& halfspace)
    : line(halfspace.line)
{
}



R2Halfspace::
R2Halfspace(RNScalar a, RNScalar b, RNScalar c)
    : line(a, b, c)
{
}



R2Halfspace::
R2Halfspace(const RNScalar a[3])
    : line(a)
{
}



R2Halfspace::
R2Halfspace(const R2Point& point, const R2Vector& normal)
    : line(point, R2Vector(-normal.Y(), normal.X()))
{
}



R2Halfspace::
R2Halfspace(const R2Point& point1, const R2Point& point2)
    : line(point1, point2)
{
}



R2Halfspace::
R2Halfspace(const R2Line& line, int )
    : line(line)
{
}



void R2Halfspace::
Mirror (const R2Line& line)
{
    // Mirror halfspace with respect to line
    this->line.Mirror(line);
}



void R2Halfspace::
Transform (const R2Transformation& transformation)
{
    // Transform halfspace
    line.Transform(transformation);
}



void R2Halfspace::
InverseTransform (const R2Transformation& transformation)
{
    // Transform halfspace
    line.InverseTransform(transformation);
}




