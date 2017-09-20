/* Source file for the GAPS diad class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Diad R2xy_diad(R2Vector(1.0, 0.0), R2Vector(0.0, 1.0));



/* Public functions */

int R2InitDiad()
{
    /* Return success */
    return TRUE;
}



void R2StopDiad()
{
}



R2Diad::
R2Diad(void)
{
}



R2Diad::
R2Diad(const R2Diad& diad)
{
    axis[0] = diad.axis[0]; 
    axis[1] = diad.axis[1]; 
}



R2Diad::
R2Diad(const R2Vector& xaxis, const R2Vector& yaxis)
{
    // Just checking ...
    assert(xaxis.IsNormalized());
    assert(yaxis.IsNormalized());
    assert(R2Perpendicular(xaxis, yaxis));

    // Assign axes
    axis[0] = xaxis;
    axis[1] = yaxis;
}



const R3Matrix R2Diad::
Matrix(void) const
{
    // Return change of basis matrix (std -> diad)
    return R3Matrix(axis[0].X(), axis[1].X(), 0.0,
		    axis[0].Y(), axis[1].Y(), 0.0,
		    0.0,         0.0,         1.0);
}



const R3Matrix R2Diad::
InverseMatrix(void) const
{
    // Return change of basis matrix (diad -> std)
    return R3Matrix(axis[0].X(), axis[0].Y(), 0.0,
		    axis[1].X(), axis[1].Y(), 0.0,
		    0.0,         0.0,         1.0);
}



void R2Diad::
Normalize(void)
{
    // Normalize each axis 
    axis[0].Normalize();
    axis[1].Normalize();
}



void R2Diad:: 
Rotate(RNAngle radians)
{
    // Mirror each axis 
    axis[0].Rotate(radians);
    axis[1].Rotate(radians);
}



void R2Diad::
Mirror(const R2Line& line)
{
    // Mirror each axis 
    axis[0].Mirror(line);
    axis[1].Mirror(line);
}



void R2Diad::
Transform(const R2Transformation& transformation)
{
    // Transform each axis 
    axis[0].Transform(transformation);
    axis[1].Transform(transformation);
}



void R2Diad::
InverseTransform(const R2Transformation& transformation)
{
    // Transform each axis 
    axis[0].InverseTransform(transformation);
    axis[1].InverseTransform(transformation);
}









