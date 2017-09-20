/* Source file for the GAPS triad class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Triad R3xyz_triad(
    R3Vector(1.0, 0.0, 0.0),
    R3Vector(0.0, 1.0, 0.0),
    R3Vector(0.0, 0.0, 1.0));



/* Public functions */

int R3InitTriad()
{
    /* Return success */
    return TRUE;
}



void R3StopTriad()
{
}



R3Triad::
R3Triad(void)
{
}



R3Triad::
R3Triad(const R3Triad& triad)
{
    axis[0] = triad.axis[0]; 
    axis[1] = triad.axis[1]; 
    axis[2] = triad.axis[2]; 
}



R3Triad::
R3Triad(const R3Vector& xaxis, const R3Vector& yaxis, const R3Vector& zaxis)
{
    // Just checking ...
    assert(xaxis.IsNormalized());
    assert(yaxis.IsNormalized());
    assert(zaxis.IsNormalized());
    assert(R3Perpendicular(xaxis, yaxis));
    assert(R3Contains(xaxis % yaxis, zaxis));

    // Assign axes
    axis[0] = xaxis;
    axis[1] = yaxis;
    axis[2] = zaxis;
}



R3Triad::
R3Triad(const R3Vector& towards, const R3Vector& up)
{
    axis[2] = -towards;
    axis[2].Normalize();
    axis[0] = up % axis[2];
    axis[0].Normalize();
    axis[1] = axis[2] % axis[0];
    axis[1].Normalize();
}



const R4Matrix R3Triad::
Matrix(void) const
{
    // Return change of basis matrix (std -> triad)
    return R4Matrix(axis[0].X(), axis[1].X(), axis[2].X(), 0.0,
		    axis[0].Y(), axis[1].Y(), axis[2].Y(), 0.0,
		    axis[0].Z(), axis[1].Z(), axis[2].Z(), 0.0,
		    0.0,         0.0,         0.0,         1.0);
}



const R4Matrix R3Triad::
InverseMatrix(void) const
{
    // Return change of basis matrix (triad -> std)
    return R4Matrix(axis[0].X(), axis[0].Y(), axis[0].Z(), 0.0,
		    axis[1].X(), axis[1].Y(), axis[1].Z(), 0.0,
		    axis[2].X(), axis[2].Y(), axis[2].Z(), 0.0,
		    0.0,         0.0,         0.0,         1.0);
}



void R3Triad::
Normalize(void)
{
    // Normalize each axis 
    axis[0].Normalize();
    axis[1].Normalize();
    axis[2].Normalize();
}



void R3Triad:: 
Rotate(RNAxis axis, RNAngle radians)
{
    // Rotate triad around an axis counterclockwise
    switch (axis) {
    case RN_XAXIS: 
	Rotate(R3posx_vector, radians); 
	break;

    case RN_YAXIS: 
	Rotate(R3posy_vector, radians); 
	break;

    case RN_ZAXIS: 
	Rotate(R3posz_vector, radians); 
	break;

    default: 
	RNWarning("Triad rotation around undefined axis");
	break;
    }
}



void R3Triad::
Rotate(const R3Vector& rotaxis, RNAngle radians)
{
    // Rotate each axis
    axis[0].Rotate(rotaxis, radians);
    axis[1].Rotate(rotaxis, radians);
    axis[2].Rotate(rotaxis, radians);
    axis[2].Normalize();
    axis[0] = axis[1] % axis[2];
    axis[0].Normalize();
    axis[1] = axis[2] % axis[0];
    axis[1].Normalize();
}



void R3Triad::
Rotate(const R3Vector& from, const R3Vector& to)
{
    // Rotate each axis
    RNAngle angle = R3InteriorAngle(from, to);
    R3Vector rotaxis = from % to;
    rotaxis.Normalize();
    Rotate(rotaxis, angle);
    axis[2].Normalize();
    axis[0] = axis[1] % axis[2];
    axis[0].Normalize();
    axis[1] = axis[2] % axis[0];
    axis[1].Normalize();
}



void R3Triad::
Mirror(const R3Plane& plane)
{
    // Mirror each axis 
    axis[0].Mirror(plane);
    axis[1].Mirror(plane);
    axis[2].Mirror(plane);
    axis[2].Normalize();
    axis[0] = axis[1] % axis[2];
    axis[0].Normalize();
    axis[1] = axis[2] % axis[0];
    axis[1].Normalize();
}



void R3Triad::
Transform(const R3Transformation& transformation)
{
    // Transform each axis 
    axis[0].Transform(transformation);
    axis[1].Transform(transformation);
    axis[2].Transform(transformation);
    axis[2].Normalize();
    axis[0] = axis[1] % axis[2];
    axis[0].Normalize();
    axis[1] = axis[2] % axis[0];
    axis[1].Normalize();
}



void R3Triad::
InverseTransform(const R3Transformation& transformation)
{
    // Transform each axis 
    axis[0].InverseTransform(transformation);
    axis[1].InverseTransform(transformation);
    axis[2].InverseTransform(transformation);
    axis[2].Normalize();
    axis[0] = axis[1] % axis[2];
    axis[0].Normalize();
    axis[1] = axis[2] % axis[0];
    axis[1].Normalize();
}









