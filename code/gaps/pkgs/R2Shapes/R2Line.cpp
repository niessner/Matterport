/* Source file for the GAPS line class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Line R2null_line(0.0, 0.0, 0.0);
const R2Line R2posx_line(0.0, -1.0, 0.0);
const R2Line R2posy_line(1.0, 0.0, 0.0);
const R2Line R2negx_line(0.0, 1.0, 0.0);
const R2Line R2negy_line(-1.0, 0.0, 0.0);



/* Public functions */

int R2InitLine()
{
    /* Return success */
    return TRUE;
}



void R2StopLine()
{
}



R2Line::
R2Line(void)
{
}



R2Line::
R2Line(const R2Line& line)
    : vector(line.vector),
      normal(line.normal),
      c(line.c)
{
}



R2Line::
R2Line(const RNScalar a, const RNScalar b, const RNScalar c)
    : vector(-b, a),
      normal(a, b),
      c(c)
{
}



R2Line::
R2Line(const RNScalar array[3])
    : vector(-array[1], array[0]),
      normal(array[0], array[1]),
      c(array[2])
{
}



R2Line::
R2Line(const R2Point& point, const R2Vector& vector, RNBoolean normalized)
    : vector(vector)
{
    if (!normalized) this->vector.Normalize();
    normal = R2Vector(this->vector.Y(), -(this->vector.X()));
    c = -(normal.X()*point.X() + normal.Y()*point.Y());
}



R2Line::
R2Line(const R2Point& point1, const R2Point& point2)
     : vector(point2 - point1)
{
    this->vector.Normalize();
    normal = R2Vector(this->vector.Y(), -(this->vector.X()));
    c = -(normal.X()*point1.X() + normal.Y()*point1.Y());
}



R2Line::
R2Line(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2)
    : vector(x2-x1, y2-y1)
{
    this->vector.Normalize();
    normal = R2Vector(this->vector.Y(), -(this->vector.X()));
    c = -(normal.X()*x1 + normal.Y()*y1);
}



void R2Line::
Mirror(const R2Line& line)
{
    // Mirror line over another line
    R2Point p = (normal * -c).Point();
    p.Mirror(line);
    vector.Mirror(line);
    normal = R2Vector(vector.Y(), -(vector.X()));
    c = -(normal.X()*p.X() + normal.Y()*p.Y());
}



void R2Line::
Project(const R2Line& line)
{
    // Project line onto another line
    *this = line;
    if (Vector().Dot(line.Vector()) < 0.0) Flip();
}



void R2Line::
Reposition(const R2Point& point)
{
    // Set point on line
    c = -(normal.X()*point.X() + normal.Y()*point.Y());
}



void R2Line::
Align(const R2Vector& vector, RNBoolean normalized)
{
    // Set vector of line
    this->vector = vector;
    if (!normalized) this->vector.Normalize();
    this->normal = R2Vector(this->vector.Y(), -this->vector.X());
}



void R2Line::
Transform (const R2Transformation& transformation)
{
    // Transform line
    vector.Transform(transformation);
    vector.Normalize();
    normal = R2Vector(vector.Y(), -vector.X());
    R2Point p(Point());
    p.Transform(transformation);
    Reposition(p);
}



void R2Line::
InverseTransform (const R2Transformation& transformation)
{
    // Inverse Transform line
    vector.InverseTransform(transformation);
    vector.Normalize();
    normal = R2Vector(vector.Y(), -vector.X());
    R2Point p(Point());
    p.InverseTransform(transformation);
    Reposition(p);
}



const RNBoolean R2Line::
operator==(const R2Line& line) const
{
    // Check if vectors are equal
    if (normal != line.normal) return FALSE;
    if (c != line.c) return FALSE;
    return TRUE;
}














