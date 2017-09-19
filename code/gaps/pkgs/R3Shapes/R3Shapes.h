/* Include file for R3 shapes module */

#ifndef __R3__SHAPES__H__
#define __R3__SHAPES__H__



/* Dependency include files */

#include "R2Shapes/R2Shapes.h"



/* Class declarations */

class R3Vector;
class R3Point;
class R3Line;
class R3Ray;
class R3Span;
class R3Plane;
class R3Halfspace;
class R4Matrix;
class R3Quaternion;
class R3Triad;
class R3CoordSystem;
class R3Transformation;
class R3Affine;
class R3Shape;
class R3Solid;
class R3Box;
class R3OrientedBox;
class R3Cylinder;
class R3Cone;
class R3Sphere;
class R3Ellipsoid;
class R3Surface;
class R3Triangle;
class R3TriangleArray;
class R3Circle;
class R3Ellipse;
class R3Rectangle;
class R3Mesh;
class R3Curve;
class R3Polyline;
class R3CatmullRomSpline;
class R3PlanarGrid;
class R3Grid;



/* Geometry basics include files */

#include "R3Shapes/R3Base.h"



/* Primitive include files */

#include "R3Shapes/R3Vector.h"
#include "R3Shapes/R3Point.h"
#include "R3Shapes/R3Line.h"
#include "R3Shapes/R3Ray.h"
#include "R3Shapes/R3Span.h"
#include "R3Shapes/R3Plane.h"
#include "R3Shapes/R3Halfspace.h"



/* Transformation include files */

#include "R3Shapes/R4Matrix.h"
#include "R3Shapes/R3Quaternion.h"
#include "R3Shapes/R3Triad.h"
#include "R3Shapes/R3Crdsys.h"
#include "R3Shapes/R3Xform.h"
#include "R3Shapes/R3Affine.h"



/* Abstract shape include files */

#include "R3Shapes/R3Shape.h"



/* Some solid shapes include files */

#include "R3Shapes/R3Solid.h"
#include "R3Shapes/R3Box.h"        



/* Surface shapes include files */

#include "R3Shapes/R3Surface.h"
#include "R3Shapes/R3Triangle.h"
#include "R3Shapes/R3TriangleArray.h"
#include "R3Shapes/R3Circle.h"
#include "R3Shapes/R3Ellipse.h"
#include "R3Shapes/R3Rectangle.h"
#include "R3Shapes/R3Mesh.h"
#include "R3Shapes/R3PlanarGrid.h"        



/* Surface shapes include files */

#include "R3Shapes/R3Curve.h"
#include "R3Shapes/R3Polyline.h"
#include "R3Shapes/R3CatmullRomSpline.h"



/* More solid shapes include files */

#include "R3Shapes/R3OrientedBox.h"        
#include "R3Shapes/R3Cylinder.h"
#include "R3Shapes/R3Cone.h"
#include "R3Shapes/R3Sphere.h"
#include "R3Shapes/R3Ellipsoid.h"
#include "R3Shapes/R3Grid.h"        



/* Shape relationship include files */

#include "R3Shapes/R3Perp.h"
#include "R3Shapes/R3Parall.h"
#include "R3Shapes/R3Dist.h"
#include "R3Shapes/R3Cont.h"
#include "R3Shapes/R3Isect.h"
#include "R3Shapes/R3Relate.h"
#include "R3Shapes/R3Align.h"
#include "R3Shapes/R3Kdtree.h"



/* Mesh utility include files */

#include "R3Shapes/R3MeshSearchTree.h"
#include "R3Shapes/R3MeshProperty.h"
#include "R3Shapes/R3MeshPropertySet.h"



/* Shape utility include files */

#include "R3Shapes/R3Draw.h"



/* Initialization functions */

int R3InitShapes(void);
void R3StopShapes(void);



#endif








