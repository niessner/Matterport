/* Include file for R2 shapes module */

#ifndef __R2__SHAPES__H__
#define __R2__SHAPES__H__



/* Dependency include files */

#include "RNBasics/RNBasics.h"



/* Class declarations */

class R2Vector;
class R2Point;
class R2Line;
class R2Ray;
class R2Span;
class R2Halfspace;
class R3Matrix;
class R2Diad;
class R2CoordSystem;
class R2Transformation;
class R2Affine;
class R2Arc;
class R2Polyline;
class R2Box;
class R2Circle;
class R2Polygon;
class R2Grid;



/* Image include files */

#include "R2Shapes/R2Image.h"



/* Primitive shape include files */

#include "R2Shapes/R2Vector.h"
#include "R2Shapes/R2Point.h"
#include "R2Shapes/R2Line.h"
#include "R2Shapes/R2Ray.h"
#include "R2Shapes/R2Span.h"
#include "R2Shapes/R2Halfspace.h"



/* Transformation include files */

#include "R2Shapes/R3Matrix.h"
#include "R2Shapes/R2Diad.h"
#include "R2Shapes/R2Crdsys.h"
#include "R2Shapes/R2Xform.h"
#include "R2Shapes/R2Affine.h"



/* Abstract shape include file */

#include "R2Shapes/R2Shape.h"



/* Solid shapes include files */

#include "R2Shapes/R2Solid.h"
#include "R2Shapes/R2Box.h"
#include "R2Shapes/R2Circle.h"
#include "R2Shapes/R2Polygon.h"



/* Curve shapes include files */

#include "R2Shapes/R2Curve.h"
#include "R2Shapes/R2Arc.h"
#include "R2Shapes/R2Polyline.h"



/* Solid shapes include files */

#include "R2Shapes/R2Grid.h"



/* Shape relationship include files */

#include "R2Shapes/R2Perp.h"
#include "R2Shapes/R2Parall.h"
#include "R2Shapes/R2Dist.h"
#include "R2Shapes/R2Cont.h"
#include "R2Shapes/R2Isect.h"
#include "R2Shapes/R2Relate.h"
#include "R2Shapes/R2Align.h"



/* Closest point search include files */

#include "R2Shapes/R2Kdtree.h"



/* Shape utility include files */

#include "R2Shapes/R2Draw.h"
#include "R2Shapes/R2Io.h"



/* Initialization functions */

int R2InitShapes(void);
void R2StopShapes(void);



#endif







