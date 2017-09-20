/* Source file for the Catmull Rom Splineclass */



/* Include files */

#include "R3Shapes/R3Shapes.h"



R3CatmullRomSpline::
R3CatmullRomSpline(void)
  : R3Polyline(),
    tao(0)
{
}



R3CatmullRomSpline::
R3CatmullRomSpline(const R3CatmullRomSpline& curve)
  : R3Polyline(curve),
    tao(curve.tao)
{
}



R3CatmullRomSpline::
R3CatmullRomSpline(const RNArray<R3Point *>& points, RNScalar tao)
  : R3Polyline(points),
    tao(tao)
{
}



R3CatmullRomSpline::
R3CatmullRomSpline(const R3Point *points, int npoints, RNScalar tao)
  : R3Polyline(points, npoints),
    tao(tao)
{
}



R3CatmullRomSpline::
R3CatmullRomSpline(const R3Point *points, const RNScalar *parameters, int npoints, RNScalar tao)
  : R3Polyline(points, parameters, npoints),
    tao(tao)
{
}



R3CatmullRomSpline::
~R3CatmullRomSpline(void)
{
}



R3Vector R3CatmullRomSpline::
VertexDerivative(int i) const
{
  // Return derivative of curve at ith vertex
  if (NVertices() < 2) return R3zero_vector;
  if (i == NVertices()-1) return Order() * (VertexPosition(NVertices()-1) - PhantomVertexPosition(NVertices()-2, 1));
  else return Order() * (PhantomVertexPosition(i, 0) - VertexPosition(i));
}



R3Point R3CatmullRomSpline::
PointPosition(RNScalar u) const
{
  // Return position at parametric value u

  // Get index variables
  int i = VertexIndex(u);
  double denom = VertexParameter(i+1) - VertexParameter(i);
  if (RNIsZero(denom)) return VertexPosition(i);
  double t = (u - VertexParameter(i)) / denom;

  // Get control vertex_positions
  R3Point P0 = (i > 0) ? VertexPosition(i-1) : VertexPosition(0);
  R3Point P1 = VertexPosition(i);
  R3Point P2 = VertexPosition(i+1);
  R3Point P3 = (i < NVertices()-2) ? VertexPosition(i+2) : VertexPosition(NVertices()-1);

  // Compute position
  R3Point p = R3zero_point;
  p += BlendingWeight(t, -1) * P0.Vector();
  p += BlendingWeight(t,  0) * P1.Vector();
  p += BlendingWeight(t,  1) * P2.Vector();
  p += BlendingWeight(t,  2) * P3.Vector();

  // Return position
  return p;
}



R3Vector R3CatmullRomSpline::
PointDerivative(RNScalar u) const
{
  // Return tangent vector at parametric value u

  // Get index variables
  int i = VertexIndex(u);
  double denom = VertexParameter(i+1) - VertexParameter(i);
  if (RNIsZero(denom)) return R3zero_vector;
  double t = (u - VertexParameter(i)) / denom;

  // Get control vertex_positions
  R3Point P0 = (i > 0) ? VertexPosition(i-1) : VertexPosition(0);
  R3Point P1 = VertexPosition(i);
  R3Point P2 = VertexPosition(i+1);
  R3Point P3 = (i < NVertices()-2) ? VertexPosition(i+2) : VertexPosition(NVertices()-1);

  // Compute derivative vector
  R3Vector v = R3zero_vector;
  v += BlendingWeightDerivative(t, -1) * P0.Vector();
  v += BlendingWeightDerivative(t,  0) * P1.Vector();
  v += BlendingWeightDerivative(t,  1) * P2.Vector();
  v += BlendingWeightDerivative(t,  2) * P3.Vector();

  // Return derivative vector
  return v;
}



R3Point R3CatmullRomSpline::
PhantomVertexPosition(int i, int j) const
{
  // Return position of phantom Bezier curve control point within ith segment
  // j=0 means first phantom point, and j=1 means second phantom point
  if ((i == 0) && (j == 0)) {
    // Return first phantom point on first curve segment
    return tao * (VertexPosition(0) + PhantomVertexPosition(0, 1));
  }
  else if ((i == NVertices()-1) && (j == 1)) {
    // Return last phantom point on last curve segment
    return tao * (VertexPosition(NVertices()-1) + PhantomVertexPosition(NVertices()-1, 0));
  }
  else {
    // Return phantom point near vertex in middle of spline
    int vertex_index = (j == 0) ? i : i+1;
    R3Point p = VertexPosition(vertex_index);
    R3Vector v = VertexPosition(vertex_index+1) - VertexPosition(vertex_index-1);
    if (j == 0) return p + tao * v;
    else return p - tao * v;
  }
}



void R3CatmullRomSpline::
UpdateBBox(void) 
{
  // Initialize bounding box
  R3Polyline::UpdateBBox();

  // Union phantom vertex_positions
  for (int i = 0; i < NVertices()-1; i++) {
    bbox.Union(PhantomVertexPosition(i, 0));
    bbox.Union(PhantomVertexPosition(i, 1));
  }
}



RNScalar R3CatmullRomSpline::
BlendingWeight(RNScalar t, int k) const
{
  // In interval i through i+1, 
  // k=-1: V(i-1)
  // k=0: V(i)
  // k=1: V(i+1)
  // k=2: V(i+2)
  // t \in [0,1]

  // Compute powers of t
  double t2 = t * t;
  double t3 = t * t2;

  // Compute blending function
  switch (k) {
  case -1: return -tao*t + 2*tao*t2 -tao*t3;
  case 0: return 1 + (tao-3)*t2 + (2-tao)*t3;
  case 1: return tao*t + (3-2*tao)*t2 + (tao-2)*t3;
  case 2: return -tao*t2 + tao*t3;
  }

  // All other vertex_positions
  return 0;
}



RNScalar R3CatmullRomSpline::
BlendingWeightDerivative(RNScalar t, int k) const
{
  // In interval i through i+1, 
  // k=-1: V(i-1)
  // k=0: V(i)
  // k=1: V(i+1)
  // k=2: V(i+2)
  // t \in [0,1]

  // Compute powers of t
  double t2 = t * t;

  // Compute blending function
  switch (k) {
  case -1: return -tao + 4*tao*t -3*tao*t2;
  case 0: return 2*(tao-3)*t + 3*(2-tao)*t2;
  case 1: return tao + 2*(3-2*tao)*t + 3*(tao-2)*t2;
  case 2: return -2*tao*t + 3*tao*t2;
  }

  // All other vertex_positions
  return 0;
}



void R3CatmullRomSpline::
Outline(const R3DrawFlags flags) const
{
  // Outline with default spacing
  Outline(0);
}



void R3CatmullRomSpline::
Outline(RNScalar sample_spacing) const
{
  // Determine sample spacing
  const int default_nsegments = 100;
  if (sample_spacing == 0) sample_spacing =  (EndParameter() - StartParameter()) / default_nsegments;
  if (sample_spacing == 0) return;
  
  // Draw spline
  glBegin(GL_LINE_STRIP);
  for (RNScalar u = StartParameter(); u <= EndParameter(); u += sample_spacing) 
    R3LoadPoint(PointPosition(u));
  glEnd();
}
