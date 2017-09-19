/* Source file for the R3 oriented box class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3OrientedBox R3null_oriented_box(R3Point(0.0, 0.0, 0.0), R3Vector(1.0, 0.0, 0.0), R3Vector(0.0, 1.0, 0.0), -1.0, -1.0, -1.0);
const R3OrientedBox R3zero_oriented_box(R3Point(0.0, 0.0, 0.0), R3Vector(1.0, 0.0, 0.0), R3Vector(0.0, 1.0, 0.0),  0.0,  0.0,  0.0);



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3OrientedBox);



/* Public functions */

int 
R3InitOrientedBox()
{
    /* Return success */
    return TRUE;
}



void 
R3StopOrientedBox()
{
}



R3OrientedBox::
R3OrientedBox(void)
{
}



R3OrientedBox::
R3OrientedBox(const R3OrientedBox& box)
    : cs(box.cs)
{
    // Copy radii from box
    this->radius[0] = box.radius[0];
    this->radius[1] = box.radius[1];
    this->radius[2] = box.radius[2];
}



R3OrientedBox::
R3OrientedBox(const R3Point& center, const R3Vector& axis0, const R3Vector& axis1, const R3Vector& axis2)
{
    // Initialize everything
    Reset(center, axis0, axis1, axis0.Length(), axis1.Length(), axis2.Length());
}



R3OrientedBox::
R3OrientedBox(const R3Point& center, 
    const R3Vector& axis0, const R3Vector& axis1, 
    RNLength radius0, RNLength radius1, RNLength radius2)
{
    // Initialize everything
    Reset(center, axis0, axis1, radius0, radius1, radius2);
}



R3OrientedBox::
R3OrientedBox(const RNArray<R3Point *>& points)
{
#if 0
    // This is correct, but produces bboxes that are too big (diagonal)
  
    // Compute coordinate system
    R3Point centroid = R3Centroid(points);
    R3Triad axes = R3PrincipleAxes(centroid, points);

    // Compute ranges
    R3Box r = R3null_box;
    for (int i = 0; i < points.NEntries(); i++) {
      R3Point *point = points.Kth(i);
      R3Vector v = *point - centroid;
      RNLength r0 = v.Dot(axes[0]);
      RNLength r1 = v.Dot(axes[1]);
      RNLength r2 = v.Dot(axes[2]);
      r.Union(R3Point(r0, r1, r2));
    }

    // Compute center
    R3Point c = r.Centroid();
    R3Point center = centroid + axes[0]*c[0] + axes[1]*c[1] + axes[2]*c[2];

    // Initialize everything
    Reset(center, axes[0], axes[1], r.XRadius(), r.YRadius(), r.ZRadius());
#else
    // Get centroid 
    R3Point centroid = R3Centroid(points);

    // Search for best axes
    R3Box best_range = R3null_box;
    R3Triad best_axes = R3xyz_triad;
    int nxsteps = 2;
    int nysteps = 2;
    int nzsteps = 4;
    for (int ix = 0; ix < nxsteps; ix++) {
      for (int iy = 0; iy < nysteps; iy++) {
        for (int iz = 0; iz < nzsteps; iz++) {
          R3Triad axes = R3xyz_triad;
          axes.Rotate(RN_Z, ix*0.5*RN_PI/nxsteps);
          axes.Rotate(RN_Z, iy*0.5*RN_PI/nysteps);
          axes.Rotate(RN_Z, iz*0.5*RN_PI/nzsteps);
    
          // Compute ranges
          R3Box range = R3null_box;
          for (int i = 0; i < points.NEntries(); i++) {
            R3Point *point = points.Kth(i);
            R3Vector v = *point - centroid;
            RNLength r0 = v.Dot(axes[0]);
            RNLength r1 = v.Dot(axes[1]);
            RNLength r2 = v.Dot(axes[2]);
            range.Union(R3Point(r0, r1, r2));
          }

          // Check if range is best
          if (best_range.IsEmpty() || (range.Volume() < best_range.Volume())) {
            best_range = range;
            best_axes = axes;
          }
        }
      }
    }
      
    // Compute center
    if (best_range.IsEmpty()) best_range = R3zero_box;
    R3Point c = best_range.Centroid();
    R3Point center = centroid + best_axes[0]*c[0] + best_axes[1]*c[1] + best_axes[2]*c[2];
           
    // Initialize everything
    Reset(center, best_axes[0], best_axes[1], best_range.XRadius(), best_range.YRadius(), best_range.ZRadius());
#endif
}



const R3Plane R3OrientedBox::
Plane(RNDirection dir, RNDimension dim) const
{
    if (IsEmpty()) return R3null_plane;
    R3Vector normal = (dir == 0) ? Axis(dim) : -1.0 * Axis(dim);
    R3Point point = cs.Origin() + Radius(dim)*normal;
    return R3Plane(point, normal);
}



const R3Plane R3OrientedBox::
Plane(RNSide side) const
{
    // Return plane for appropriate side
    switch (side) {
    case RN_LX_SIDE: return Plane(RN_LO, RN_X);
    case RN_HX_SIDE: return Plane(RN_HI, RN_X);
    case RN_LY_SIDE: return Plane(RN_LO, RN_Y);
    case RN_HY_SIDE: return Plane(RN_HI, RN_Y);
    case RN_LZ_SIDE: return Plane(RN_LO, RN_Z);
    case RN_HZ_SIDE: return Plane(RN_HI, RN_Z);
    }
    
    // Should never reach here
    return R3null_plane;
}



const R3Point R3OrientedBox::
ClosestPoint(const R3Point& point) const
{
    // Check if box is empty
    if (IsEmpty()) return R3zero_point;

    // Start with original point
    R3Point closest_point = point;

    // Iteratively project onto each plane, if outside
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Plane plane = Plane(dir, dim);
        RNScalar d = R3SignedDistance(plane, closest_point);
        if (d > 0) closest_point -= d * plane.Normal();
      }
    }

    // Return closest point
    return closest_point;
}



const R3Point R3OrientedBox::
FurthestPoint(const R3Point& point) const
{
    // Check if box is empty
    if (IsEmpty()) return R3zero_point;

    // Start with original point
    R3Point closest_point = point;

    // Iteratively project onto each plane, if inside
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Plane plane = Plane(dir, dim);
        RNScalar d = R3SignedDistance(plane, closest_point);
        if (d < 0) closest_point += d * plane.Normal();
      }
    }

    // Return closest point
    return closest_point;
}



const RNBoolean R3OrientedBox::
IsPoint(void) const
{
    // A box lies on a single point if its maximum radius is zero
    return (RNIsZero(MaxRadius()));
}



const RNBoolean R3OrientedBox::
IsLinear(void) const
{
    // Check if box is empty
    if (IsEmpty()) return FALSE;

    // A box lies on a line if its second highest radius is zero
    if (Radius(0) > Radius(1)) {
        if (Radius(1) > Radius(2)) return RNIsZero(Radius(1));
        else return RNIsZero(Radius(2));
    }
    else {
        if (Radius(0) > Radius(2)) return RNIsZero(Radius(0));
        else return RNIsZero(Radius(2));
    }
}



const RNBoolean R3OrientedBox::
IsPlanar(void) const
{
    // A box lies within a plane if it has zero min radius
    return (RNIsZero(MinRadius()));
}



const RNBoolean R3OrientedBox::
IsConvex(void) const
{
    // Check if box is empty
    if (IsEmpty()) return FALSE;

    // All boxes are convex
    return TRUE;
}



RNBoolean R3OrientedBox::
operator==(const R3OrientedBox& box) const
{
    // Return whether box is equal
    return ((Center() == Center()) &&
            (Axis(0) == box.Axis(0)) && 
	    (Axis(1) == box.Axis(1)) && 
	    (Axis(2) == box.Axis(2)) && 
	    (Radius(0) == box.Radius(0)) && 
	    (Radius(1) == box.Radius(1)) && 
	    (Radius(2) == box.Radius(2)));
}



const RNInterval R3OrientedBox::
NFacets(void) const
{
    // Check if box is empty
    if (IsEmpty()) return RNInterval(0, 0);

    // Return number of facets (polygons)
    return RNInterval(6, 6);
}



const RNArea R3OrientedBox::
Area(void) const
{
    // Check if box is empty
    if (IsEmpty()) return 0.0;

    // Return surface area of box
    RNScalar area = 0;
    area += 2*4*Radius(0)*Radius(1);
    area += 2*4*Radius(1)*Radius(2);
    area += 2*4*Radius(2)*Radius(0);
    return area;
}



const RNVolume R3OrientedBox::
Volume(void) const
{
    // Check if box is empty
    if (IsEmpty()) return 0.0;

    // Return volume of box
    return 8*Radius(0)*Radius(1)*Radius(2);
}



const R3Point R3OrientedBox::
Centroid(void) const
{
    // Check if box is empty
    if (IsEmpty()) return R3zero_point;

    // Return centroid of box
    return Center();
}



const R3Shape& R3OrientedBox::
BShape(void) const
{
    // Return this oriented box
    return *this;
}



const R3Box R3OrientedBox::
BBox(void) const
{
    // Check if box is empty
    if (IsEmpty()) return R3null_box;

    // Return bounding box 
    R3Box bbox = R3null_box;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      bbox.Union(Corner(i));
    return bbox;
}



const R3Sphere R3OrientedBox::
BSphere(void) const
{
    // Check if box is empty
    if (IsEmpty()) return R3null_sphere;

    // Return BSphere
    return R3Sphere(Center(), DiagonalRadius());
}



void R3OrientedBox::
Translate(const R3Vector& vector)
{
    // Check if box is empty
    if (IsEmpty()) return;

    // Move box coordinate system
    cs.Translate(vector);
}



void R3OrientedBox::
Reposition(const R3Point& center)
{
    // Check if box is empty
    if (IsEmpty()) return;

    // Set box center
    cs.SetOrigin(center);
}



void R3OrientedBox::
Resize(RNLength radius0, RNLength radius1, RNLength radius2) 
{
    // Check if box is empty
    if (IsEmpty()) return;

    // Set box radius
    this->radius[0] = radius0;
    this->radius[1] = radius1;
    this->radius[2] = radius2;
}



void R3OrientedBox::
Reorient(const R3Vector& axis0, const R3Vector& axis1) 
{
    // Check if box is empty
    if (IsEmpty()) return;

    // Set box orientation
    Reset(cs.Origin(), axis0, axis1, radius[0], radius[1], radius[2]);
}



void R3OrientedBox::
Transform (const R3Transformation& transformation)
{
    // Check if box is empty
    if (IsEmpty()) return;

    // Transform 
    cs.Transform(transformation);
    radius[0] *= transformation.ScaleFactor();
    radius[1] *= transformation.ScaleFactor();
    radius[2] *= transformation.ScaleFactor();
}



void R3OrientedBox::
Empty(void)
{
    // Empty box
    *this = R3null_oriented_box;
}



void R3OrientedBox::
Reset(const R3Point& center, 
    const R3Vector& axis0, const R3Vector& axis1, 
    RNLength radius0, RNLength radius1, RNLength radius2)
{
    // Check axes
    assert(RNIsPositive(axis0.Length()));
    assert(RNIsPositive(axis1.Length()));

    // Normalize and orthogonalize axes
    R3Vector axis[3];
    axis[0] = axis0;
    axis[0].Normalize();
    axis[1] = axis1;
    axis[1].Normalize();
    axis[1] = axis[0] % axis[1];
    axis[1] = axis[1] % axis[0];
    axis[1].Normalize();
    axis[2] = axis[0] % axis[1];
    axis[2].Normalize();

    // Assign radii
    radius[0] = radius0;
    radius[1] = radius1;
    radius[2] = radius2;

    // Sort axes
    R3Vector swap_axis;
    RNScalar swap_radius;
    if (radius[1] > radius[0]) {
      swap_axis = axis[0]; swap_radius = radius[0];
      axis[0] = axis[1]; radius[0] = radius[1];
      axis[1] = swap_axis; radius[1] = swap_radius;
    }
    if (radius[2] > radius[0]) {
      swap_axis = axis[0]; swap_radius = radius[0];
      axis[0] = axis[2]; radius[0] = radius[2];
      axis[2] = swap_axis; radius[2] = swap_radius;
    }
    if (radius[2] > radius[1]) {
      swap_axis = axis[1]; swap_radius = radius[1];
      axis[1] = axis[2]; radius[1] = radius[2];
      axis[2] = swap_axis; radius[2] = swap_radius;
    }

    // Ensure right-handed coordinate system
    axis[2] = axis[0] % axis[1];
    axis[2].Normalize();

    // Set coordinate system
    cs.Reset(center, R3Triad(axis[0], axis[1], axis[2]));
}


    
    










