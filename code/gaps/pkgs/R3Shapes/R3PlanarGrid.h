// Include file for planar grid class (a 2D grid on a 3D plane)



// Class definitions

class R3PlanarGrid {
public:
  // Constructor/destructor
  R3PlanarGrid(void);
  R3PlanarGrid(const R3Plane& plane, const R3Box& world_bbox, RNLength spacing);
  R3PlanarGrid(const R3Plane& plane, const R3Box& world_bbox, const R3Point& origin, const R3Vector& yaxis, RNLength spacing);
  R3PlanarGrid(const R3Plane& plane, const R3Point& origin, const R3Vector& yaxis, 
    RNLength xradius, RNLength yradius, RNLength offplane_radius, RNLength spacing);
  R3PlanarGrid(const R3Plane& plane, const R3Affine& world_to_xyplane,
    int xres, int yres, const R2Box& xybox);
  ~R3PlanarGrid(void);

  // Property functions
  int NEntries() const;
  int XResolution(void) const;
  int YResolution(void) const;
  int Resolution(RNDimension dim) const;
  RNScalar Sum(void) const;
  RNScalar Mean(void) const;
  RNScalar Median(void) const;
  RNScalar Maximum(void) const;
  RNScalar Minimum(void) const;
  RNInterval Range(void) const;
  RNScalar L1Norm(void) const;
  RNScalar L2Norm(void) const;
  RNScalar Area(void) const;
  int Cardinality(void) const;
  R2Box GridBox(void) const;
  const R3Box& WorldBox(void) const;
  const R3Plane& Plane(void) const;

  // Grid value access functions
  const R2Grid& Grid(void) const;
  RNScalar GridValue(int index) const;
  RNScalar GridValue(int i, int j) const;
  RNScalar GridValue(RNCoord x, RNCoord y) const;
  RNScalar GridValue(const R2Point& grid_point) const;
  RNScalar WorldValue(RNCoord x, RNCoord y, RNCoord z) const;
  RNScalar WorldValue(const R3Point& world_point) const;
  RNScalar& operator()(int i, int j);

  // Transformation functions
  R3Affine WorldToGridTransformation(void) const;
  R3Affine GridToWorldTransformation(void) const;
  RNScalar WorldToGridScaleFactor(void) const;
  RNScalar GridToWorldScaleFactor(void) const;

  // Grid manipulation functions
  void Abs(void);
  void Sqrt(void);
  void Square(void);
  void Negate(void);
  void Invert(void);
  void Normalize(void);
  void Laplacian(void);
  void Laplacian(int dim);
  void Sobel(void);
  void DetectEdges(void);
  void FillHoles(void);
  void GradientAngle(void);
  void GradientMagnitude(void);
  void Gradient(RNDimension dim);
  void Hessian(RNDimension dim1, RNDimension dim2);
  void Clear(RNScalar value = 0);
  void Dilate(RNScalar grid_distance);
  void Erode(RNScalar grid_distance);
  void Blur(RNScalar grid_sigma = 2);
  void BilateralFilter(RNLength grid_sigma = 2, RNScalar value_sigma = -1);
  void AnisotropicDiffusion(RNLength grid_sigma = 2, RNScalar gradient_sigma = -1);
  void PercentileFilter(RNLength grid_radius, RNScalar percentile);
  void MinFilter(RNLength grid_radius);
  void MaxFilter(RNLength grid_radius);
  void MedianFilter(RNLength grid_radius);
  void MaskNonMinima(RNLength grid_radius = 0);
  void MaskNonMaxima(RNLength grid_radius = 0);
  void Convolve(const RNScalar filter[3][3]);
  void Substitute(RNScalar old_value, RNScalar new_value);
  void Add(RNScalar value);
  void Add(const R3PlanarGrid& grid);
  void Subtract(RNScalar value);
  void Subtract(const R3PlanarGrid& grid);
  void Multiply(RNScalar value);
  void Multiply(const R3PlanarGrid& grid);
  void Divide(RNScalar value);
  void Divide(const R3PlanarGrid& grid);
  void Pow(RNScalar exponent);
  void Mask(const R3PlanarGrid& grid);
  void Threshold(RNScalar threshold, RNScalar low, RNScalar high);
  void Threshold(const R3PlanarGrid& threshold, RNScalar low, RNScalar high);
  void SignedDistanceTransform(void);
  void SquaredDistanceTransform(void);
  void Voronoi(R3PlanarGrid *squared_distance_grid = NULL);
  void PointSymmetryTransform(int radius = -1);
  void Gauss(RNLength sigma = sqrt(8.0), RNBoolean square = TRUE);
  void Resample(int xres, int yres);
  void SetGridValue(int index, RNScalar value);
  void SetGridValue(int i, int j, RNScalar value);
  void AddGridValue(int i, int j, RNScalar value);
  void Reset(const R3Plane& plane, const R3Box& world_bbox, RNLength spacing);
  void Reset(const R3Plane& plane, const R3Box& world_bbox, const R3Point& origin, const R3Vector& yaxis, RNLength spacing);
  void Reset(const R3Plane& plane, const R3Point& origin, const R3Vector& yaxis, 
    RNLength xradius, RNLength yradius, RNLength offplane_radius, RNLength spacing);
  void Reset(const R3Plane& plane, const R3Affine& world_to_xyplane,
    int xres, int yres, const R2Box& xybox);
  
  // Rasterization functions
  void RasterizeGridPoint(RNCoord x, RNCoord y, RNScalar value);
  void RasterizeWorldPoint(RNCoord x, RNCoord y, RNCoord z, RNScalar value);
  void RasterizeGridPoint(const R2Point& point, RNScalar value);
  void RasterizeWorldPoint(const R3Point& point, RNScalar value);
  void RasterizeGridSpan(const R2Point& p1, const R2Point& p2, RNScalar value);
  void RasterizeGridSpan(const R2Point& p1, const R2Point& p2, RNScalar value1, RNScalar value2);
  void RasterizeWorldSpan(const R3Point& p1, const R3Point& p2, RNScalar value);
  void RasterizeWorldSpan(const R3Point& p1, const R3Point& p2, RNScalar value1, RNScalar value2);
  void RasterizeGridTriangle(const R2Point& p1, const R2Point& p2, const R2Point& p3, RNScalar value);
  void RasterizeGridTriangle(const R2Point& p1, const R2Point& p2, const R2Point& p3, RNScalar value1, RNScalar value2, RNScalar value3);
  void RasterizeWorldTriangle(const R3Point& p1, const R3Point& p2, const R3Point& p3, RNScalar value);
  void RasterizeWorldTriangle(const R3Point& p1, const R3Point& p2, const R3Point& p3, RNScalar value1, RNScalar value2, RNScalar value3);
  void RasterizeGridCircle(const R2Point& center, RNLength radius, RNScalar value);
  void RasterizeWorldCircle(const R3Point& center, RNLength radius, RNScalar value);
  void RasterizeGridPolygon(const R2Polygon& polygon, RNScalar value);

  // Transformation utility functions
  R3Point WorldPosition(const R2Point& grid_point) const;
  R2Point GridPosition(const R3Point& world_point) const;
  R3Point WorldPosition(RNCoord x, RNCoord y) const;
  R2Point GridPosition(RNCoord x, RNCoord y, RNCoord z) const;

  // Reading/writing
  int ReadFile(const char *filename);
  int ReadGridFile(const char *filename);
  int ReadGrid(FILE *fp = NULL);
  int WriteFile(const char *filename) const;
  int WriteGridFile(const char *filename) const;
  int WriteGrid(FILE *fp = NULL) const;

  // Draw functions
  void Draw(void) const;

  // Utility functions
  void ConnectedComponentLabelFilter(RNScalar isolevel);
  void ConnectedComponentSizeFilter(RNScalar isolevel);
  void ConnectedComponentCentroidFilter(RNScalar isolevel);
  void ConnectedComponentFilter(RNScalar isolevel, RNArea min_grid_area, RNArea max_grid_area, 
    RNScalar under_isolevel_value = 0, RNScalar too_small_value = 0, RNScalar too_large_value = 0);
  int ConnectedComponents(RNScalar isolevel = 0, int max_components = 0, int *seeds = NULL, int *sizes = NULL, int *grid_components = NULL);

public:
  R3Plane plane;
  R3Affine transformation; // 3D world to 2D plane
  R2Grid grid;
  R3Box bbox;
  int texture_id;
};



// Inline functions taken directly from R2Grid

inline int R3PlanarGrid::NEntries() const 
  { return grid.NEntries(); }
inline int R3PlanarGrid::XResolution(void) const  
  { return grid.XResolution(); }
inline int R3PlanarGrid::YResolution(void) const  
  { return grid.YResolution(); }
inline int R3PlanarGrid::Resolution(RNDimension dim) const  
  { return grid.Resolution(dim); }
inline RNScalar R3PlanarGrid::Sum(void) const  
  { return grid.Sum(); }
inline RNScalar R3PlanarGrid::Mean(void) const  
  { return grid.Mean(); }
inline RNScalar R3PlanarGrid::Median(void) const  
  { return grid.Median(); }
inline RNScalar R3PlanarGrid::Maximum(void) const  
  { return grid.Maximum(); }
inline RNScalar R3PlanarGrid::Minimum(void) const 
  { return grid.Minimum(); }
inline RNInterval R3PlanarGrid::Range(void) const 
  { return grid.Range(); }
inline RNScalar R3PlanarGrid::L1Norm(void) const 
  { return grid.L1Norm(); }
inline RNScalar R3PlanarGrid::L2Norm(void) const 
  { return grid.L2Norm(); }
inline RNScalar R3PlanarGrid::Area(void) const 
  { return grid.Area(); }
inline int R3PlanarGrid::Cardinality(void) const 
  { return grid.Cardinality(); }
inline R2Box R3PlanarGrid::GridBox(void) const 
  { return grid.GridBox(); }
inline const R3Box& R3PlanarGrid::WorldBox(void) const 
  { return bbox; }
inline const R2Grid& R3PlanarGrid::Grid(void) const 
  { return grid; }
inline RNScalar R3PlanarGrid::GridValue(int index) const 
  { return grid.GridValue(index); }
inline RNScalar R3PlanarGrid::GridValue(int i, int j) const 
  { return grid.GridValue(i, j); }
inline RNScalar R3PlanarGrid::GridValue(RNCoord x, RNCoord y) const 
  { return grid.GridValue(x, y); }
inline RNScalar R3PlanarGrid::GridValue(const R2Point& grid_point) const 
  { return grid.GridValue(grid_point); }
inline RNScalar& R3PlanarGrid::operator()(int i, int j)  
  { return grid(i, j); }
inline void R3PlanarGrid::Abs(void)  
  { grid.Abs(); }
inline void R3PlanarGrid::Sqrt(void)  
  { grid.Sqrt(); }
inline void R3PlanarGrid::Square(void)  
  { grid.Square(); }
inline void R3PlanarGrid::Negate(void)  
  { grid.Negate(); }
inline void R3PlanarGrid::Invert(void)  
  { grid.Invert(); }
inline void R3PlanarGrid::Normalize(void)  
  { grid.Normalize(); }
inline void R3PlanarGrid::Laplacian(void)  
  { grid.Laplacian(); }
inline void R3PlanarGrid::Laplacian(int dim)  
  { grid.Laplacian(dim); }
inline void R3PlanarGrid::Sobel(void)  
  { grid.Sobel(); }
inline void R3PlanarGrid::DetectEdges(void)  
  { grid.DetectEdges(); }
inline void R3PlanarGrid::FillHoles(void)  
  { grid.FillHoles(); }
inline void R3PlanarGrid::GradientAngle(void)  
  { grid.GradientAngle(); }
inline void R3PlanarGrid::GradientMagnitude(void)  
  { grid.GradientMagnitude(); }
inline void R3PlanarGrid::Gradient(RNDimension dim)  
  { grid.Gradient(dim); }
inline void R3PlanarGrid::Hessian(RNDimension dim1, RNDimension dim2)  
  { grid.Hessian(dim1, dim2); }
inline void R3PlanarGrid::Clear(RNScalar value)  
  { grid.Clear(value); }
inline void R3PlanarGrid::Dilate(RNScalar grid_distance)  
  { grid.Dilate(grid_distance); }
inline void R3PlanarGrid::Erode(RNScalar grid_distance)  
  { grid.Erode(grid_distance); }
inline void R3PlanarGrid::Blur(RNScalar grid_sigma)  
  { grid.Blur(grid_sigma); }
inline void R3PlanarGrid::BilateralFilter(RNLength grid_sigma, RNScalar value_sigma)  
  { grid.BilateralFilter(grid_sigma, value_sigma); }
inline void R3PlanarGrid::AnisotropicDiffusion(RNLength grid_sigma, RNScalar gradient_sigma)  
  { grid.AnisotropicDiffusion(grid_sigma, gradient_sigma); }
inline void R3PlanarGrid::PercentileFilter(RNLength grid_radius, RNScalar percentile)  
  { grid.PercentileFilter(grid_radius, percentile); }
inline void R3PlanarGrid::MinFilter(RNLength grid_radius)  
  { grid.MinFilter(grid_radius); }
inline void R3PlanarGrid::MaxFilter(RNLength grid_radius)  
  { grid.MaxFilter(grid_radius); }
inline void R3PlanarGrid::MedianFilter(RNLength grid_radius)  
  { grid.MedianFilter(grid_radius); }
inline void R3PlanarGrid::MaskNonMinima(RNLength grid_radius)  
  { grid.MaskNonMinima(grid_radius); }
inline void R3PlanarGrid::MaskNonMaxima(RNLength grid_radius)  
  { grid.MaskNonMaxima(grid_radius); }
inline void R3PlanarGrid::Convolve(const RNScalar filter[3][3])  
  { grid.Convolve(filter); }
inline void R3PlanarGrid::Substitute(RNScalar old_value, RNScalar new_value)  
  { grid.Substitute(old_value, new_value); }
inline void R3PlanarGrid::Add(RNScalar value)  
  { grid.Add(value); }
inline void R3PlanarGrid::Add(const R3PlanarGrid& g)  
  { grid.Add(g.Grid()); }
inline void R3PlanarGrid::Subtract(RNScalar value)  
  { grid.Subtract(value); }
inline void R3PlanarGrid::Subtract(const R3PlanarGrid& g)  
  { grid.Subtract(g.Grid()); }
inline void R3PlanarGrid::Multiply(RNScalar value)  
  { grid.Multiply(value); }
inline void R3PlanarGrid::Multiply(const R3PlanarGrid& g)  
  { grid.Multiply(g.Grid()); }
inline void R3PlanarGrid::Divide(RNScalar value)  
  { grid.Divide(value); }
inline void R3PlanarGrid::Divide(const R3PlanarGrid& g)  
  { grid.Divide(g.Grid()); }
inline void R3PlanarGrid::Pow(RNScalar exponent)  
  { grid.Pow(exponent); }
inline void R3PlanarGrid::Mask(const R3PlanarGrid& g)  
  { grid.Mask(g.Grid()); }
inline void R3PlanarGrid::Threshold(RNScalar threshold, RNScalar low, RNScalar high)  
  { grid.Threshold(threshold, low, high); }
inline void R3PlanarGrid::Threshold(const R3PlanarGrid& g, RNScalar low, RNScalar high)  
  { grid.Threshold(g.Grid(), low, high); }
inline void R3PlanarGrid::SignedDistanceTransform(void)  
  { grid.SignedDistanceTransform(); }
inline void R3PlanarGrid::SquaredDistanceTransform(void)  
  { grid.SquaredDistanceTransform(); }
inline void R3PlanarGrid::Voronoi(R3PlanarGrid *sd)  
  { grid.Voronoi((sd) ? &(sd->grid) : NULL); }
inline void R3PlanarGrid::Gauss(RNLength sigma, RNBoolean square)  
  { grid.Gauss(sigma, square); }
inline void R3PlanarGrid::Resample(int xres, int yres)  
  { grid.Resample(xres, yres); }
inline void R3PlanarGrid::SetGridValue(int index, RNScalar value)  
  { grid.SetGridValue(index, value); }
inline void R3PlanarGrid::SetGridValue(int i, int j, RNScalar value)  
  { grid.SetGridValue(i, j, value); }
inline void R3PlanarGrid::AddGridValue(int i, int j, RNScalar value)  
  { grid.AddGridValue(i, j, value); }
inline void R3PlanarGrid::RasterizeGridPoint(RNCoord x, RNCoord y, RNScalar value)  
  { grid.RasterizeGridPoint(x, y, value); }
inline void R3PlanarGrid::RasterizeGridPoint(const R2Point& point, RNScalar value)  
  { grid.RasterizeGridPoint(point, value); }
inline void R3PlanarGrid::RasterizeGridSpan(const R2Point& p1, const R2Point& p2, RNScalar value)  
  { grid.RasterizeGridSpan(p1, p2, value); }
inline void R3PlanarGrid::RasterizeGridSpan(const R2Point& p1, const R2Point& p2, RNScalar value1, RNScalar value2)  
  { grid.RasterizeGridSpan(p1, p2, value1, value2); }
inline void R3PlanarGrid::RasterizeGridTriangle(const R2Point& p1, const R2Point& p2, const R2Point& p3, RNScalar value)  
  { grid.RasterizeGridTriangle(p1, p2, p3, value); }
inline void R3PlanarGrid::RasterizeGridTriangle(const R2Point& p1, const R2Point& p2, const R2Point& p3, RNScalar value1, RNScalar value2, RNScalar value3)  
  { grid.RasterizeGridTriangle(p1, p2, p3, value1, value2, value3); }
inline void R3PlanarGrid::RasterizeGridCircle(const R2Point& center, RNLength radius, RNScalar value)  
  { grid.RasterizeGridCircle(center, radius, value); }
inline void R3PlanarGrid::RasterizeGridPolygon(const R2Polygon& polygon, RNScalar value)  
  { grid.RasterizeGridPolygon(polygon, value); }
inline void R3PlanarGrid::ConnectedComponentLabelFilter(RNScalar isolevel)
  { grid.ConnectedComponentLabelFilter(isolevel); }
inline void R3PlanarGrid::ConnectedComponentSizeFilter(RNScalar isolevel)
  { grid.ConnectedComponentSizeFilter(isolevel); }
inline void R3PlanarGrid::ConnectedComponentCentroidFilter(RNScalar isolevel)
  { grid.ConnectedComponentCentroidFilter(isolevel); }
inline void R3PlanarGrid::ConnectedComponentFilter(RNScalar isolevel, 
  RNArea min_grid_area, RNArea max_grid_area,  RNScalar under_isolevel_value, RNScalar too_small_value, RNScalar too_large_value)
  { grid.ConnectedComponentFilter(isolevel, min_grid_area, max_grid_area, under_isolevel_value, too_small_value, too_large_value); }
inline int R3PlanarGrid::ConnectedComponents(RNScalar isolevel, int max_components, int *seeds, int *sizes, int *grid_components)
  { return grid.ConnectedComponents(isolevel, max_components, seeds, sizes, grid_components); }



// Inline functions for R3PlanarGrid

inline const R3Plane& R3PlanarGrid::
Plane(void) const
{ 
  // Return plane
  return plane; 
}



inline RNScalar R3PlanarGrid::
WorldValue(RNCoord x, RNCoord y, RNCoord z) const
{ 
  // Return value at position on the plane closest to (x, y, z) given in 3D world coordinates
  R3Point p(x, y, z);
  p.Transform(transformation);
  return grid.WorldValue(p.X(), p.Y()); 
}



inline RNScalar R3PlanarGrid::
WorldValue(const R3Point& world_point) const
{ 
  // Return value at position on the plane closest to point given in 3D world coordinates
  R3Point p(world_point);
  p.Transform(transformation);
  return grid.WorldValue(p.X(), p.Y()); 
}



inline R3Affine R3PlanarGrid::
WorldToGridTransformation(void) const
{
  // Create 3D transformation for plane->grid
  R3Matrix m2 = grid.WorldToGridTransformation().Matrix();
  R4Matrix m3(m2[0][0], m2[0][1], 0, m2[0][2],
              m2[1][0], m2[1][1], 0, m2[1][2],
              0,        0,        1, 0,
              0,        0,        0, 1);

  // Return transformation from 3D world to 2D grid
  R3Affine result = R3identity_affine;
  result.Transform(R3Affine(m3));
  result.Transform(transformation);
  return result;
}



inline R3Affine R3PlanarGrid::
GridToWorldTransformation(void) const
{
  // Create 3D transformation for grid->plane
  R3Matrix m2 = grid.GridToWorldTransformation().Matrix();
  R4Matrix m3(m2[0][0], m2[0][1], 0, m2[0][2],
              m2[1][0], m2[1][1], 0, m2[1][2],
              0,        0,        1, 0,
              0,        0,        0, 1);

  // Return transformation from 2D grid to 3D world
  R3Affine result = R3identity_affine;
  result.Transform(transformation.Inverse());
  result.Transform(R3Affine(m3));
  return result;
}



inline RNScalar R3PlanarGrid::
WorldToGridScaleFactor(void) const 
{
  // Return scale from world->grid
  return transformation.ScaleFactor() * grid.WorldToGridScaleFactor();
}




inline RNScalar R3PlanarGrid::
GridToWorldScaleFactor(void) const 
{
  // Return scale from grid->world
  RNScalar s = transformation.ScaleFactor();
  RNScalar plane_to_world_scale = (s != 0) ? 1.0 / s : 1.0;
  RNScalar grid_to_plane_scale = grid.GridToWorldScaleFactor();
  return grid_to_plane_scale * plane_to_world_scale;
}



inline void R3PlanarGrid::
RasterizeWorldPoint(RNCoord x, RNCoord y, RNCoord z, RNScalar value) 
{ 
  // Rasterize point at position on the plane closest to (x, y, z) given in 3D world coordinates
  R3Point p(x, y, z);
  p.Transform(transformation);
  return grid.RasterizeWorldPoint(p.X(), p.Y(), value); 
}



inline void R3PlanarGrid::
RasterizeWorldPoint(const R3Point& point, RNScalar value) 
{ 
  // Rasterize point at position on the plane closest to point given in 3D world coordinates
  R3Point p(point);
  p.Transform(transformation);
  return grid.RasterizeWorldPoint(p.X(), p.Y(), value); 
}



inline void R3PlanarGrid::
RasterizeWorldSpan(const R3Point& point1, const R3Point& point2, RNScalar value1, RNScalar value2)
{ 
  // Rasterize point at position on the plane closest to point given in 3D world coordinates
  R3Point p1(point1);
  R3Point p2(point2);
  p1.Transform(transformation);
  p2.Transform(transformation);
  return grid.RasterizeWorldSpan(R2Point(p1.X(), p1.Y()), R2Point(p2.X(), p2.Y()), value1, value2); 
}



inline void R3PlanarGrid::
RasterizeWorldSpan(const R3Point& point1, const R3Point& point2, RNScalar value)
{ 
  // Rasterize span at position on the plane closest to span given in 3D world coordinates
  return RasterizeWorldSpan(point1, point2, value, value); 
}



inline void R3PlanarGrid::
RasterizeWorldTriangle(const R3Point& point1, const R3Point& point2, const R3Point& point3, RNScalar value1, RNScalar value2, RNScalar value3)
{ 
  // Rasterize triangle at position on the plane closest to triangle given in 3D world coordinates
  R3Point p1(point1);
  R3Point p2(point2);
  R3Point p3(point3);
  p1.Transform(transformation);
  p2.Transform(transformation);
  p3.Transform(transformation);
  return grid.RasterizeWorldTriangle(R2Point(p1.X(), p1.Y()), R2Point(p2.X(), p2.Y()), R2Point(p3.X(), p3.Y()), value1, value2, value3); 
}



inline void R3PlanarGrid::
RasterizeWorldTriangle(const R3Point& point1, const R3Point& point2, const R3Point& point3, RNScalar value)
{ 
  // Rasterize triangle at position on the plane closest to triangle given in 3D world coordinates
  return RasterizeWorldTriangle(point1, point2, point3, value, value, value); 
}



inline void R3PlanarGrid::
RasterizeWorldCircle(const R3Point& center, RNLength radius, RNScalar value)
{ 
  // Rasterize circle at position on the plane closest to circle given in 3D world coordinates
  R3Point p(center);
  p.Transform(transformation);
  RNLength r = radius * transformation.ScaleFactor();
  return grid.RasterizeWorldCircle(R2Point(p.X(), p.Y()), r, value); 
}



inline R3Point R3PlanarGrid::
WorldPosition(RNCoord x, RNCoord y) const
{
  // Return world position of grid point
  R3Point p(x, y, 0);
  p.Transform(GridToWorldTransformation());
  return p;
}




inline R3Point R3PlanarGrid::
WorldPosition(const R2Point& grid_point) const
{
  // Return world position of grid point
  return WorldPosition(grid_point[0], grid_point[1]);
}



inline R2Point R3PlanarGrid::
GridPosition(const R3Point& world_point) const
{
  // Return grid position of world point
  R3Point p(world_point);
  p.Transform(WorldToGridTransformation());
  return R2Point(p[0], p[1]);
}



inline R2Point R3PlanarGrid::
GridPosition(RNCoord x, RNCoord y, RNCoord z) const
{
  // Return grid position of world point
  return GridPosition(R3Point(x, y, z));
}



