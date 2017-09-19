/* Include file for the Catmull Rom spline class */



/* Class definition */

class R3Polyline : public R3Curve {
public:
  // Constructor functions
  R3Polyline(void);
  R3Polyline(const R3Polyline& polyline);
  R3Polyline(const RNArray<R3Point *>& points);
  R3Polyline(const R3Point *points, int npoints);
  R3Polyline(const R3Point *points, const RNScalar *parameters, int npoints);
  virtual ~R3Polyline(void);

  // Curve properties
  virtual const RNScalar StartParameter(void) const;
  virtual const RNScalar EndParameter(void) const;
  virtual const R3Point StartPosition(void) const;
  virtual const R3Point EndPosition(void) const;

  // Shape property functions/operators
  virtual const RNBoolean IsPoint(void) const;
  virtual const RNBoolean IsLinear(void) const;
  virtual const RNBoolean IsPlanar(void) const;
  virtual const RNLength Length(void) const;
  virtual const R3Point Centroid(void) const;
  virtual const R3Shape& BShape(void) const;
  virtual const R3Box BBox(void) const;

  // Vertex properties
  const int NVertices(void) const;
  R3Point VertexPosition(int k) const;
  RNScalar VertexParameter(int k) const;
  virtual R3Vector VertexDirection(int k) const;
  virtual R3Vector VertexDerivative(int k) const;
  virtual void *VertexData(int k) const;
  virtual int VertexIndex(RNScalar u) const;

  // Point access
  virtual R3Point PointPosition(RNScalar u) const;
  virtual R3Vector PointDerivative(RNScalar u) const;
  
  // Manipulation functions
  virtual void Transform(const R3Transformation& transformation);
  virtual void SetVertexPosition(int k, const R3Point& position);
  virtual void SetVertexParameter(int k, RNScalar u);
  virtual void SetVertexData(int k, void *data);

  // Draw functions
  virtual void Draw(const R3DrawFlags flags = R3_DEFAULT_DRAW_FLAGS) const;
  virtual void Outline(const R3DrawFlags flags = R3_DEFAULT_DRAW_FLAGS) const;

private:
  R3Point *vertex_positions;
  RNScalar *vertex_parameters;
  void **vertex_datas;
  int nvertices;

protected:
  virtual void UpdateBBox(void);
  R3Box bbox;
};



/* Inline functions */

inline const int R3Polyline::
NVertices(void) const
{
  // Return number of vertices
  return nvertices;
}



inline R3Point R3Polyline::
VertexPosition(int i) const
{
  // Return position of ith vertex
  assert((i >= 0) && (i < nvertices));
  return vertex_positions[i];
}



inline RNScalar R3Polyline::
VertexParameter(int i) const
{
  // Return parameter of ith vertex
  assert((i >= 0) && (i < nvertices));
  if (vertex_parameters) return vertex_parameters[i];
  else return i;
}



