/* Include file for the Catmull Rom spline class */



/* Class definition */

class R3CatmullRomSpline : public R3Polyline {
public:
  // Constructor functions
  R3CatmullRomSpline(void);
  R3CatmullRomSpline(const R3CatmullRomSpline& curve);
  R3CatmullRomSpline(const RNArray<R3Point *>& points, RNScalar tao = 0.5);
  R3CatmullRomSpline(const R3Point *points, int npoints, RNScalar tao = 0.5);
  R3CatmullRomSpline(const R3Point *points, const RNScalar *parameters, int npoints, RNScalar tao = 0.5);
  virtual ~R3CatmullRomSpline(void);

  // Spline properties
  int Order(void) const;

  // Vertex properties
  virtual R3Vector VertexDerivative(int k) const;

  // Point access
  virtual R3Point PointPosition(RNScalar u) const;
  virtual R3Vector PointDerivative(RNScalar u) const;

  // Draw functions
  virtual void Outline(const R3DrawFlags flags = R3_DEFAULT_DRAW_FLAGS) const;
  virtual void Outline(RNScalar sample_spacing) const;

public:
  // Utility functions
  RNScalar BlendingWeight(RNScalar t, int k) const;
  RNScalar BlendingWeightDerivative(RNScalar t, int k) const;
  R3Point PhantomVertexPosition(int i, int j) const;

private:
  // Upkeep functions
  virtual void UpdateBBox(void);

private:
  RNScalar tao;
};



/* Inline functions */

inline int R3CatmullRomSpline::
Order(void) const
{
  // Return order of polynomial
  return 3;
}



