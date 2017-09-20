// Include file for R3 view frustum class



////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

struct R3Frustum {
public:
  // Constructor
  R3Frustum(void);
  R3Frustum(const R3Camera& camera);
  R3Frustum(const R3Point& viewpoint, const R3Vector& towards, 
    const R3Vector& up, RNAngle xfov, RNAngle yfov, 
    RNLength neardist, RNLength fardist);

  // Property functions
  const R3Camera& Camera(void) const;
  const R3Halfspace& Halfspace(int dir, int dim) const;

  // Manipulation functions
  void SetCamera(const R3Camera& camera);

  // Intersection/containment functions
  RNBoolean Intersects(const R3Point& point) const;
  RNBoolean Intersects(const R3Span& span) const;
  RNBoolean Intersects(const R3Box& box) const;
  RNBoolean Contains(const R3Point& point) const;
  RNBoolean Contains(const R3Span& span) const;
  RNBoolean Contains(const R3Box& box) const;

  // Draw function
  void Draw(void) const;

public:
  R3Camera camera;
  R3Halfspace halfspaces[2][3];
};



// Inline functions

inline const R3Camera& R3Frustum::
Camera(void) const
{
  // Return camera
  return camera;
}



inline const R3Halfspace& R3Frustum::
Halfspace(int dir, int dim) const
{
  // Return halfspace
  return halfspaces[dir][dim];
}

