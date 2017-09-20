/* Include file for the R3 scene element class */



/* Class definitions */

class R3SceneElement {
public:
  // Constructor functions
  R3SceneElement(R3Material *material = NULL);
  R3SceneElement(const R3SceneElement& element);
  ~R3SceneElement(void);

  // Property functions
  const R3Box& BBox(void) const;
  const R3Point Centroid(void) const;
  const RNInterval NFacets(void) const;
  const RNLength Length(void) const;
  const RNArea Area(void) const;
  const RNVolume Volume(void) const;

  // Access functions
  R3Material *Material(void) const;
  int NShapes(void) const;
  R3Shape *Shape(int k) const;
  R3SceneNode *Node(void) const;

  // Manipulation function
  void SetMaterial(R3Material *material);
  void InsertShape(R3Shape *shape);
  void RemoveShape(R3Shape *shape);
  void Transform(const R3Transformation& transformation);

  // Query functions
  RNLength Distance(const R3Point& point) const;
  RNBoolean FindClosest(const R3Point& point, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNLength *hit_d = NULL,
    RNLength min_distance = 0, RNLength max_distance = RN_INFINITY) const;
  RNBoolean Intersects(const R3Ray& ray, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL,
    RNScalar min_t = 0.0, RNScalar max_t = RN_INFINITY) const;

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS,
    const RNArray<R3Material *> *materials = NULL) const;

public:
  // Internal update functions
  void InvalidateBBox(void);
  void UpdateBBox(void);

private:
  friend class R3SceneNode;
  R3SceneNode *node;
  R3Material *material;
  RNArray<R3Shape *> shapes;
  unsigned int opengl_id;
  R3Box bbox;
};



/* Inline functions */

inline R3Material *R3SceneElement::
Material(void) const
{
  // Return material
  return material;
}



inline int R3SceneElement::
NShapes(void) const
{
  // Return number of shapes
  return shapes.NEntries();
}



inline R3Shape *R3SceneElement::
Shape(int k) const
{
  // Return shape
  return shapes[k];
}



inline const R3Box& R3SceneElement::
BBox(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3SceneElement *) this)->UpdateBBox();
  return bbox;
}



inline R3SceneNode *R3SceneElement::
Node(void) const
{
  // Return node
  return node;
}



