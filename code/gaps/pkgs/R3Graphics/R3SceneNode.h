/* Include file for the R3 scene node class */



/* Initialization functions */

int R3InitSceneNode();
void R3StopSceneNode();



/* Class definition */

class R3SceneNode {
public:
  // Constructor functions
  R3SceneNode(R3Scene *scene = NULL);
  R3SceneNode(const R3SceneNode& node);
  ~R3SceneNode(void);

  // Property functions
  const R3Shape& BShape(void) const;
  const R3Box& BBox(void) const;
  const R3Point Centroid(void) const;
  const RNInterval NFacets(void) const;
  const RNLength Length(void) const;
  const RNArea Area(void) const;
  const RNVolume Volume(void) const;
  const R3Affine& Transformation(void) const;
  const char *Name(void) const;
  const char *Info(const char *key) const;
  void *Data(void) const;

  // More property functions
  RNBoolean IsAncestor(const R3SceneNode *node) const;
  RNBoolean IsDecendent(const R3SceneNode *node) const;
  R3Affine CumulativeTransformation(void) const;
  R3Affine CumulativeParentTransformation(void) const;
  R3Point ClosestPoint(const R3Point& point) const;

  // Access functions
  R3Scene *Scene(void) const;
  int SceneIndex(void) const;
  R3SceneNode *Parent(void) const;
  int ParentIndex(void) const;
  int NChildren(void) const;
  R3SceneNode *Child(int k) const;
  int NElements(void) const;
  R3SceneElement *Element(int k) const;
  int NReferences(void) const;
  R3SceneReference *Reference(int k) const;

  // Manipulation functions
  void InsertChild(R3SceneNode *node);
  void RemoveChild(R3SceneNode *node);
  void InsertElement(R3SceneElement *element);
  void RemoveElement(R3SceneElement *element);
  void InsertReference(R3SceneReference *reference);
  void RemoveReference(R3SceneReference *reference);
  void InsertInfo(const char *key, const char *info);
  void ReplaceInfo(const char *key, const char *info);
  void RemoveInfo(const char *key);
  void SetTransformation(const R3Affine& transformation);
  void Transform(const R3Affine& transformation);
  void SetName(const char *name);
  void SetData(void *data);

  // Query functions
  RNLength Distance(const R3Point& point) const;
  RNBoolean FindClosest(const R3Point& point,
    R3SceneNode **hit_node = NULL, R3Material **hit_material = NULL, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNLength *hit_d = NULL,
    RNLength min_d = 0.0, RNLength max_d = RN_INFINITY) const;
  RNBoolean Intersects(const R3Ray& ray,
    R3SceneNode **hit_node = NULL, R3Material **hit_material = NULL, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL,
    RNScalar min_t = 0.0, RNScalar max_t = RN_INFINITY) const;

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS,
    const RNArray<R3Material *> *materials = NULL) const;

public:
  // Useful bbox functions
  R3Box WorldBBox(void) const;
  R3Box LocalBBox(void) const;
  
  // Internal update functions
  void UpdateBBox(void);
  void InvalidateBBox(void);

private:
  friend class R3Scene;
  R3Scene *scene;
  int scene_index;
  R3SceneNode *parent;
  int parent_index;
  RNArray<R3SceneNode *> children;
  RNArray<R3SceneElement *> elements;
  RNArray<R3SceneReference *> references;
  RNSymbolTable<const char *> info;
  R3Affine transformation;
  R3Box bbox;
  char *name;
  void *data;
};



/* Inline functions */

inline R3Scene *R3SceneNode::
Scene(void) const
{
  // Return scene 
  return scene;
}



inline int R3SceneNode::
SceneIndex(void) const
{
  // Return index of node in scene (can be used with scene->Node(xxx))
  return scene_index;
}



inline R3SceneNode *R3SceneNode::
Parent(void) const
{
  // Return parent node
  return parent;
}



inline int R3SceneNode::
ParentIndex(void) const
{
  // Return index of node in parent (can be used with parent->Child(xxx))
  return parent_index;
}



inline int R3SceneNode::
NChildren(void) const
{
  // Return number of children nodes
  return children.NEntries();
}



inline R3SceneNode *R3SceneNode::
Child(int k) const
{
  // Return kth child node
  return children.Kth(k);
}



inline int R3SceneNode::
NElements(void) const
{
  // Return number of elements
  return elements.NEntries();
}



inline R3SceneElement *R3SceneNode::
Element(int k) const
{
  // Return kth element
  return elements.Kth(k);
}



inline int R3SceneNode::
NReferences(void) const
{
  // Return number of references
  return references.NEntries();
}



inline R3SceneReference *R3SceneNode::
Reference(int k) const
{
  // Return kth reference
  return references.Kth(k);
}



inline const R3Affine& R3SceneNode::
Transformation(void) const
{
  // Return transformation
  return transformation;
}



inline const char *R3SceneNode::
Name(void) const
{
  // Return name
  return name;
}



inline const char *R3SceneNode::
Info(const char *key) const
{
  // Return info associated with key
  const char *value = NULL;
  if (!info.Find(key, &value)) return NULL;
  return value;
}



inline void *R3SceneNode::
Data(void) const
{
  // Return user-defined data
  return data;
}



inline RNBoolean R3SceneNode::
IsAncestor(const R3SceneNode *node) const
{
  // Return whether this is an ancestor of node
  return node->IsDecendent(this);
}




