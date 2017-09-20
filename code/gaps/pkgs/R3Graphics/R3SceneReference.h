/* Include file for the R3 scene reference class */



/* Class definitions */

class R3SceneReference {
public:
  // Constructor functions
  R3SceneReference(R3Scene *referenced_scene = NULL);
  R3SceneReference(R3Scene *referenced_scene, const RNArray<R3Material *>& materials);
  ~R3SceneReference(void);

  // Access functions
  int NMaterials(void) const;
  R3Material *Material(int k) const;
  const RNArray<R3Material *>& Materials(void) const;
  R3Scene *ReferencedScene(void) const;

  // Property functions
  const R3Box& BBox(void) const;
  const R3Point Centroid(void) const;
  const RNInterval NFacets(void) const;
  const RNLength Length(void) const;
  const RNArea Area(void) const;
  const RNVolume Volume(void) const;
  const R3Point ClosestPoint(const R3Point& point) const;
  const char *Name(void) const;
  const char *Info(const char *key) const;
  void *Data(void) const;

  // Manipulation function
  void SetReferencedScene(R3Scene *scene);
  void InsertMaterial(R3Material *material);
  void ReplaceMaterial(int k, R3Material *material);
  void InsertInfo(const char *key, const char *info);
  void ReplaceInfo(const char *key, const char *info);
  void RemoveInfo(const char *key);
  void SetName(const char *name);
  void SetData(void *data);

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS,
    const RNArray<R3Material *> *materials = NULL) const;

private:
  R3Scene *referenced_scene;
  RNArray<R3Material *> materials;
  RNSymbolTable<std::string> info;
  char *name;
  void *data;
};



/* Inline functions */

inline R3Scene *R3SceneReference::
ReferencedScene(void) const
{
  // Return referenced scene
  return referenced_scene;
}



inline int R3SceneReference::
NMaterials(void) const
{
  // Return number of materials
  return materials.NEntries();
}



inline R3Material *R3SceneReference::
Material(int k) const
{
  // Return material
  return materials.Kth(k);
}



inline const RNArray<R3Material *>& R3SceneReference::
Materials(void) const
{
  // Return materials
  return materials;
}



inline const char *R3SceneReference::
Info(const char *key) const
{
  // Return info associated with key
  std::string value;
  if (!info.Find(key, &value)) return NULL;
  return value.c_str();
}



inline const char *R3SceneReference::
Name(void) const
{
  // Return name
  return name;
}



inline void *R3SceneReference::
Data(void) const
{
  // Return user-defined data
  return data;
}



inline void R3SceneReference::
SetReferencedScene(R3Scene *scene)
{
  // Set referenced scene
  this->referenced_scene = scene;
}



inline void R3SceneReference::
InsertMaterial(R3Material *material)
{
  // Insert material
  materials.Insert(material);
}



inline void R3SceneReference::
ReplaceMaterial(int k, R3Material *material)
{
  // Replace material
  RNArrayEntry *entry = materials.KthEntry(k);
  materials.EntryContents(entry) = material;
}



inline void R3SceneReference::
SetName(const char *name)
{
  // Set name
  if (this->name) free(this->name);
  if (name) this->name = strdup(name);
  else this->name = NULL;
}



inline void R3SceneReference::
SetData(void *data)
{
  // Set data
  this->data = data;
}



