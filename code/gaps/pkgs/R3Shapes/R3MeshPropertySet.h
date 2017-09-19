// Include file for mesh property set



// Class definition

struct R3MeshPropertySet {
public:
  // Constructor/destructor
  R3MeshPropertySet(R3Mesh *mesh = NULL);
  ~R3MeshPropertySet(void);

  // Property set info
  R3Mesh *Mesh(void) const;

  // Property access functions
  int NProperties(void) const;
  R3MeshProperty *Property(int k) const;
  R3MeshProperty *Property(const char *name) const;
  int PropertyIndex(R3MeshProperty *property) const;
  int PropertyIndex(const char *property_name) const;
  R3MeshProperty *operator[](int k) const;

  // Property manipulation functions
  void Empty(void);
  void SortByName(void);
  void Insert(R3MeshProperty *property);
  void Insert(R3MeshPropertySet *set);
  void Remove(R3MeshProperty *property);
  void Remove(R3MeshPropertySet *set);
  void Remove(const char *name);
  void Remove(int k);
  void Reset(R3Mesh *mesh);
  
  // Input/output functions
  int Read(const char *filename);
  int ReadARFF(const char *filename);
  int ReadToronto(const char *filename);
  int ReadProperty(const char *filename);
  int Write(const char *filename) const;
  int WriteARFF(const char *filename) const;
  int WriteValues(const char *filename) const;

public:
  R3Mesh *mesh;
  RNArray<R3MeshProperty *> properties;
};



// Inline functions

inline R3Mesh *R3MeshPropertySet::
Mesh(void) const
{
  // Return mesh
  return mesh;
}



inline int R3MeshPropertySet::
NProperties(void) const
{
  // Return the number of properties
  return properties.NEntries();
}



inline R3MeshProperty *R3MeshPropertySet::
Property(int k) const
{
  // Return Kth property 
  return properties.Kth(k);
}



inline R3MeshProperty *R3MeshPropertySet::
Property(const char *property_name) const
{
  // Return Kth property 
  int index = PropertyIndex(property_name);
  if (index < 0) return NULL;
  else return Property(index);
}



inline R3MeshProperty *R3MeshPropertySet::
operator[](int k) const
{
  // Return kth property
  return Property(k);
}



inline void R3MeshPropertySet::
Insert(R3MeshPropertySet *set) 
{
  // Insert all properties in set
  for (int i = 0; i < set->NProperties(); i++) {
    Insert(set->Property(i));
  }
}


inline void R3MeshPropertySet::
Remove(R3MeshPropertySet *set) 
{
  // Remove all properties in set
  for (int i = 0; i < set->NProperties(); i++) {
    Remove(set->Property(i));
  }
}




