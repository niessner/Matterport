// Include file for mesh property class



// Class definition

struct R3MeshProperty {
public:
  // Constructors/destructors
  R3MeshProperty(R3Mesh *mesh, const char *name = NULL, const RNScalar *vertex_values = NULL);
  R3MeshProperty(const R3MeshProperty& property);
  ~R3MeshProperty(void);

  // Property info functions
  R3Mesh *Mesh(void) const;
  const char *Name(void) const;

  // Value access functions
  int NVertexValues(void) const;
  RNScalar VertexValue(R3MeshVertex *vertex) const;
  RNScalar VertexValue(int vertex_index) const;

  // Vertex property functions
  RNScalar Gaussian(R3MeshVertex *vertex, RNLength sigma) const;
  RNScalar Laplacian(R3MeshVertex *vertex) const;
  RNBoolean IsLocalMaximum(R3MeshVertex *vertex) const;
  RNBoolean IsLocalMaximum(int vertex_index) const;
  RNBoolean IsLocalMinimum(R3MeshVertex *vertex) const;
  RNBoolean IsLocalMinimum(int vertex_index) const;
  RNBoolean IsLocalExtremum(R3MeshVertex *vertex) const;
  RNBoolean IsLocalExtremum(int vertex_index) const;

  // Statistics functions
  RNScalar Mean(void) const;
  RNScalar Variance(void) const;
  RNScalar StandardDeviation(void) const;
  RNScalar Minimum(void) const;
  RNScalar Maximum(void) const;
  RNScalar Entropy(void) const;
  RNScalar Median(void) const;
  RNScalar Percentile(RNScalar percentile) const;
  int LocalMaximumCount(void) const;
  int LocalMinimumCount(void) const;
  int NonZeroCount(void) const;
  RNScalar L1Norm(void) const;
  RNScalar L2Norm(void) const;

  // Manipulation functions
  void Abs(void);
  void Sqrt(void);
  void Square(void);
  void Negate(void);
  void Invert(void);
  void Clear(RNScalar value = 0);
  void Copy(const R3MeshProperty& property);
  void Mask(const R3MeshProperty& property);
  void Substitute(RNScalar current, RNScalar replacement);
  void Add(RNScalar value);
  void Add(const R3MeshProperty& property);
  void Subtract(RNScalar value);
  void Subtract(const R3MeshProperty& property);
  void Multiply(RNScalar value);
  void Multiply(const R3MeshProperty& property);
  void Divide(RNScalar value);
  void Divide(const R3MeshProperty& property);
  void Pow(RNScalar exponent);
  void Pow(const R3MeshProperty& property);
  void DistanceTransform(void);
  void NonExtremumSuppression(RNScalar radius = 0, RNBoolean keep_local_minima = TRUE, RNBoolean keep_local_maxima = TRUE);
  void NonMaximumSuppression(RNScalar radius = 0);
  void NonMinimumSuppression(RNScalar radius = 0);
  void Threshold(RNScalar threshold, RNScalar low, RNScalar high);
  void AddVertexValue(R3MeshVertex *vertex, RNScalar value);
  void AddVertexValue(int vertex_index, RNScalar value);
  void SetVertexValue(R3MeshVertex *vertex, RNScalar value);
  void SetVertexValue(int vertex_index, RNScalar value);
  void SetName(const char *name);
  void Reset(R3Mesh *mesh);

  // More property manipulation functions
  void Normalize(void);
  void Percentilize(void);
  void Laplace(void);
  void Dilate(RNScalar radius);
  void Erode(RNScalar radius);
  void Blur(RNScalar sigma);
  void DoG(RNScalar sigma);
  void Strength(RNScalar sigma);
  void Prominence(RNScalar sigma);

  // Arithmetic functions
  RNScalar Dot(const R3MeshProperty& property) const;

  // Arithmetic operators
  R3MeshProperty& operator=(const R3MeshProperty& property);
  R3MeshProperty& operator+=(RNScalar scale);
  R3MeshProperty& operator+=(const R3MeshProperty& property);
  R3MeshProperty& operator-=(RNScalar scale);
  R3MeshProperty& operator-=(const R3MeshProperty& property);
  R3MeshProperty& operator*=(RNScalar scale);
  R3MeshProperty& operator*=(const R3MeshProperty& property);
  R3MeshProperty& operator/=(RNScalar scale);
  R3MeshProperty& operator/=(const R3MeshProperty& property);

  // Input/output functions
  int Read(const char *filename);
  int ReadValues(const char *filename);
  int ReadARFF(const char *filename);
  int ReadPoints(const char *filename);
  int Write(const char *filename) const;
  int WriteValues(const char *filename) const;
  int WriteARFF(const char *filename) const;

public:
  // Internal functions
  void ResetStatistics(void);

public:
  R3Mesh *mesh;
  char name[1024];
  int nvalues;
  RNScalar *values;
  RNScalar mean;
  RNScalar stddev;
  RNScalar minimum;
  RNScalar maximum;
  RNScalar median;
  RNScalar l2norm;
};



// Useful constants

extern const RNScalar R3_MESH_PROPERTY_KEEP_VALUE;



// Inline functions

inline R3Mesh *R3MeshProperty::
Mesh(void) const
{
  // Return mesh
  return mesh;
}



inline const char *R3MeshProperty::
Name(void) const
{
  // Return property name
  return name;
}



inline int R3MeshProperty::
NVertexValues(void) const
{
  // Return number of vertex values
  return nvalues;
}




inline RNScalar R3MeshProperty::
VertexValue(int vertex_index) const
{
  // Return value at vertex
  return values[vertex_index];
}




inline RNScalar R3MeshProperty::
VertexValue(R3MeshVertex *vertex) const
{
  // Return value at vertex
  int vertex_index = mesh->VertexID(vertex);
  return VertexValue(vertex_index);
}




inline RNBoolean R3MeshProperty::
IsLocalMinimum(int vertex_index) const
{
  // Return whether vertex value is local minimum
  R3MeshVertex *vertex = mesh->Vertex(vertex_index);
  return IsLocalMinimum(vertex);
}




inline RNBoolean R3MeshProperty::
IsLocalMaximum(int vertex_index) const
{
  // Return whether vertex value is local maximum
  R3MeshVertex *vertex = mesh->Vertex(vertex_index);
  return IsLocalMaximum(vertex);
}




inline RNBoolean R3MeshProperty::
IsLocalExtremum(R3MeshVertex *vertex) const
{
  // Return whether vertex value is local extremum
  return (IsLocalMinimum(vertex) || IsLocalMaximum(vertex));
}




inline RNBoolean R3MeshProperty::
IsLocalExtremum(int vertex_index) const
{
  // Return whether vertex value is local minimum
  R3MeshVertex *vertex = mesh->Vertex(vertex_index);
  return IsLocalExtremum(vertex);
}




inline R3MeshProperty& R3MeshProperty::
operator+=(RNScalar value) 
{
  // Add value to all property values 
  Add(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator+=(const R3MeshProperty& property) 
{
  // Add passed property values to corresponding entries of this property
  Add(property);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator-=(RNScalar value) 
{
  // Subtract value from all property values 
  Subtract(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator-=(const R3MeshProperty& property) 
{
  // Subtract passed property values from corresponding entries of this property
  Subtract(property);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator*=(RNScalar value) 
{
  // Multiply property values by value
  Multiply(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator*=(const R3MeshProperty& property) 
{
  // Multiply passed property values by corresponding entries of this property
  Multiply(property);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator/=(RNScalar value) 
{
  // Divide property values by value
  Divide(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator/=(const R3MeshProperty& property) 
{
  // Divide passed property values by corresponding entries of this property
  Divide(property);
  return *this;
}



