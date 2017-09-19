// Source file for the mesh property class



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"



// Private variables

static int next_property_name = 1;



// Useful constants

const RNScalar R3_MESH_PROPERTY_KEEP_VALUE = -987654321;



////////////////////////////////////////////////////////////////////////
// Constructors/destructors
////////////////////////////////////////////////////////////////////////

R3MeshProperty::
R3MeshProperty(R3Mesh *mesh, const char *name, const RNScalar *vertex_values)
  : mesh(mesh),
    nvalues(0),
    values(NULL)
{
  // Assign property name
  if (name) strncpy(this->name, name, 1024); 
  else sprintf(this->name, "Property%d", next_property_name++);
  this->name[1023] = '\0';
  
  // Assign vertex values
  this->nvalues = mesh->NVertices();
  this->values = new RNScalar [ this->nvalues ];
  if (vertex_values) {
    for (int i = 0; i < mesh->NVertices(); i++) {
      this->values[i] = vertex_values[i];
    }
  }
  else {
    for (int i = 0; i < mesh->NVertices(); i++) {
      this->values[i] = 0;
    }
  }

  // Reset statistics
  ResetStatistics();
}



R3MeshProperty::
R3MeshProperty(const R3MeshProperty& property)
  : mesh(property.mesh),
    nvalues(0),
    values(NULL),
    mean(property.mean),
    stddev(property.stddev),
    minimum(property.minimum),
    maximum(property.maximum),
    median(property.median),
    l2norm(property.l2norm)
{
  // Copy property name
  strcpy(this->name, property.name);
  
  // Copy vertex values
  this->nvalues = mesh->NVertices();
  this->values = new RNScalar [ this->nvalues ];
  for (int i = 0; i < mesh->NVertices(); i++) {
    this->values[i] = property.values[i];
  }
}



R3MeshProperty::
~R3MeshProperty(void)
{
  // Delete vertex values
  if (values) delete [] values;
}



////////////////////////////////////////////////////////////////////////
// Vertex access/property functions
////////////////////////////////////////////////////////////////////////

RNScalar R3MeshProperty::
Gaussian(R3MeshVertex *vertex, RNScalar sigma) const
{
  // Check vertex value
  RNScalar vertex_value = VertexValue(vertex);
  if (vertex_value == RN_UNKNOWN) return RN_UNKNOWN;
  if (sigma == 0) return vertex_value;

  // Get convenient variables
  RNScalar radius = 3 * sigma;
  RNScalar denom = -2 * sigma * sigma;

  // Compute array of vertices within blur filter
  RNArray<R3MeshVertex *> neighbor_vertices;
  RNLength *distances = mesh->DijkstraDistances(vertex, radius, neighbor_vertices);
  delete [] distances;

  // Compute blurred value
  RNScalar total_value = 0;
  RNScalar total_weight = 0;
  for (int j = 0; j < neighbor_vertices.NEntries(); j++) {
    R3MeshVertex *neighbor_vertex = neighbor_vertices.Kth(j);
    int neighbor_id = mesh->VertexID(neighbor_vertex);
    RNScalar value = values[neighbor_id];
    if (value == RN_UNKNOWN) continue;
    RNLength distance = distances[neighbor_id];
    if (distance > radius) continue;
    RNScalar weight = exp(distance * distance / denom); 
    total_value += weight * value;
    total_weight += weight;
  }

  // Check total weight
  if (total_weight == 0) return vertex_value;

  // Return value 
  return total_value / total_weight;
}



RNScalar R3MeshProperty::
Laplacian(R3MeshVertex *vertex) const
{
  // Get vertex value
  R3MeshVertex *v1 = vertex;
  RNScalar value1 = VertexValue(v1);
  if (value1 == RN_UNKNOWN) return RN_UNKNOWN;
  const R3Point& p1 = mesh->VertexPosition(v1);

  // Compute vector from vertex to cotan weighted average of neighbors 
  RNScalar total_sum = 0;
  RNScalar total_weight = 0;
  for (int j = 0; j < mesh->VertexValence(v1); j++) {
    R3MeshEdge *e = mesh->EdgeOnVertex(v1, j);
    R3MeshVertex *v2 = mesh->VertexAcrossEdge(e, v1);
    RNScalar value2 = VertexValue(v2);
    if (value2 == RN_UNKNOWN) continue;
    const R3Point& p2 = mesh->VertexPosition(v2);

    // Compute cotan weight
    double weight = 0;
    for (int k = 0; k < 2; k++) {
      R3MeshFace *f = mesh->FaceOnEdge(e, k);
      if (!f) continue;
      R3MeshVertex *v3 = mesh->VertexAcrossFace(f, e);
      const R3Point& p3 = mesh->VertexPosition(v3);
      R3Vector vec1 = p1 - p3; vec1.Normalize();
      R3Vector vec2 = p2 - p3; vec2.Normalize();
      RNAngle angle = R3InteriorAngle(vec1, vec2);
      if (angle == 0) continue; 
      double tan_angle = tan(angle);
      if (tan_angle == 0) continue; 
      weight += 1.0 / tan_angle;
    }

    // Add weighted position
    total_sum += weight * value2;
    total_weight += weight;
  }

  // Check total weight
  if (total_weight == 0) return 0;

  // Compute weighted average
  RNScalar weighted_average = total_sum / total_weight;

  // Compute laplacian 
  RNScalar laplacian = value1 - weighted_average;

  // Return laplacian
  return laplacian;
}



RNBoolean R3MeshProperty::
IsLocalMinimum(R3MeshVertex *vertex) const
{
  // Return whether value at vertex is local minimum
  RNScalar value = VertexValue(mesh->VertexID(vertex));
  if (value == RN_UNKNOWN) return FALSE;
  for (int i = 0; i < mesh->VertexValence(vertex); i++) {
    R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
    R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
    RNScalar neighbor_value = VertexValue(mesh->VertexID(neighbor_vertex));
    if (neighbor_value == RN_UNKNOWN) continue;
    if (neighbor_value <= value) return FALSE;
  }

  // Vertex value is less than all neighbors
  return TRUE;
}



RNBoolean R3MeshProperty::
IsLocalMaximum(R3MeshVertex *vertex) const
{
  // Return whether value at vertex is local maximum
  RNScalar value = VertexValue(mesh->VertexID(vertex));
  if (value == RN_UNKNOWN) return FALSE;
  for (int i = 0; i < mesh->VertexValence(vertex); i++) {
    R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
    R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
    RNScalar neighbor_value = VertexValue(mesh->VertexID(neighbor_vertex));
    if (neighbor_value == RN_UNKNOWN) continue;
    if (neighbor_value >= value) return FALSE;
  }

  // Vertex value is greater than all neighbors
  return TRUE;
}



////////////////////////////////////////////////////////////////////////
// Statistics functions
////////////////////////////////////////////////////////////////////////

RNScalar R3MeshProperty::
Mean(void) const
{
  // Compute mean of property values
  if (mean == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->mean = 0;
    if (nvalues > 0) {
      // Compute sum 
      int count = 0;
      RNScalar sum = 0;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        sum += values[i];
        count++;
      }

      // Compute mean
      if (count > 0) property->mean = sum / count;
    }
  }

  // Return mean
  return mean;
}



RNScalar R3MeshProperty::
Minimum(void) const
{
  // Compute minimum of property values
  if (minimum == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->minimum = 0;
    if (nvalues > 0) {
      property->minimum = FLT_MAX;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        if (values[i] < property->minimum) property->minimum = values[i];
      }
    }
  }

  // Return minimum
  return minimum;
}



RNScalar R3MeshProperty::
Maximum(void) const
{
  // Compute maximum of property values
  if (maximum == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->maximum = 0;
    if (nvalues > 0) {
      property->maximum = -FLT_MAX;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        if (values[i] > property->maximum) property->maximum = values[i];
      }
    }
  }

  // Return maximum
  return maximum;
}



RNScalar R3MeshProperty::
StandardDeviation(void) const
{
  // Compute standard deviation of property values
  if (stddev == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->stddev = 0;
    if (nvalues > 0) {
      // Comute sum of squared residuals
      int count = 0;
      RNScalar ssd = 0;
      RNScalar avg = Mean();
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        RNScalar delta = values[i] - avg;
        ssd += delta * delta;
        count++;
      }

      // Compute standard deviation
      if (count > 0) property->stddev = sqrt(ssd / count);
    }
  }

  // Return standard deviation
  return stddev;
}



RNScalar R3MeshProperty::
Variance(void) const
{
  // Return variance of property values
  return StandardDeviation() * StandardDeviation();
}



RNScalar R3MeshProperty::
Entropy(void) const
{
  // Get useful variables
  const int nbins = 50;
  RNScalar minimum = Percentile(1);
  RNScalar maximum = Percentile(99);
  RNScalar width = maximum - minimum;
  RNScalar scale = nbins / width;

  // Check nvalues
  if (nvalues == 0) return 0;

  // Make a histogram of values
  int count = 0;
  RNScalar bins[nbins];
  for (int i = 0; i < nbins; i++) bins[i] = 0;
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    RNScalar bin = scale * (values[i] - minimum);
    int bin1 = (int) bin;
    if (bin1 < 0) bin1 = 0;
    if (bin1 >= nbins) bin1 = nbins - 1;
    int bin2 = bin1 + 1;
    if (bin2 >= nbins) bin2 = nbins - 1;
    RNScalar t = bin - bin1;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    bins[bin1] += 1 - t;
    bins[bin2] += t;
    count++;
  }

  // Normalize histogram 
  if (count > 0) {
    for (int i = 0; i < nbins; i++) {
      bins[i] /= count;
    }
  }

  // Compute entropy
  RNScalar entropy = 0;
  for (int i = 0; i < nbins; i++) {
    entropy -= bins[i] * (log(bins[i]) / log(2.0));
  }

  // Return entropy
  return entropy;
}



RNScalar R3MeshProperty::
Median(void) const
{
  // Compute mean of property values
  if (median == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->median = Percentile(50);
  }

  // Return median
  return median;
}



RNScalar R3MeshProperty::
Percentile(RNScalar percentile) const
{
  // Compute value at given percentile (100=max, 0=min)
  if (nvalues == 0) return 0;

  // Copy values
  int count = 0;
  RNScalar *copy = new RNScalar [ nvalues ];
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    copy[count++] = values[i];
  }

  // Check count of values
  if (count == 0) return 0;

  // Compute value at given percentile
  qsort(copy, count, sizeof(RNScalar), RNCompareScalars);
  int index = (int) (count * percentile / 100.0);
  if (index >= count) index = count - 1;
  RNScalar result = copy[index];

  // Delete copy
  delete [] copy;

  // Return result
  return result;
}



int R3MeshProperty::
LocalMinimumCount(void) const
{
  // Count values that are minima
  int count = 0;
  for (int i = 0; i < mesh->NVertices(); i++) 
    if (IsLocalMinimum(i)) count++;
  return count;
}



int R3MeshProperty::
LocalMaximumCount(void) const
{
  // Count values that are maxima
  int count = 0;
  for (int i = 0; i < mesh->NVertices(); i++) 
    if (IsLocalMaximum(i)) count++;
  return count;
}



int R3MeshProperty::
NonZeroCount(void) const
{
  // Count values that are non-zero
  int count = 0;
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (values[i] == 0) continue;
    count++;
  }
  return count;
}



RNScalar R3MeshProperty::
L1Norm(void) const
{
  // Compute sum 
  RNScalar sum = 0;
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    sum += fabs(values[i]);
  }

  // Return sum
  return sum;
}



RNScalar R3MeshProperty::
L2Norm(void) const
{
  // Compute L2 norm of property values
  if (l2norm == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->l2norm = 0;
    if (nvalues > 0) {
      // Compute sum 
      RNScalar sum = 0;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        sum += values[i] * values[i];
      }

      // Compute l2norm
      property->l2norm = sqrt( sum );
    }
  }

  // Return l2norm
  return l2norm;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void R3MeshProperty::
Abs(void)
{
  // Replace every value by its absolute value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = fabs(values[i]);
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Sqrt(void)
{
  // Replace every value by its sqrt
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (values[i] >= 0) values[i] = sqrt(values[i]);
    else values[i] = -1.0 * sqrt( -values[i] );
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Square(void)
{
  // Replace every value by its square
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = values[i] * values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Negate(void)
{
  // Negate every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = -values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Invert(void)
{
  // Invert every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (values[i] == 0) continue;
    values[i] = 1.0 / values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Clear(RNScalar value)
{
  // Set every value
  for (int i = 0; i < nvalues; i++) {
    values[i] = value;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Copy(const R3MeshProperty& property)
{
  // Copy every value
  for (int i = 0; i < nvalues; i++) {
    values[i] = property.values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Mask(const R3MeshProperty& property)
{
  // Zero values where property is zero
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == 0) values[i] = 0;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Substitute(RNScalar current, RNScalar replacement)
{
  // Substitute values that match
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == current) values[i] = replacement;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Add(RNScalar value)
{
  // Add scalar to every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] += value;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Add(const R3MeshProperty& property)
{
  // Add every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    values[i] += property.values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Subtract(RNScalar value)
{
  // Subtract scalar from every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] -= value;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Subtract(const R3MeshProperty& property)
{
  // Subtract every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    values[i] -= property.values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Multiply(RNScalar value)
{
  // Multiply every value by a scalar
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] *= value;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Multiply(const R3MeshProperty& property)
{
  // Multiply every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    values[i] *= property.values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Divide(RNScalar value)
{
  // Check value
  if (value == 0) return;

  // Divide every value by a scalar
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] /= value;
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Divide(const R3MeshProperty& property)
{
  // Divide every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == 0) continue;
    values[i] /= property.values[i];
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Pow(RNScalar exponent)
{
  // Raise every value to an exponent
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = pow(values[i], exponent);
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Pow(const R3MeshProperty& property)
{
  // Raise every value to an exponent
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == 0) continue;
    values[i] = pow(values[i], property.values[i]);
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
DistanceTransform(void)
{
  // Create array of source vertices
  RNArray<R3MeshVertex *> vertices;
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (values[i] == 0) continue;
    if (values[i] == RN_UNKNOWN) continue;
    R3MeshVertex *vertex = mesh->Vertex(i);
    vertices.Insert(vertex);
  }

  // Compute distances to those vertices
  RNLength *distances = mesh->DijkstraDistances(vertices);

  // Fill property
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = distances[i];
  }

  // Delete distances
  delete [] distances;

  // Reset statistics
  ResetStatistics();
}



#if 1

void R3MeshProperty::
NonExtremumSuppression(RNLength radius, RNBoolean keep_local_minima, RNBoolean keep_local_maxima)
{
  // Set every value to unknown if it has a neighbor with a more extreme value
  R3MeshProperty copy(*this);
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (copy.VertexValue(i) == RN_UNKNOWN) continue;
    R3MeshVertex *vertex = mesh->Vertex(i);
    RNScalar vertex_value = copy.VertexValue(i);

    // Find neighborhood of vertex
    RNArray<R3MeshVertex *> neighbor_vertices;
    RNLength *distances = mesh->DijkstraDistances(vertex, radius, neighbor_vertices); 
    for (int j = 0; j < mesh->VertexValence(vertex); j++) neighbor_vertices.Insert(mesh->VertexOnVertex(vertex, j));
    delete [] distances;

    // Check value of vertices within neighborhood
    RNBoolean local_minimum = TRUE;
    RNBoolean local_maximum = TRUE;
    for (int j = 0; j < neighbor_vertices.NEntries(); j++) {
      R3MeshVertex *neighbor_vertex = neighbor_vertices.Kth(j);
      RNScalar neighbor_value = copy.VertexValue(neighbor_vertex);
      if (neighbor_value == RN_UNKNOWN) continue;
      if (neighbor_value < vertex_value) local_minimum = FALSE;
      if (neighbor_value > vertex_value) local_maximum = FALSE;
    }

    // Mask this vertex if not extremum
    if (keep_local_minima && local_minimum) continue;
    if (keep_local_maxima && local_maximum) continue;
    SetVertexValue(i, RN_UNKNOWN);
  }

  // Reset statistics
  ResetStatistics();
}

#elif 0

void R3MeshProperty::
NonExtremumSuppression(RNLength radius, RNBoolean keep_local_minima, RNBoolean keep_local_maxima)
{
  // Get range of values
  RNScalar global_minimum = Minimum();
  RNScalar global_maximum = Maximum();

  // Set every value to unknown if it has a neighbor with a more extreme value
  R3MeshProperty copy(*this);
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (copy.VertexValue(i) == RN_UNKNOWN) continue;
    R3MeshVertex *vertex = mesh->Vertex(i);

    // Find neighborhood of vertex
    static RNArray<R3MeshVertex *> neighbor_vertices; neighbor_vertices.Empty();
    RNLength *distances = mesh->DijkstraDistances(vertex, radius, neighbor_vertices); 
    for (int j = 0; j < mesh->VertexValence(vertex); j++) neighbor_vertices.Insert(mesh->VertexOnVertex(vertex, j));
    delete [] distances;

    // Iteratively blur to resolve ties
    RNScalar sigma = 0;
    while (sigma < mesh->BBox().DiagonalRadius()) {
      // Get vertex value
      if (VertexValue(i) == RN_UNKNOWN) break;
      RNScalar vertex_value = copy.Gaussian(vertex, sigma);

      // Check if global extrema
      if (keep_local_minima && (vertex_value == global_minimum)) break;
      if (keep_local_maxima && (vertex_value == global_maximum)) break;

      // Check relative values of vertices within neighborhood
      int greater = 0, less = 0;
      static RNArray<R3MeshVertex *> equal; equal.Empty();
      for (int j = 0; j < neighbor_vertices.NEntries(); j++) {
        R3MeshVertex *neighbor_vertex = neighbor_vertices.Kth(j);
        if (neighbor_vertex == vertex) continue;
        RNScalar neighbor_value = copy.Gaussian(neighbor_vertex, sigma);
        if (neighbor_value == RN_UNKNOWN) continue;
        if (neighbor_value < vertex_value) less++;
        else if (neighbor_value > vertex_value) greater++;
        else equal.Insert(neighbor_vertex);
      }

      // Check if can resolve vertex
      if (greater > 0) {
        if (less > 0) {
          // On slope
          SetVertexValue(i, RN_UNKNOWN);
        }
        else {
          if (equal.IsEmpty()) {
            // Local minimum
            if (keep_local_minima) {
              break;
            }
            else {
              SetVertexValue(i, RN_UNKNOWN);
            }
          }
          else {
            // Some more, some equal
            if (keep_local_minima) {
              neighbor_vertices = equal;
            }
            else {
              SetVertexValue(i, RN_UNKNOWN);
            }
          }
        }
      }
      else {
        if (less > 0) {
          if (equal.IsEmpty()) { 
            // Local maximum 
            if (keep_local_maxima) {
              break;
            }
            else {
              SetVertexValue(i, RN_UNKNOWN);
            }
          }
          else { 
            // Some less, Some equal
            if (keep_local_maxima) {
              neighbor_vertices = equal;
            }
            else {
              SetVertexValue(i, RN_UNKNOWN);
            }
          }
        }
        else {
          if (equal.IsEmpty()) {
            // No neighbors!
            SetVertexValue(i, RN_UNKNOWN);
          }
          else {
            // All equal
            neighbor_vertices = equal;
          }
        }
      }

      // Update sigma
      if (sigma <= 0) sigma = 0.01 *  mesh->BBox().DiagonalRadius();
      else sigma *= 2;
    }
  }

  // Reset statistics
  ResetStatistics();
}

#elif 0

void R3MeshProperty::
NonExtremumSuppression(RNLength radius, RNBoolean keep_local_minima, RNBoolean keep_local_maxima)
{
  // Create temporary property
  R3MeshProperty copy(*this);

  // Compute initial blur sigma
  RNScalar sigma = 0.25 * radius;
  RNScalar max_sigma = mesh->BBox().DiagonalRadius();

  // Create array of vertices to process
  int swap_index = 0;
  RNArray<R3MeshVertex *> swap_buffer[2];
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    if (VertexValue(vertex) == RN_UNKNOWN) continue;
    swap_buffer[swap_index].Insert(vertex);
  }

  // Iteratively mask non-extrema
  while (!swap_buffer[swap_index].IsEmpty()) {
    // Empty next buffer
    swap_buffer[1-swap_index].Empty();

    // Check remaining vertices
    for (int i = 0; i < swap_buffer[swap_index].NEntries(); i++) {
      R3MeshVertex *vertex = swap_buffer[swap_index].Kth(i);
      RNScalar vertex_value = copy.VertexValue(vertex);
      assert(vertex_value != RN_UNKNOWN);

      // Find neighborhood of vertex
      RNArray<R3MeshVertex *> neighbor_vertices;
      RNLength *distances = mesh->DijkstraDistances(vertex, radius, neighbor_vertices); 
      for (int j = 0; j < mesh->VertexValence(vertex); j++) neighbor_vertices.Insert(mesh->VertexOnVertex(vertex, j));
      delete [] distances;

      // Check relationship to neighbors
      int less_value_count = 0;
      int greater_value_count = 0;
      int neighbor_count = 0;
      for (int j = 0; j < neighbor_vertices.NEntries(); j++) {
        R3MeshVertex *neighbor_vertex = neighbor_vertices.Kth(j);
        if (neighbor_vertex == vertex) continue;
        RNScalar neighbor_value = copy.VertexValue(neighbor_vertex);
        if (neighbor_value == RN_UNKNOWN) continue;
        if (neighbor_value < vertex_value) less_value_count++;
        else if (neighbor_value > vertex_value) greater_value_count++;
        neighbor_count++;
      }

      // Mask non-extrema vertices
      if (neighbor_count > 0) {
        if (less_value_count == neighbor_count) {
          // Final maximum
          if (!keep_local_maxima) SetVertexValue(vertex, RN_UNKNOWN); 
        }
        else if (greater_value_count == neighbor_count) {
          // Final minimum
          if (!keep_local_minima) SetVertexValue(vertex, RN_UNKNOWN); 
        }
        else if ((greater_value_count > 0) && (less_value_count > 0)) {
          // Not extremum
          SetVertexValue(vertex, RN_UNKNOWN); 
        }
        else if (!keep_local_maxima && (less_value_count > 0)) {
          // Not local minimum
          SetVertexValue(vertex, RN_UNKNOWN); 
        }
        else if (!keep_local_minima && (greater_value_count > 0)) {
          // Not local maximum
          SetVertexValue(vertex, RN_UNKNOWN); 
        }
        else {
          // Don't know yet -- visit again in next iteration
          swap_buffer[1-swap_index].Insert(vertex);
        }
      }
    }

    // Blur copy
    if (sigma > max_sigma) break;
    copy.Blur(sigma);
    sigma *= 2;

    // Go to next swap buffer
    swap_index = 1 - swap_index;
  } 

  // Reset statistics
  ResetStatistics();
}

#endif



void R3MeshProperty::
NonMaximumSuppression(RNScalar radius)
{
  // Keep only local maxima
  NonExtremumSuppression(radius, FALSE, TRUE);
}



void R3MeshProperty::
NonMinimumSuppression(RNScalar radius)
{
  // Keep only local minima
  NonExtremumSuppression(radius, TRUE, FALSE);
}



void R3MeshProperty::
Threshold(RNScalar threshold, RNScalar low, RNScalar high)
{
  // Threshold every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (values[i] <= threshold) {
      if (low != R3_MESH_PROPERTY_KEEP_VALUE) values[i] = low;
    }
    else {
      if (high != R3_MESH_PROPERTY_KEEP_VALUE) values[i] = high;
    }
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
AddVertexValue(int vertex_index, RNScalar value)
{
  // Add value to a vertex
  if (values[vertex_index] == RN_UNKNOWN) values[vertex_index] = value;
  else values[vertex_index] += value;

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
AddVertexValue(R3MeshVertex *vertex, RNScalar value)
{
  // Add value to a vertex
  AddVertexValue(mesh->VertexID(vertex), value);

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
SetVertexValue(int vertex_index, RNScalar value)
{
  // Set vertex value
  values[vertex_index] = value;

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
SetVertexValue(R3MeshVertex *vertex, RNScalar value)
{
  // Set vertex value
  int vertex_index = mesh->VertexID(vertex);
  SetVertexValue(vertex_index, value);

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
SetName(const char *name)
{
  // Set name
  strncpy(this->name, name, 1024);
  this->name[0123] = '\0';
}



void R3MeshProperty::
Normalize(void)
{
  // Get mean and stddev
  RNScalar m = Mean();
  RNScalar s = StandardDeviation();

  // Check stddev
  if (s == 0) {
    // Set everything to zero
    for (int i = 0; i < nvalues; i++) {
      if (values[i] == RN_UNKNOWN) continue;
      values[i] = 0;
    }
  }
  else {
    // Normalize by mean and stddev
    for (int i = 0; i < nvalues; i++) {
      if (values[i] == RN_UNKNOWN) continue;
      values[i] = (values[i] - m) / s;
    }
  }

  // Reset statistics
  ResetStatistics();
}



RNScalar
VertexPropertyValueCallback(R3MeshVertex *vertex, void *property_data)
{
  // Return property value of vertex
  R3MeshProperty *property = (R3MeshProperty *) property_data;
  return property->VertexValue(vertex);
}



void R3MeshProperty::
Percentilize(void)
{
  // Push vertices onto heap sorted by property value
  RNHeap<R3MeshVertex *> heap(VertexPropertyValueCallback, NULL, this);
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (values[i] == RN_UNKNOWN) continue;
    heap.Push(mesh->Vertex(i));
  }

  // Pop vertices in order, storing perecentiles
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (heap.IsEmpty()) break;
    R3MeshVertex *vertex = heap.Pop();
    int vertex_index = mesh->VertexID(vertex);
    values[vertex_index] = 100.0 * i / mesh->NVertices();
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Laplace(void)
{
  // Assign laplacian  at every vertex
  R3MeshProperty copy(*this);
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    if (values[i] == RN_UNKNOWN) continue;
    SetVertexValue(vertex, copy.Laplacian(vertex));
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Dilate(RNScalar radius)
{
  // Copy property
  R3MeshProperty copy(*this);
  
  // Clear property
  Clear(0);
 
  // For every vertex with non-zero value
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (copy.VertexValue(i) == 0) continue;
    if (copy.VertexValue(i) == RN_UNKNOWN) continue;
    R3MeshVertex *vertex = mesh->Vertex(i);

    // Find neighborhood of vertex
    RNArray<R3MeshVertex *> neighbor_vertices;
    RNLength *distances = mesh->DijkstraDistances(vertex, radius, neighbor_vertices);
    delete [] distances;

    // Set value of vertices within neighborhood
    for (int j = 0; j < neighbor_vertices.NEntries(); j++) {
      R3MeshVertex *neighbor_vertex = neighbor_vertices.Kth(j);
      if (VertexValue(neighbor_vertex) == RN_UNKNOWN) continue;
      SetVertexValue(neighbor_vertex, 1);
    }
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Erode(RNScalar radius)
{
  // Flip values
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (VertexValue(i) == RN_UNKNOWN) continue;
    if (VertexValue(i) == 0) SetVertexValue(i, 1);
    else SetVertexValue(i, 0);
  }

  // Dilate
  Dilate(radius);

  // Flip values
  for (int i = 0; i < mesh->NVertices(); i++) {
    if (VertexValue(i) == RN_UNKNOWN) continue;
    if (VertexValue(i) == 0) SetVertexValue(i, 1);
    else SetVertexValue(i, 0);
  }

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Blur(RNScalar sigma)
{
  // Blur across mesh

  // Get convenient variables
  if (sigma == 0) return;
  RNScalar radius = 3 * sigma;
  RNScalar denom = -2 * sigma * sigma;

  // Make copy of values
  RNScalar *old_values = new RNScalar [ mesh->NVertices() ];
  for (int i = 0; i < mesh->NVertices(); i++) old_values[i] = values[i];

  // Allocate array of vertices within blur filter
  RNArray<R3MeshVertex *> neighbor_vertices;

  // Blur value at every vertex
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    if (values[i] == RN_UNKNOWN) continue;

    // Compute distances
    RNLength *distances = mesh->DijkstraDistances(vertex, radius, neighbor_vertices);

    // Compute blurred value
    RNScalar total_value = 0;
    RNScalar total_weight = 0;
    for (int j = 0; j < neighbor_vertices.NEntries(); j++) {
      R3MeshVertex *neighbor_vertex = neighbor_vertices.Kth(j);
      int neighbor_id = mesh->VertexID(neighbor_vertex);
      RNScalar value = old_values[neighbor_id];
      if (value == RN_UNKNOWN) continue;
      RNLength distance = distances[neighbor_id];
      if (distance > radius) continue;
      RNScalar weight = exp(distance * distance / denom); 
      total_value += weight * value;
      total_weight += weight;
    }

    // Assign blurred value (normalized by total weight)
    if (total_weight > 0) {
      values[i] = total_value / total_weight;
    }

    // Delete distances
    delete [] distances;
  }

  // Delete copy of values
  delete [] old_values;

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
DoG(RNScalar sigma)
{
  // Difference of Gaussians
  // Replace values with difference from blurred version
  R3MeshProperty blur(*this); 
  blur.Blur(sigma);
  Subtract(blur);

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Strength(RNScalar sigma)
{
  // Replace values with significance score
  // For now, compute product of value and DoG
  R3MeshProperty dog(*this); 
  dog.DoG(sigma);
  dog.Abs();
  Multiply(dog);

  // Reset statistics
  ResetStatistics();
}



void R3MeshProperty::
Reset(R3Mesh *mesh)
{
  // Set mesh
  this->mesh = mesh;

  // Assign vertex values
  if (values) delete [] values;
  this->nvalues = mesh->NVertices();
  this->values = new RNScalar [ this->nvalues ];
  for (int i = 0; i < mesh->NVertices(); i++) {
    this->values[i] = 0;
  }

  // Reset statistics
  ResetStatistics();
}



////////////////////////////////////////////////////////////////////////
// Arithmetic functions
////////////////////////////////////////////////////////////////////////

RNScalar R3MeshProperty::
Dot(const R3MeshProperty& property) const
{
  // Check if number of values is consistent
  if (nvalues != property.nvalues) return 0;

  // Compute dot product of property "vectors"
  RNScalar dot = 0;
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    dot += values[i] * property.values[i];
  }

  // Return dot product
  return dot;
}



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

int R3MeshProperty::
Read(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .arff)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strcmp(extension, ".arff")) return ReadARFF(filename);
  else if (!strcmp(extension, ".txt")) return ReadValues(filename);
  else if (!strcmp(extension, ".val")) return ReadValues(filename);
  else if (!strcmp(extension, ".dat")) return ReadValues(filename);
  else if (!strcmp(extension, ".pid")) return ReadPoints(filename);

  // Reset statistics
  ResetStatistics();

  // None of the extensions matched
  fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
  return 0;
}



int R3MeshProperty::
ReadARFF(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open ARFF file: %s\n", filename);
    return 0;
  }

  // Read attribute names and set property name for first one
  int property_count = 0;
  const int max_line_size = 32 * 1024;
  char buffer[max_line_size], token[1024], name[1024];
  while (fgets(buffer, max_line_size, fp)) {
    if (strstr(buffer, "@data")) break;
    else if (strstr(buffer, "@attribute")) {
      if (sscanf(buffer, "%s%s", token, name) == (unsigned int) 2) {
        if (!strcmp(token, "@attribute")) {
          if (property_count == 0) SetName(name);
          property_count++;
        }
      }
    }
  }

  // Check number of properties
  if (property_count == 0) {
    fprintf(stderr, "No properties in %s\n", filename);
    return 0;
  }

  // Read data and assign property values
  int vertex_index = 0;
  RNScalar property_value = 0;
  while (fgets(buffer, max_line_size, fp)) {
    // Check for header line
    if (buffer[0] == '#') continue;
    if (buffer[0] == '@') continue;
    if (buffer[0] == '\0') continue;
    if (buffer[0] == '\n') continue;
    char *bufferp = strtok(buffer, "\n\t ");
    if (!bufferp) { 
      fprintf(stderr, "Unable to read property value for vertex %d\n", vertex_index);
      return 0;
    }
    if (sscanf(bufferp, "%lf", &property_value) != (unsigned int) 1) {
      fprintf(stderr, "Unable to read property value for vertex %d\n", vertex_index);
      return 0;
    }
    if (vertex_index >= mesh->NVertices()) {
      fprintf(stderr, "Too many data lines in %s\n", filename);
      return 0;
    }
    SetVertexValue(vertex_index, property_value);
    vertex_index++;
  }

  // Check if read value for every vertex
  // if (vertex_index != mesh->NVertices()) {
  //   fprintf(stderr, "Mismatching number of data lines (%d) and vertices (%d) in %s\n", vertex_index, mesh->NVertices(), filename);
  //   return 0;
  // }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3MeshProperty::
ReadValues(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open ASCII file: %s\n", filename);
    return 0;
  }

  // Create property name
  char name[1024];
  strcpy(name, filename);
  char *start = strrchr(name, '/');
  if (start) start++; 
  else start = name;
  char *end = strrchr(start, '.');
  if (end) *end = '\0';

  // Read data 
  double value;
  int vertex_index = 0;
  while (fscanf(fp, "%lf", &value) == (unsigned int) 1) {
    if (vertex_index >= mesh->NVertices()) {
      fprintf(stderr, "Too many data lines in %s\n", filename);
      return 0;
    }
    SetVertexValue(vertex_index, value);
    vertex_index++;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3MeshProperty::
ReadPoints(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open ASCII file: %s\n", filename);
    return 0;
  }

  // Create property name
  SetName("Points");

  // Initialize all values to zero
  Clear(0);

  // Assign value 1 to vertices in pid file
  int vertex_index;
  while (fscanf(fp, "%d", &vertex_index) == (unsigned int) 1) {
    if ((vertex_index > 0) && (vertex_index < mesh->NVertices())) {
      SetVertexValue(vertex_index, 1);
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3MeshProperty::
Write(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .arff)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strcmp(extension, ".arff")) return WriteARFF(filename);
  else if (!strcmp(extension, ".txt")) return WriteValues(filename);
  else if (!strcmp(extension, ".val")) return WriteValues(filename);
  else if (!strcmp(extension, ".dat")) return WriteValues(filename);

  // None of the extensions matched
  fprintf(stderr, "Unable to write file %s (unrecognized extension: %s)\n", filename, extension);
  return 0;
}



int R3MeshProperty::
WriteARFF(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open ARFF file: %s\n", filename);
    return 0;
  }

  // Write header
  const char *mesh_name = (mesh->Name()) ? mesh->Name() : "Mesh";
  fprintf(fp, "@relation %s\n\n", mesh_name);

  // Write names
  const char *property_name = (Name()) ? Name() : "Property";
  fprintf(fp, "@attribute %s real\n", property_name);

  // Write values
  fprintf(fp, "\n@data\n\n");
  for (int i = 0; i < mesh->NVertices(); i++) {
    fprintf(fp, "%g\n", VertexValue(i));
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3MeshProperty::
WriteValues(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open values file: %s\n", filename);
    return 0;
  }

  // Write values
  for (int i = 0; i < mesh->NVertices(); i++) {
    fprintf(fp, "%g\n", values[i]);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



void R3MeshProperty::
ResetStatistics(void) 
{
  // Reset all the statistics derived from values
  mean = RN_UNKNOWN;
  stddev = RN_UNKNOWN;
  minimum = RN_UNKNOWN;
  maximum = RN_UNKNOWN;
  median = RN_UNKNOWN;
  l2norm = RN_UNKNOWN;
}
