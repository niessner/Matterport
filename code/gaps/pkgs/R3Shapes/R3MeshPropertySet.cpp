// Source file for the GAPS mesh property set



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor
////////////////////////////////////////////////////////////////////////

R3MeshPropertySet::
R3MeshPropertySet(R3Mesh *mesh)
  : mesh(mesh),
    properties()
{
}



R3MeshPropertySet::
~R3MeshPropertySet(void)
{
  // Empty all properties
  Empty();
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

int R3MeshPropertySet::
PropertyIndex(R3MeshProperty *query_property) const
{
  // Search through properties for matching name
  for (int i = 0; i < NProperties(); i++) {
    R3MeshProperty *property = Property(i);
    if (property == query_property) return i;
  }

  // Not found
  return -1;
}


int R3MeshPropertySet::
PropertyIndex(const char *property_name) const
{
  // Search through properties for matching name
  for (int i = 0; i < NProperties(); i++) {
    R3MeshProperty *property = Property(i);
    if (!strcmp(property->Name(), property_name)) return i;
  }

  // Not found
  return -1;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void R3MeshPropertySet::
Empty(void)
{
  // Empty array of properties
  properties.Empty();
}



void R3MeshPropertySet::
SortByName(void)
{
  // Sort properties
  for (int i = 0; i < properties.NEntries(); i++) {
    for (int j = i+1; j < properties.NEntries(); j++) {
      if (strcmp(properties[i]->Name(), properties[j]->Name()) > 0) {
        R3MeshProperty *swap = properties[i];
        properties.EntryContents(properties.KthEntry(i)) = properties[j];
        properties.EntryContents(properties.KthEntry(j)) = swap;
      }
    }
  }
}



void R3MeshPropertySet::
Insert(R3MeshProperty *property) 
{
  // Insert reference to property 
  properties.Insert(property);
}



void R3MeshPropertySet::
Remove(R3MeshProperty *property) 
{
  // Remove reference to property 
  properties.Remove(property);
}



void R3MeshPropertySet::
Remove(int k) 
{
  // Remove kth property
  properties.RemoveKth(k);
}



void R3MeshPropertySet::
Reset(R3Mesh *mesh)
{
  // Empty array of properties
  properties.Empty();

  // Set mesh
  this->mesh = mesh;
}



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

int R3MeshPropertySet::
Read(const char *filename)
{
  // Check mesh
  if (!mesh) {
    fprintf(stderr, "Property set must be associated with mesh before file can be read: %s\n", filename);
    return 0;
  }
  
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .arff)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strcmp(extension, ".arff")) return ReadARFF(filename);
  else if (!strcmp(extension, ".trt")) return ReadToronto(filename);
  else return ReadProperty(filename);
}



int R3MeshPropertySet::
ReadARFF(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open ARFF file: %s\n", filename);
    return 0;
  }

  // Read attribute names and allocate property for each one
  int property_count = 0;
  const int max_line_size = 32 * 1024;
  char buffer[max_line_size], token[1024], name[1024];
  while (fgets(buffer, max_line_size, fp)) {
    if (strstr(buffer, "@data")) break;
    else if (strstr(buffer, "@attribute")) {
      if (sscanf(buffer, "%s%s", token, name) == (unsigned int) 2) {
        if (!strcmp(token, "@attribute")) {
          R3MeshProperty *property = new R3MeshProperty(mesh, name);
          Insert(property);
          property_count++;
        }
      }
    }
  }

  // Check number of properties
  if (property_count == 0) return 1;

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
    for (int i = 0; i < property_count; i++) {
      R3MeshProperty *property = Property(NProperties()-property_count+i);
      if (!bufferp) { 
        fprintf(stderr, "Unable to read property value %d for vertex %d\n", i, vertex_index);
        return 0;
      }
      if (sscanf(bufferp, "%lf", &property_value) != (unsigned int) 1) {
        fprintf(stderr, "Unable to read property value %d for vertex %d\n", i, vertex_index);
        return 0;
      }
      if (vertex_index >= mesh->NVertices()) {
        fprintf(stderr, "Too many data lines in %s\n", filename);
        return 0;
      }
      property->SetVertexValue(vertex_index, property_value);
      bufferp = strtok(NULL, "\n\t ");
    }
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



int R3MeshPropertySet::
ReadToronto(const char *filename)
{
  // Get base property name
  char *name = NULL;
  char buffer[1024];
  strcpy(buffer, filename);
  name = strrchr(buffer, '/');
  if (name) name++; 
  else name = buffer;
  char *end = strrchr(name, '.');
  if (end) *end = '\0';

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open toronto file: %s\n", filename);
    return 0;
  }

  // Read header
  int nprops, nfaces, dummy;
  if (fscanf(fp, "%d%d%d\n", &nprops, &nfaces, &dummy) != (unsigned int) 3) {
    fprintf(stderr, "Unable to read toronto file: %s\n", filename);
    return 0;
  }

  // Check header
  if ((nprops == 0) || (nfaces == 0) || (nfaces != mesh->NFaces())) {
    fprintf(stderr, "Invalid header in toronto file: %s\n", filename);
    return 0;
  }

  // Allocate memory for face values
  double *face_values = new double [ nprops * nfaces ];
  if (!face_values) {
    fprintf(stderr, "Unable to allocate memory for toronto file: %s\n", filename);
    return 0;
  }

  // Read face values
  for (int j = 0; j < nfaces; j++) {
    for (int i = 0; i < nprops; i++) {
      if (fscanf(fp, "%lf", &face_values[i*nfaces+j]) != (unsigned int) 1) {
        fprintf(stderr, "Unable to read value %d for face %d in %s\n", i, j, filename);
        return 0;
      }
    }
  }

  // Create properties
  for (int i = 0; i < nprops; i++) {
    // Allocate property
    char property_name[1024];
    sprintf(property_name, "%s%d", name, i);
    R3MeshProperty *property = new R3MeshProperty(mesh, property_name);

    // Compute property values by interpolation from attached faces
    for (int j = 0; j < mesh->NVertices(); j++) {
      R3MeshVertex *vertex = mesh->Vertex(j);
      RNScalar total_value = 0;
      RNScalar total_weight = 0;
      for (int k = 0; k < mesh->VertexValence(vertex); k++) {
        R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, k);
        R3MeshFace *face = mesh->FaceOnVertex(vertex, edge);
        if (!face) continue;
        int face_index = mesh->FaceID(face);
        RNScalar value = face_values[i*nfaces + face_index];
        RNScalar weight = mesh->FaceArea(face);
        total_value += weight * value;
        total_weight += weight;
      }
      if (total_weight > 0) {
        property->SetVertexValue(j, total_value / total_weight);
      }
    }

    // Insert property
    Insert(property);
  }

  // Close file
  fclose(fp);

  // Delete memory for face values
  delete [] face_values;

  // Return success
  return 1;
}



int R3MeshPropertySet::
ReadProperty(const char *filename)
{
  // Create property name
  char name[1024];
  strcpy(name, filename);
  char *start = strrchr(name, '/');
  if (start) start++; 
  else start = name;
  char *end = strrchr(start, '.');
  if (end) *end = '\0';

  // Create mesh property
  R3MeshProperty *property = new R3MeshProperty(mesh, start);
  if (!property) {
    fprintf(stderr, "Unable to create mesh property for %s\n", filename);
    return 0;
  }

  // Read data 
  if (!property->Read(filename)) {
    fprintf(stderr, "Unable to read values from %s\n", filename);
    return 0;
  }

  // Add property to set
  Insert(property);

  // Return success
  return 1;
}



int R3MeshPropertySet::
Write(const char *filename) const
{
  // Check mesh
  if (!mesh) {
    fprintf(stderr, "Property set must be associated with mesh before file can be written: %s\n", filename);
    return 0;
  }
  
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



int R3MeshPropertySet::
WriteARFF(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open ARFF file: %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "@relation %s\n\n", mesh->Name());

  // Write names
  for (int i = 0; i < NProperties(); i++) {
    R3MeshProperty *property = Property(i);
    fprintf(fp, "@attribute %s real\n", property->Name());
  }

  // Write values
  fprintf(fp, "\n@data\n\n");
  for (int i = 0; i < mesh->NVertices(); i++) {
    for (int j = 0; j < NProperties(); j++) {
      fprintf(fp, "%g ", properties[j]->VertexValue(i));
    }
    fprintf(fp, "\n");
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3MeshPropertySet::
WriteValues(const char *filename) const
{
  // Check number of properties
  for (int i = 0; i < NProperties(); i++) {
    R3MeshProperty *property = Property(i);

    // Determine file name
    char name[1024];
    strcpy(name, filename);
    if (NProperties() > 1) {
      char *extension = strrchr(name, '.');
      if (extension) { *extension = '\0'; }
      else extension = (char *) ".val";
      sprintf(name, "%s_%d.%s", name, i, extension);
    }

    // Open file
    FILE *fp = fopen(name, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open values file: %s\n", name);
      return 0;
    }

    // Write values
    for (int i = 0; i < mesh->NVertices(); i++) {
      fprintf(fp, "%g\n", property->VertexValue(i));
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



