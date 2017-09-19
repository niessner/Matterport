// Source file for planar grid class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"



////////////////////////////////////////////////////////////////////////
// Member functions
////////////////////////////////////////////////////////////////////////

R3PlanarGrid::
R3PlanarGrid(void)
  : plane(0,0,0,0),
    transformation(R3identity_affine),
    grid(),
    bbox(R3null_box),
    texture_id(-1)
{
}



R3PlanarGrid::
R3PlanarGrid(const R3Plane& plane, const R3Box& world_bbox, RNLength spacing)
  : plane(0,0,0,0),
    transformation(),
    grid(),
    bbox(R3null_box),
    texture_id(-1)
{
  Reset(plane, world_bbox, spacing);
}



R3PlanarGrid::
R3PlanarGrid(const R3Plane& plane, const R3Box& world_bbox, 
  const R3Point& origin, const R3Vector& up, 
  RNLength spacing)
  : plane(plane),
    transformation(),
    grid(),
    bbox(world_bbox),
    texture_id(-1)
{
  Reset(plane, world_bbox, origin, up, spacing);
}



R3PlanarGrid::
R3PlanarGrid(const R3Plane& plane, 
  const R3Point& origin, const R3Vector& up, 
  RNLength xradius, RNLength yradius, RNLength offplane_radius,
  RNLength spacing)
  : plane(plane),
    transformation(),
    grid(),
    bbox(R3null_box),
    texture_id(-1)
{
  Reset(plane, origin, up, xradius, yradius, offplane_radius, spacing);
}



R3PlanarGrid::
R3PlanarGrid(const R3Plane& plane, 
  const R3Affine& world_to_xyplane, 
  int xres, int yres, const R2Box& xyplane_bbox)
  : plane(plane),
    transformation(world_to_xyplane),
    grid(),
    bbox(R3null_box),
    texture_id(-1)
{
  Reset(plane, world_to_xyplane, xres, yres, xyplane_bbox);
}



R3PlanarGrid::
~R3PlanarGrid(void)
{
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void R3PlanarGrid::
Reset(const R3Plane& plane, const R3Box& world_bbox, RNLength spacing)
{
  // Update plane and bbox
  this->plane = plane;
  this->bbox = world_bbox;

  // Determine transformation from 3D to 2D
  R3Point centroid = world_bbox.Centroid();
  centroid.Project(plane);
  RNDimension dim = plane.Normal().MinDimension();
  R3Vector zaxis = plane.Normal();
  R3Vector xaxis = R3xyz_triad[dim] % zaxis;
  xaxis.Normalize();
  R3Vector yaxis = zaxis % xaxis;
  yaxis.Normalize();
  R3Triad triad(xaxis, yaxis, zaxis);
  R3CoordSystem cs(centroid, triad);
  transformation.Reset(cs.InverseMatrix());

  // Determine 2D bounding box
  R2Box planar_bbox = R2null_box;
  for (int i = 0; i < 8; i++) {
    R3Point corner = world_bbox.Corner(i);
    corner.Transform(transformation);
    planar_bbox.Union(R2Point(corner[0], corner[1]));
  }

  // Determine grid resolution
  int xres = (int) (planar_bbox.XLength() / spacing + 0.5);
  int yres = (int) (planar_bbox.YLength() / spacing + 0.5);
  if (xres < 2) xres = 2;
  if (yres < 2) yres = 2;
  grid.Resample(xres, yres);

  // Determine transformation from 2D to grid
  grid.SetWorldToGridTransformation(planar_bbox);
}



void R3PlanarGrid::
Reset(const R3Plane& plane, const R3Box& world_bbox, 
  const R3Point& origin, const R3Vector& up, 
  RNLength spacing)
{
  // Set plane and bbox
  this->plane = plane;
  this->bbox = world_bbox;

  // Determine transformation from 3D to 2D
  R3Point centroid = origin;
  centroid.Project(plane);
  R3Vector zaxis = plane.Normal();
  R3Vector xaxis = up % zaxis;
  xaxis.Normalize();
  R3Vector yaxis = zaxis % xaxis;
  yaxis.Normalize();
  R3Triad triad(xaxis, yaxis, zaxis);
  R3CoordSystem cs(centroid, triad);
  transformation.Reset(cs.InverseMatrix());

  // Determine 2D bounding box
  R2Box planar_bbox = R2null_box;
  for (int i = 0; i < 8; i++) {
    R3Point corner = world_bbox.Corner(i);
    corner.Transform(transformation);
    planar_bbox.Union(R2Point(corner[0], corner[1]));
  }

  // Determine grid resolution
  int xres = (int) (planar_bbox.XLength() / spacing + 0.5);
  int yres = (int) (planar_bbox.YLength() / spacing + 0.5);
  if (xres < 2) xres = 2;
  if (yres < 2) yres = 2;
  grid.Resample(xres, yres);

  // Determine transformation from 2D to grid
  grid.SetWorldToGridTransformation(planar_bbox);
}



void R3PlanarGrid::
Reset(const R3Plane& plane, 
  const R3Point& origin, const R3Vector& up, 
  RNLength xradius, RNLength yradius, RNLength offplane_radius,
  RNLength spacing)
{
  // Set plane
  this->plane = plane;

  // Determine transformation from 3D to 2D
  R3Point centroid = origin;
  centroid.Project(plane);
  R3Vector zaxis = plane.Normal();
  R3Vector xaxis = up % zaxis;
  xaxis.Normalize();
  R3Vector yaxis = zaxis % xaxis;
  yaxis.Normalize();
  R3Triad triad(xaxis, yaxis, zaxis);
  R3CoordSystem cs(centroid, triad);
  transformation.Reset(cs.InverseMatrix());

  // Determine 3D bounding box
  R3Vector dx = xradius * cs.Axes()[RN_X];
  R3Vector dy = yradius * cs.Axes()[RN_Y];
  R3Vector dz = offplane_radius * cs.Axes()[RN_Z];
  this->bbox = R3null_box;
  bbox.Union(centroid - dx - dy - dz);
  bbox.Union(centroid - dx - dy + dz);
  bbox.Union(centroid - dx + dy - dz);
  bbox.Union(centroid - dx + dy + dz);
  bbox.Union(centroid + dx - dy - dz);
  bbox.Union(centroid + dx - dy + dz);
  bbox.Union(centroid + dx + dy - dz);
  bbox.Union(centroid + dx + dy + dz);

  // Determine 2D bounding box
  R2Box planar_bbox(-xradius, -yradius, xradius, yradius);

  // Determine grid resolution
  int xres = (int) (planar_bbox.XLength() / spacing + 0.5);
  int yres = (int) (planar_bbox.YLength() / spacing + 0.5);
  if (xres < 2) xres = 2;
  if (yres < 2) yres = 2;
  grid.Resample(xres, yres);

  // Determine transformation from 2D to grid
  grid.SetWorldToGridTransformation(planar_bbox);
}



void R3PlanarGrid::
Reset(const R3Plane& plane, 
  const R3Affine& world_to_xyplane, 
  int xres, int yres, const R2Box& xyplane_bbox)
{
  // Update stuff
  this->plane = plane;
  this->transformation = world_to_xyplane;
  this->grid.Resample(xres, yres);
  this->grid.SetWorldToGridTransformation(xyplane_bbox);
  
  // Update world bbox
  bbox = R3null_box;
  R2Box xybox = grid.WorldBox();
  R3Point p00(xybox[0][0], xybox[0][1], 0); p00.InverseTransform(world_to_xyplane); bbox.Union(p00);
  R3Point p01(xybox[0][0], xybox[1][1], 0); p01.InverseTransform(world_to_xyplane); bbox.Union(p01);
  R3Point p10(xybox[1][0], xybox[0][1], 0); p10.InverseTransform(world_to_xyplane); bbox.Union(p10);
  R3Point p11(xybox[1][0], xybox[1][1], 0); p11.InverseTransform(world_to_xyplane); bbox.Union(p11);
}



////////////////////////////////////////////////////////////////////////
// Reading/writing functions
////////////////////////////////////////////////////////////////////////

int R3PlanarGrid::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *input_extension;
  if (!(input_extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .pfm).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".grd", 4)) return ReadGridFile(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R3PlanarGrid::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *input_extension;
  if (!(input_extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Output file has no extension (e.g., .pfm).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".grd", 4)) return WriteGridFile(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// GRID FILE READ/WRITE
////////////////////////////////////////////////////////////////////////

int R3PlanarGrid::
ReadGridFile(const char *filename)
{
  // Open file
  FILE *fp = stdin;
  if (filename) {
    fp = fopen(filename, "rb");
    if (!fp) {
      RNFail("Unable to open file %s", filename);
      return 0;
    }
  }

  // Read file
  int status = ReadGrid(fp);
  if (!status) return 0;

  // Close file
  fclose(fp);

  // Return number of grid values read
  return status;
}



int R3PlanarGrid::
WriteGridFile(const char *filename) const
{
  // Open file
  FILE *fp = stdout;
  if (filename) {
    fp = fopen(filename, "wb");
    if (!fp) {
      RNFail("Unable to open file %s", filename);
      return 0;
    }
  }

  // Write file
  int status = WriteGrid(fp);
  if (!status) return 0;

  // Close file
  fclose(fp);

  // Return number of grid values written
  return status;
}



int R3PlanarGrid::
ReadGrid(FILE *fp)
{
  // Read grid 
  if (!grid.ReadGrid(fp)) return 0;

  // Read plane from file
  RNScalar p[4];
  if (fread(p, sizeof(RNScalar), 4, fp) != 4) {
    RNFail("Unable to read plane from file");
    return 0;
  }

  // Read bounding box from file
  RNScalar b[6];
  if (fread(b, sizeof(RNScalar), 6, fp) != 6) {
    RNFail("Unable to read bounding box from file");
    return 0;
  }

  // Read transformation from file
  RNScalar m[16];
  if (fread(m, sizeof(RNScalar), 16, fp) != 16) {
    RNFail("Unable to read transformation matrix from file");
    return 0;
  }

  // Update variables
  transformation.Reset(R4Matrix(m));
  plane = R3Plane(p);
  bbox.Reset(R3Point(b[0], b[1], b[2]), R3Point(b[3], b[4], b[5]));

  // Return success
  return 1;
}



int R3PlanarGrid::
WriteGrid(FILE *fp) const
{
  // Write grid 
  if (!grid.WriteGrid(fp)) return 0;

  // Write plane to file
  for (int i = 0; i < 4; i++) {
    RNScalar a = plane[i];
    if (fwrite(&a, sizeof(RNScalar), 1, fp) != 1) {
      RNFail("Unable to write plane to file");
      return 0;
    }
  }

  // Write bounding box to file
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      RNScalar a = bbox[i][j];
      if (fwrite(&a, sizeof(RNScalar), 1, fp) != 1) {
        RNFail("Unable to write plane to file");
        return 0;
      }
    }
  }

  // Write transformation to file
  const R4Matrix& matrix = transformation.Matrix();
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      RNScalar a = matrix[i][j];
      if (fwrite(&a, sizeof(RNScalar), 1, fp) != 1) {
        RNFail("Unable to write plane to file");
        return 0;
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Draw functions
////////////////////////////////////////////////////////////////////////

void R3PlanarGrid::
Draw(void) const
{
  // Just checking
  if (R3Contains(plane, R3null_plane)) return;

  // Define texture
  if (texture_id <= 0) {
    // Create texture id
    GLuint i;
    glGenTextures(1, &i);
    ((R3PlanarGrid *) this)->texture_id = i;
    assert(texture_id > 0);

    // Create texture pixels
    float *pixels = new float [ 4 * grid.NEntries() ];
    for (int i = 0; i < grid.NEntries(); i++) {
      pixels[4*i + 0] = grid.GridValue(i);
      pixels[4*i + 1] = grid.GridValue(i);
      pixels[4*i + 2] = grid.GridValue(i);
      pixels[4*i + 3] = grid.GridValue(i);
    }

    // Begin texture definition
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, grid.XResolution(), grid.YResolution(), 0, GL_RGBA, GL_FLOAT, pixels);

    // Delete texture pixels
    delete [] pixels;
  }

  // Set texture
  glBindTexture(GL_TEXTURE_2D, texture_id);

  // Push transformation
  R3Affine inverse(transformation);
  inverse.Invert();
  inverse.Push();

  // Set rendering modes
  glEnable(GL_TEXTURE_2D);
  glAlphaFunc(GL_GREATER, 0.1);
  glEnable(GL_ALPHA_TEST);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glDisable(GL_CULL_FACE);

  // Draw polygon
  R3BeginPolygon();
  const R2Box& bbox = grid.WorldBox();
  R3LoadTextureCoords(0.0, 0.0);  R2LoadPoint(bbox[0][0], bbox[0][1]);
  R3LoadTextureCoords(1.0, 0.0);  R2LoadPoint(bbox[1][0], bbox[0][1]);
  R3LoadTextureCoords(1.0, 1.0);  R2LoadPoint(bbox[1][0], bbox[1][1]);
  R3LoadTextureCoords(0.0, 1.0);  R2LoadPoint(bbox[0][0], bbox[1][1]);
  R3EndPolygon();

  // Reset rendering modes
  glEnable(GL_CULL_FACE);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_TEXTURE_2D);

  // Pop transformation
  inverse.Pop();
}



