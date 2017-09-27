////////////////////////////////////////////////////////////////////////
// Source file for mp structures
////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include "RGBD/RGBD.h"
#include "mp.h"



////////////////////////////////////////////////////////////////////////
// UTILITY DRAWING FUNCTIONS
////////////////////////////////////////////////////////////////////////

static void 
DrawText(const R3Point& p, const char *s)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}



#if 0
static void 
DrawText(const R3Point& p, RNScalar value)
{
  // Draw text string s and position p
  char buffer[4096];
  sprintf(buffer, "%g", value);
  DrawText(p, buffer);
}
#endif



////////////////////////////////////////////////////////////////////////
// UTILITY COLORING FUNCTIONS
////////////////////////////////////////////////////////////////////////

static void
LoadIndex(int index, int tag)
{
  // Set color to represent an integer (24 bits)
  int k = index + 1;
  unsigned char color[4];
  color[0] = (k >> 16) & 0xFF;
  color[1] = (k >>  8) & 0xFF;
  color[2] = (k      ) & 0xFF;
  color[3] = tag;
  glColor4ubv(color);
}



static void
LoadColor(int k)
{
  // Make array of colors
  const int ncolors = 72;
  const RNRgb colors[ncolors] = {
    RNRgb(0.5, 0.5, 0.5), RNRgb(1, 0, 0), RNRgb(0, 0, 1), 
    RNRgb(0, 1, 0), RNRgb(0, 1, 1), RNRgb(1, 0, 1), 
    RNRgb(1, 0.5, 0), RNRgb(0, 1, 0.5), RNRgb(0.5, 0, 1), 
    RNRgb(0.5, 1, 0), RNRgb(0, 0.5, 1), RNRgb(1, 0, 0.5), 
    RNRgb(0.5, 0, 0), RNRgb(0, 0.5, 0), RNRgb(0, 0, 0.5), 
    RNRgb(0.5, 0.5, 0), RNRgb(0, 0.5, 0.5), RNRgb(0.5, 0, 0.5),
    RNRgb(0.7, 0, 0), RNRgb(0, 0.7, 0), RNRgb(0, 0, 0.7), 
    RNRgb(0.7, 0.7, 0), RNRgb(0, 0.7, 0.7), RNRgb(0.7, 0, 0.7), 
    RNRgb(0.7, 0.3, 0), RNRgb(0, 0.7, 0.3), RNRgb(0.3, 0, 0.7), 
    RNRgb(0.3, 0.7, 0), RNRgb(0, 0.3, 0.7), RNRgb(0.7, 0, 0.3), 
    RNRgb(0.3, 0, 0), RNRgb(0, 0.3, 0), RNRgb(0, 0, 0.3), 
    RNRgb(0.3, 0.3, 0), RNRgb(0, 0.3, 0.3), RNRgb(0.3, 0, 0.3),
    RNRgb(1, 0.3, 0.3), RNRgb(0.3, 1, 0.3), RNRgb(0.3, 0.3, 1), 
    RNRgb(1, 1, 0.3), RNRgb(0.3, 1, 1), RNRgb(1, 0.3, 1), 
    RNRgb(1, 0.5, 0.3), RNRgb(0.3, 1, 0.5), RNRgb(0.5, 0.3, 1), 
    RNRgb(0.5, 1, 0.3), RNRgb(0.3, 0.5, 1), RNRgb(1, 0.3, 0.5), 
    RNRgb(0.5, 0.3, 0.3), RNRgb(0.3, 0.5, 0.3), RNRgb(0.3, 0.3, 0.5), 
    RNRgb(0.5, 0.5, 0.3), RNRgb(0.3, 0.5, 0.5), RNRgb(0.5, 0.3, 0.5),
    RNRgb(0.3, 0.5, 0.5), RNRgb(0.5, 0.3, 0.5), RNRgb(0.5, 0.5, 0.3), 
    RNRgb(0.3, 0.3, 0.5), RNRgb(0.5, 0.3, 0.3), RNRgb(0.3, 0.5, 0.3), 
    RNRgb(0.3, 0.8, 0.5), RNRgb(0.5, 0.3, 0.8), RNRgb(0.8, 0.5, 0.3), 
    RNRgb(0.8, 0.3, 0.5), RNRgb(0.5, 0.8, 0.3), RNRgb(0.3, 0.5, 0.8), 
    RNRgb(0.8, 0.5, 0.5), RNRgb(0.5, 0.8, 0.5), RNRgb(0.5, 0.5, 0.8), 
    RNRgb(0.8, 0.8, 0.5), RNRgb(0.5, 0.8, 0.8), RNRgb(0.8, 0.5, 0.8)
  };

  // Load color
  if (k == -1) glColor3d(0.8, 0.8, 0.8);
  else if (k == 0) RNLoadRgb(colors[0]);
  else RNLoadRgb(colors[1 + (k % (ncolors-1))]);
}



static void
LoadColor(const char *s)
{
  // Load color based on first character of string
  LoadColor((int) ((s) ? s[0] : 0));
}


    
////////////////////////////////////////////////////////////////////////
// IMAGE MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPImage::
MPImage(void)
  : house(NULL),
    house_index(-1),
    panorama(NULL),
    panorama_index(-1),
    name(NULL),
    camera_index(-1),
    yaw_index(-1),
    rgbd(),
    extrinsics(1,0,0,0,0,1,0,0, 0,0,1,0,0,0,0,1),
    intrinsics(1,0,0,0,1,0,0,0,1),
    width(0), height(0),
    position(0,0,0)
{
}



MPImage::
~MPImage(void)
{
  // Remove from panorama and house
  if (panorama) panorama->RemoveImage(this);
  if (house) house->RemoveImage(this);

  // Delete names
  if (name) free(name);
}



void MPImage::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && panorama && panorama->region)
    LoadColor(panorama->region->label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_IMAGE))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_PANORAMA) && panorama)
    LoadColor(panorama->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION) && panorama->region)
    LoadColor(panorama->region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && panorama && panorama->region && panorama->region->level)
    LoadColor(panorama->region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_IMAGE_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_IMAGES] && draw_flags[MP_DRAW_DEPICTIONS]) DrawCamera(draw_flags);
  if (draw_flags[MP_SHOW_IMAGES] && draw_flags[MP_DRAW_BBOXES]) DrawBBox(draw_flags);
  if (draw_flags[MP_SHOW_IMAGES] && draw_flags[MP_DRAW_FACES]) DrawQuads(draw_flags);
  if (draw_flags[MP_SHOW_IMAGES] && draw_flags[MP_DRAW_VERTICES]) DrawPoints(draw_flags);
  if (draw_flags[MP_SHOW_IMAGES] && draw_flags[MP_DRAW_IMAGES]) DrawImage(draw_flags);
}



void MPImage::
DrawCamera(RNFlags draw_flags) const
{
#if 0
  // Determine camera parameters in world coordinates
  R4Matrix camera_to_world = extrinsics.Inverse();
  R3Point eye = camera_to_world * R3zero_point;
  R3Vector towards = camera_to_world * R3negz_vector;
  R3Vector up = camera_to_world * R3posy_vector;

  // Draw vectors along view directions
  glDisable(GL_LIGHTING);
  R3Span(eye, eye+0.5*towards).Draw();
  R3Span(eye, eye+0.3*up).Draw();
#else
  rgbd.DrawCamera(0);
#endif
}



void MPImage::
DrawBBox(RNFlags draw_flags) const
{
  // Draw bounding box
  glDisable(GL_LIGHTING);
  rgbd.DrawBBox(0);
}



void MPImage::
DrawPoints(RNFlags draw_flags) const
{
  // Draw points
  glDisable(GL_LIGHTING);
  rgbd.DrawPoints(RGBD_PHOTO_COLOR_SCHEME);
}



void MPImage::
DrawQuads(RNFlags draw_flags) const
{
  // Draw points
  glEnable(GL_LIGHTING);
  rgbd.DrawQuads(RGBD_RENDER_COLOR_SCHEME);
}



void MPImage::
DrawImage(RNFlags draw_flags) const
{
  // Draw image
  rgbd.DrawImage(RGBD_PHOTO_COLOR_SCHEME, 0.25);
}



////////////////////////////////////////////////////////////////////////
// PANORAMA MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPPanorama::
MPPanorama(void)
  : house(NULL),
    house_index(-1),
    region(NULL),
    region_index(-1),
    name(NULL),
    images(),
    position(0,0,0)
{
}



MPPanorama::
~MPPanorama(void)
{
  // Remove all images
  while (!images.IsEmpty()) RemoveImage(images.Tail());
  
  // Remove from region and house
  if (region) region->RemovePanorama(this);
  if (house) house->RemovePanorama(this);
}



void MPPanorama::
InsertImage(MPImage *image)
{
  // Insert image
  image->panorama = this;
  image->panorama_index = images.NEntries();
  images.Insert(image);
}



void MPPanorama::
RemoveImage(MPImage *image)
{
  // Remove image
  MPImage *tail = images.Tail();
  tail->panorama_index = image->panorama_index;
  images[image->panorama_index] = tail;
  images.RemoveTail();
  image->panorama = NULL;
  image->panorama_index = -1;
}



void MPPanorama::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && region)
    LoadColor(region->label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_PANORAMA))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION) && region)
    LoadColor(region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && region && region->level)
    LoadColor(region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_PANORAMA_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_PANORAMAS] && draw_flags[MP_DRAW_DEPICTIONS]) DrawPosition(draw_flags);
  if (draw_flags[MP_SHOW_PANORAMAS] && draw_flags[MP_DRAW_LABELS]) DrawName(draw_flags);

  // Draw contents
  // DrawImages(draw_flags);
}



void MPPanorama::
DrawPosition(RNFlags draw_flags) const
{
  // Draw sphere at position
  if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
  else glEnable(GL_LIGHTING);
  R3Sphere(position, 0.2).Draw();
}


  
void MPPanorama::
DrawName(RNFlags draw_flags) const
{
  // Draw name
  if (!name) return;
  glDisable(GL_LIGHTING);
  glColor3d(1,1,1);
  DrawText(position + 0.25 * R3posz_vector, name);
}


  
void MPPanorama::
DrawImages(RNFlags draw_flags) const
{
  // Draw all images
  for (int i = 0; i < images.NEntries(); i++) {
    MPImage *image = images.Kth(i);
    image->Draw(draw_flags);
  }
}



////////////////////////////////////////////////////////////////////////
// SEGMENT MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPSegment::
MPSegment(void)
  : house(NULL),
    house_index(-1),
    object(NULL),
    object_index(-1),
    mesh(NULL),
    faces(),
    area(0),
    position(0,0,0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
    id(-1)
{
}



MPSegment::
~MPSegment(void)
{
  // Remove from object and house
  if (object) object->RemoveSegment(this);
  if (house) house->RemoveSegment(this);
}



void MPSegment::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_OBJECT) && object && object->category)
    LoadColor(object->category->mpcat40_id);
  else if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_REGION) && object && object->region && object->region->label)
    LoadColor(object->region->label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_SEGMENT))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_OBJECT) && object)
    LoadColor(object->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION) && object->region)
    LoadColor(object->region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && object && object->region && object->region->level)
    LoadColor(object->region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_SEGMENT_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_SEGMENTS] && draw_flags[MP_DRAW_FACES | MP_DRAW_EDGES | MP_DRAW_VERTICES]) DrawMesh(draw_flags);
  if (draw_flags[MP_SHOW_SEGMENTS] && draw_flags[MP_DRAW_BBOXES]) DrawBBox(draw_flags);
}



void MPSegment::
DrawMesh(RNFlags draw_flags) const
{
  // Check mesh
  if (!mesh) return;

  // Draw faces
  if (draw_flags & MP_DRAW_FACES) {
    if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
    else glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < faces.NEntries(); i++) {
      R3MeshFace *face = faces.Kth(i);
      R3LoadNormal(mesh->FaceNormal(face));
      for (int j = 0; j < 3; j++) {
        R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
        if (draw_flags[MP_COLOR_BY_RGB]) R3LoadRgb(mesh->VertexColor(vertex));
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
    }
    glEnd();
  }

  // Draw edges
  if (draw_flags & MP_DRAW_EDGES) {
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < faces.NEntries(); i++) {
      R3MeshFace *face = faces.Kth(i);
      for (int j = 0; j < 3; j++) {
        R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
        if (draw_flags[MP_COLOR_BY_RGB]) R3LoadRgb(mesh->VertexColor(vertex));
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
    }
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // Draw vertices
  if (draw_flags & MP_DRAW_VERTICES) {
    R3mesh_mark++;
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);
    for (int i = 0; i < faces.NEntries(); i++) {
      R3MeshFace *face = faces.Kth(i);
      for (int j = 0; j < 3; j++) {
        R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
        if (mesh->VertexMark(vertex) == R3mesh_mark) continue;
        mesh->SetVertexMark(vertex, R3mesh_mark);
        if (draw_flags[MP_COLOR_BY_RGB]) R3LoadRgb(mesh->VertexColor(vertex));
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
    }
    glEnd();
  }
}



void MPSegment::
DrawBBox(RNFlags draw_flags) const
{
  // Draw bounding box
  glDisable(GL_LIGHTING);
  bbox.Outline();
}



////////////////////////////////////////////////////////////////////////
// OBJECT MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPObject::
MPObject(void)
  : house(NULL),
    house_index(-1),
    region(NULL),
    region_index(-1),
    category(NULL),
    segments(),
    position(0,0,0),
    obb(R3Point(0.0, 0.0, 0.0), R3Vector(1.0, 0.0, 0.0), R3Vector(0.0, 1.0, 0.0), -1.0, -1.0, -1.0)
{
}



MPObject::
~MPObject(void)
{
  // Remove all segments
  while (!segments.IsEmpty()) RemoveSegment(segments.Tail());
  
  // Remove from category, region and house
  if (category) category->RemoveObject(this);
  if (region) region->RemoveObject(this);
  if (house) house->RemoveObject(this);
}



void MPObject::
InsertSegment(MPSegment *segment)
{
  // Insert segment
  segment->object = this;
  segment->object_index = segments.NEntries();
  segments.Insert(segment);
}



void MPObject::
RemoveSegment(MPSegment *segment)
{
  // Remove segment
  MPSegment *tail = segments.Tail();
  tail->object_index = segment->object_index;
  segments[segment->object_index] = segments.Tail();
  segments.RemoveTail();
  segment->object = NULL;
  segment->object_index = -1;
}



void MPObject::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_OBJECT) && category)
    LoadColor(category->mpcat40_id);
  else if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_REGION) && region && region->label)
    LoadColor(region->label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_OBJECT))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION) && region)
    LoadColor(region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && region && region->level)
    LoadColor(region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_OBJECT_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_OBJECTS] && draw_flags[MP_DRAW_BBOXES]) DrawBBox(draw_flags);
  if (draw_flags[MP_SHOW_OBJECTS] && draw_flags[MP_DRAW_LABELS]) DrawLabel(draw_flags);

  // Draw contents
  DrawSegments(draw_flags & ~(MP_DRAW_BBOXES | MP_DRAW_LABELS));
}



void MPObject::
DrawBBox(RNFlags draw_flags) const
{
  // Draw the oriented box
  glDisable(GL_LIGHTING);
  obb.Outline();
}



void MPObject::
DrawLabel(RNFlags draw_flags) const
{
  // Draw sphere at position
  if (!category) return;
  if (!category->mpcat40_name) return;
  glDisable(GL_LIGHTING);
  glColor3d(1,1,1);
  DrawText(position + 0.25 * R3posz_vector, category->mpcat40_name);
}


  
void MPObject::
DrawSegments(RNFlags draw_flags) const
{
  // Draw segments
  for (int i = 0; i < segments.NEntries(); i++) {
    MPSegment *segment = segments.Kth(i);
    segment->Draw(draw_flags);
  }
}



////////////////////////////////////////////////////////////////////////
// CATEGORY MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPCategory::
MPCategory(void)
  : house(NULL),
    house_index(-1),
    objects(),
    label_id(-1),
    label_name(NULL),
    mpcat40_id(-1),
    mpcat40_name(NULL)
{
}



MPCategory::
~MPCategory(void)
{
  // Remove all objects
  while (!objects.IsEmpty()) RemoveObject(objects.Tail());
  
  // Remove from house
  if (house) house->RemoveCategory(this);
}



void MPCategory::
InsertObject(MPObject *object)
{
  // Insert object
  object->category = this;
  object->category_index = objects.NEntries();
  objects.Insert(object);
}



void MPCategory::
RemoveObject(MPObject *object)
{
  // Remove object
  MPObject *tail = objects.Tail();
  tail->category_index = object->category_index;
  objects[object->category_index] = tail;
  objects.RemoveTail();
  object->category = NULL;
  object->category_index = -1;
}



////////////////////////////////////////////////////////////////////////
// VERTEX MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPVertex::
MPVertex(void)
  : house(NULL),
    house_index(-1),
    surface(NULL),
    surface_index(-1),
    position(0,0,0),
    normal(0,0,0),
    label(NULL)
{
}



MPVertex::
~MPVertex(void)
{
  // Remove from surface and house
  if (surface) surface->RemoveVertex(this);
  if (house) house->RemoveVertex(this);

  // Delete label
  if (label) free(label);
}



void MPVertex::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_SURFACE) && surface)
    LoadColor(surface->label);
  else if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_REGION) && surface && surface->region)
    LoadColor(surface->region->label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_SURFACE) && surface)
    LoadColor(surface->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION) && surface && surface->region)
    LoadColor(surface->region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && surface && surface->region && surface->region->level)
    LoadColor(surface->region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_VERTEX_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_VERTICES] && draw_flags[MP_DRAW_DEPICTIONS]) DrawPosition(draw_flags);
}



void MPVertex::
DrawPosition(RNFlags draw_flags) const
{
  // Draw sphere at position
  if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
  else glEnable(GL_LIGHTING);
  R3Sphere(position, 0.1).Draw();
}



////////////////////////////////////////////////////////////////////////
// SURFACE MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPSurface::
MPSurface(void)
  : house(NULL),
    house_index(-1),
    region(NULL),
    region_index(-1),
    vertices(),
    position(0,0,0),
    normal(0,0,0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
    label(NULL)
{
}



MPSurface::
~MPSurface(void)
{
  // Delete all vertices
  while (vertices.NEntries() > 0) delete vertices.Tail();

  // Remove from region and house
  if (region) region->RemoveSurface(this);
  if (house) house->RemoveSurface(this);

  // Delete label
  if (label) free(label);
}



void MPSurface::
InsertVertex(MPVertex *vertex, RNBoolean search_for_best_index)
{
  // Find location for vertex
  int best_index = vertices.NEntries();
  if (search_for_best_index) {
    RNLength best_d = FLT_MAX;
    for (int i = 0; i < vertices.NEntries(); i++) {
      MPVertex *v0 = vertices.Kth(i);
      MPVertex *v1 = vertices.Kth((i+1)%vertices.NEntries());
      R3Span span(v0->position, v1->position);
      RNLength d = R3Distance(vertex->position, span);
      if (d < best_d) { best_d = d; best_index = (i+1)%vertices.NEntries(); }
    }
  }

  // Insert vertex
  vertex->surface = this;
  vertex->surface_index = best_index;
  vertices.InsertKth(vertex, best_index);

  // Recompute bounding box
  RecomputeBBox();
}



void MPSurface::
RemoveVertex(MPVertex *vertex)
{
  // Remove vertex
  vertices.Remove(vertex);
  vertex->surface = NULL;
  vertex->surface_index = -1;

  // Recompute bounding box
  RecomputeBBox();
}



void MPSurface::
FlipOrientation(void)
{
  // Reverse order of vertices
  vertices.Reverse();
  for (int i = 0; i < vertices.NEntries(); i++) {
    MPVertex *vertex = vertices.Kth(i);
    vertex->surface_index = i;
  }

  // Reverse direction of normal (normals should always be up)
  // normal.Flip();
}



void MPSurface::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_REGION) && region)
    LoadColor(region->label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_SURFACE))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION) && region)
    LoadColor(region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && region && region->level)
    LoadColor(region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_SURFACE_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_SURFACES] && draw_flags[MP_DRAW_FACES | MP_DRAW_EDGES | MP_DRAW_VERTICES]) DrawPolygon(draw_flags);
}



void MPSurface::
DrawPolygon(RNFlags draw_flags) const
{
  // Draw faces
  if (draw_flags & MP_DRAW_FACES) {
    if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
    else glEnable(GL_LIGHTING);
    R3LoadNormal(normal);
    GLUtesselator *tess = gluNewTess();
    gluTessCallback(tess, GLU_TESS_BEGIN, (void (*)()) glBegin);
    gluTessCallback(tess, GLU_TESS_VERTEX, (void (*)()) glVertex3dv);
    gluTessCallback(tess, GLU_TESS_END, (void (*)()) glEnd);
    gluTessBeginPolygon(tess, NULL);
    gluTessBeginContour(tess);
    for (int i = 0; i < vertices.NEntries(); i++) {
      GLdouble *coords = (GLdouble *) vertices[i]->position.Coords();
      gluTessVertex(tess, coords, coords);
    }
    gluTessEndContour(tess);
    gluTessEndPolygon(tess);
    gluDeleteTess(tess);
  }

  // Draw edges
  if (draw_flags & MP_DRAW_EDGES) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < vertices.NEntries(); i++) {
      MPVertex *vertex = vertices.Kth(i);
      R3LoadPoint(vertex->position);
    }
    glEnd();
  }

  // Draw vertices
  if (draw_flags & MP_DRAW_VERTICES) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      MPVertex *vertex = vertices.Kth(i);
      vertex->DrawPosition(draw_flags);
    }
  }
}



void MPSurface::
RecomputeBBox(void)
{
  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < vertices.NEntries(); i++) {
    MPVertex *vertex = vertices.Kth(i);
    vertex->surface_index = i;
    bbox.Union(vertex->position);
  }

  // Update region
  if (region) region->RecomputeBBox(TRUE);
}


////////////////////////////////////////////////////////////////////////
// REGION MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPRegion::
MPRegion(void)
  : house(NULL),
    house_index(-1),
    level(NULL),
    level_index(-1),
    panoramas(),
    surfaces(),
    objects(),
    portals(),
    position(0,0,0),
    height(0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
    label(NULL)
{
}



MPRegion::
~MPRegion(void)
{
  // Remove panoramas, surfaces, and portals
  while (panoramas.NEntries() > 0) RemovePanorama(panoramas.Tail());
  while (surfaces.NEntries() > 0) RemoveSurface(surfaces.Tail());
  while (objects.NEntries() > 0) RemoveObject(objects.Tail());
  while (portals.NEntries() > 0) RemovePortal(portals.Tail());

  // Remove from level and house
  if (level) level->RemoveRegion(this);
  if (house) house->RemoveRegion(this);

  // Delete label
  if (label) free(label);
}



R2Polygon MPRegion::
FloorPolygon(void) const
{
  // Create wall grid
  RNArray<R2Point *> points;
  for (int i = 0; i < surfaces.NEntries(); i++) {
    MPSurface *surface = surfaces.Kth(i);
    if (surface->normal.Z() < 0.9) continue;
    for (int j = 0; j < surface->vertices.NEntries(); j++) 
      points.Insert((R2Point *) &(surface->vertices[j]->position));
    return R2Polygon(points);
  }

  // If no surface found, return empty polygon
  return R2Polygon();
}



void MPRegion::
InsertPanorama(MPPanorama *panorama)
{
  // Insert panorama
  panorama->region = this;
  panorama->region_index = panoramas.NEntries();
  panoramas.Insert(panorama);

  // Update bounding box
  bbox.Union(panorama->position);
}



void MPRegion::
RemovePanorama(MPPanorama *panorama)
{
  // Remove panorama
  MPPanorama *tail = panoramas.Tail();
  tail->region_index = panorama->region_index;
  panoramas[panorama->region_index] = tail;
  panoramas.RemoveTail();
  panorama->region = NULL;
  panorama->region_index = -1;
}



void MPRegion::
InsertSurface(MPSurface *surface)
{
  // Insert surface
  surface->region = this;
  surface->region_index = surfaces.NEntries();
  surfaces.Insert(surface);

  // Update bounding box
  bbox.Union(surface->bbox);
}



void MPRegion::
RemoveSurface(MPSurface *surface)
{
  // Remove surface
  MPSurface *tail = surfaces.Tail();
  tail->region_index = surface->region_index;
  surfaces[surface->region_index] = tail;
  surfaces.RemoveTail();
  surface->region = NULL;
  surface->region_index = -1;
}



void MPRegion::
InsertObject(MPObject *object)
{
  // Insert object
  object->region = this;
  object->region_index = objects.NEntries();
  objects.Insert(object);

  // Update bounding box
  bbox.Union(object->obb.BBox());
}



void MPRegion::
RemoveObject(MPObject *object)
{
  // Remove object
  MPObject *tail = objects.Tail();
  tail->region_index = object->region_index;
  objects[object->region_index] = tail;
  objects.RemoveTail();
  object->region = NULL;
  object->region_index = -1;
}



void MPRegion::
InsertPortal(MPPortal *portal, int index)
{
  // Insert portal
  portal->regions[index] = this;
  portal->region_indices[index] = portals.NEntries();
  portals.Insert(portal);

  // Update bounding box
  bbox.Union(portal->span.BBox());
}



void MPRegion::
RemovePortal(MPPortal *portal)
{
  // Find index
  int index = -1;
  if (portal->regions[0] == this) index = 0;
  else if (portal->regions[1] == this) index = 1;
  else RNAbort("portal not in region");
  
  // Remove portal
  MPPortal *tail = portals.Tail();
  tail->region_indices[index] = portal->region_indices[index];
  portals[portal->region_indices[index]] = tail;
  portals.RemoveTail();
  portal->regions[index] = NULL;
  portal->region_indices[index] = -1;
}



void MPRegion::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_LABEL) && (draw_flags & MP_COLOR_BY_REGION))
    LoadColor(label);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_REGION))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL) && level)
    LoadColor(level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_REGION_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_REGIONS] && draw_flags[MP_DRAW_DEPICTIONS])
    DrawPosition(draw_flags);
  if (draw_flags[MP_SHOW_REGIONS] && draw_flags[MP_DRAW_BBOXES])
    DrawBBox(draw_flags);
  if (draw_flags[MP_SHOW_REGIONS] && draw_flags[MP_DRAW_LABELS])
    DrawLabel(draw_flags);

  // Draw contents
  // DrawPanoramas();
  DrawSurfaces(draw_flags & ~(MP_DRAW_BBOXES | MP_DRAW_DEPICTIONS));
  // DrawObjects();

}



void MPRegion::
DrawPosition(RNFlags draw_flags) const
{
  // Draw a sphere at position
  if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
  else glEnable(GL_LIGHTING);
  R3Sphere(position, 0.2).Draw(R3_SURFACES_DRAW_FLAG);
}


 
void MPRegion::
DrawBBox(RNFlags draw_flags) const
{
  // Draw bounding box
  glDisable(GL_LIGHTING);
  bbox.Outline();
}



void MPRegion::
DrawLabel(RNFlags draw_flags) const
{
  // Draw sphere at position
  if (!label) return;
  glDisable(GL_LIGHTING);
  glColor3d(1,1,1);
  DrawText(position + 0.25 * R3posz_vector, label);
}


  
void MPRegion::
DrawPanoramas(RNFlags draw_flags) const
{
  // Draw panoramas
  for (int i = 0; i < panoramas.NEntries(); i++) {
    MPPanorama *panorama = panoramas.Kth(i);
    panorama->Draw(draw_flags);
  }
}



void MPRegion::
DrawSurfaces(RNFlags draw_flags) const
{
  // Draw surfaces
  for (int i = 0; i < surfaces.NEntries(); i++) {
    MPSurface *surface = surfaces.Kth(i);
    surface->Draw(draw_flags);
  }
}



void MPRegion::
DrawObjects(RNFlags draw_flags) const
{
  // Draw objects
  for (int i = 0; i < objects.NEntries(); i++) {
    MPObject *object = objects.Kth(i);
    object->Draw(draw_flags);
  }
}



void MPRegion::
DrawPortals(RNFlags draw_flags) const
{
  // Draw portals
  for (int i = 0; i < portals.NEntries(); i++) {
    MPPortal *portal = portals.Kth(i);
    portal->Draw(draw_flags);
  }
}



void MPRegion::
RecomputeBBox(RNBoolean preserve_zmax)
{
  // Remember height
  RNScalar zmax = bbox[1][2];

  // Update bbox
  bbox = R3null_box;
  for (int i = 0; i < surfaces.NEntries(); i++) {
    MPSurface *surface = surfaces.Kth(i);
    bbox.Union(surface->bbox);
  }
  
  // Restore zextent
  if (preserve_zmax) bbox[1][2] = zmax;
}



////////////////////////////////////////////////////////////////////////
// PORTAL MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPPortal::
MPPortal(void)
  : house(NULL),
    house_index(-1),
    span(R3Point(0,0,0), R3Point(0,0,0)),
    label(NULL)
{
  // Initialize regions
  regions[0] = regions[1] = NULL;
  region_indices[0] = region_indices[1] = -1;
}



MPPortal::
~MPPortal(void)
{
  // Remove from region and house
  if (regions[0]) regions[0]->RemovePortal(this);
  if (regions[1]) regions[1]->RemovePortal(this);
  if (house) house->RemovePortal(this);

  // Delete label
  if (label) free(label);
}



void MPPortal::
Draw(RNFlags draw_flags) const
{
  // Set the color
  MPRegion *region = regions[0];
  if ((draw_flags & MP_COLOR_BY_PORTAL) && (draw_flags & MP_COLOR_BY_LABEL))
    LoadColor(label);
  else if ((draw_flags & MP_COLOR_BY_PORTAL) && (draw_flags & MP_COLOR_BY_INDEX))
    LoadColor(house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_REGION) && (draw_flags & MP_COLOR_BY_INDEX) && region)
    LoadColor(region->house_index + 1);
  else if ((draw_flags & MP_COLOR_BY_REGION) && (draw_flags & MP_COLOR_BY_LABEL) && region)
    LoadColor(region->label);
  else if ((draw_flags & MP_COLOR_BY_LEVEL) && (draw_flags & MP_COLOR_BY_INDEX) && region && region->level)
    LoadColor(region->level->house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_PORTAL_TAG);

  // Draw this
  if (draw_flags[MP_SHOW_PORTALS] && draw_flags[MP_DRAW_DEPICTIONS]) DrawSpan(draw_flags);
  if (draw_flags[MP_SHOW_PORTALS] && draw_flags[MP_DRAW_LABELS]) DrawLabel(draw_flags);
}



void MPPortal::
DrawSpan(RNFlags draw_flags) const
{
  // Draw span
  glDisable(GL_LIGHTING);
  span.Draw();

  // Draw spheres
  R3Sphere(span.Start(), 0.1).Draw();  
  R3Sphere(span.End(), 0.1).Draw();  
}



void MPPortal::
DrawLabel(RNFlags draw_flags) const
{
  // Draw label above midpoint
  if (!label) return;
  glDisable(GL_LIGHTING);
  glColor3d(1, 1, 1);
  DrawText(span.Midpoint() + 0.5 * R3posz_vector, label);
}


  
////////////////////////////////////////////////////////////////////////
// LEVEL MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPLevel::
MPLevel(void)
  : house(NULL),
    house_index(-1),
    regions(),
    position(0,0,0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
    label(NULL)  
{
}



MPLevel::
~MPLevel(void)
{
  // Remove regions
  while (regions.NEntries() > 0) RemoveRegion(regions.Tail());

  // Remove from house
  if (house) house->RemoveLevel(this);

  // Delete label
  if (label) free(label);
}



void MPLevel::
InsertRegion(MPRegion *region)
{
  // Insert region
  region->level = this;
  region->level_index = regions.NEntries();
  regions.Insert(region);

  // Update bounding box
  bbox.Union(region->bbox);
}



void MPLevel::
RemoveRegion(MPRegion *region)
{
  // Remove region
  MPRegion *tail = regions.Tail();
  tail->level_index = region->level_index;
  regions[region->level_index] = tail;
  regions.RemoveTail();
  region->level = NULL;
  region->level_index = -1;
}



void MPLevel::
Draw(RNFlags draw_flags) const
{
  // Set the color
  if ((draw_flags & MP_COLOR_BY_INDEX) && (draw_flags & MP_COLOR_BY_LEVEL))
    LoadColor(house_index + 1);
  else if (draw_flags & MP_COLOR_FOR_PICK)
    LoadIndex(house_index, MP_LEVEL_TAG);
  else glColor3d(0.1, 0.5, 0.9);

  // Draw this
  if (draw_flags[MP_SHOW_LEVELS] && draw_flags[MP_DRAW_DEPICTIONS])
    DrawPosition(draw_flags);

  // Draw contents
  // DrawRegions(draw_flags);
}



void MPLevel::
DrawPosition(RNFlags draw_flags) const
{
  // Draw a sphere at position
  if (draw_flags[MP_DRAW_DEPICTIONS]) {
    if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
    else glEnable(GL_LIGHTING);
    R3Sphere(position, 0.2).Draw(R3_SURFACES_DRAW_FLAG);
  }
}


 
void MPLevel::
DrawBBox(RNFlags draw_flags) const
{
  // Draw bounding box
  if (draw_flags & MP_DRAW_BBOXES) {
    glDisable(GL_LIGHTING);
    bbox.Outline();
  }
}



void MPLevel::
DrawRegions(RNFlags draw_flags) const
{
  // Draw regions
  for (int i = 0; i < regions.NEntries(); i++) {
    MPRegion *region = regions.Kth(i);
    region->Draw(draw_flags);
  }
}



////////////////////////////////////////////////////////////////////////
// HOUSE MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

MPHouse::
MPHouse(const char *name, const char *label)
  : images(),
    panoramas(),
    categories(),
    segments(),
    objects(),
    vertices(),
    surfaces(),
    regions(),
    portals(),
    levels(),
    rgbd(),
    scene(NULL),
    mesh(NULL),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),    
    label((label) ? strdup(label) : NULL),
    name((name) ? strdup(name) : NULL)    
{
  // Set RGBD configuration parameters
  rgbd.SetDatasetFormat("matterport");
  rgbd.SetDepthDirectory("undistorted_depth_images");
  rgbd.SetColorDirectory("undistorted_color_images");
}



MPHouse::
~MPHouse(void)
{
  // Delete everything
  while (!images.IsEmpty()) delete images.Tail();
  while (!panoramas.IsEmpty()) delete panoramas.Tail();
  while (!objects.IsEmpty()) delete objects.Tail();
  while (!segments.IsEmpty()) delete segments.Tail();
  while (!categories.IsEmpty()) delete categories.Tail();
  while (!vertices.IsEmpty()) delete vertices.Tail();
  while (!surfaces.IsEmpty()) delete surfaces.Tail();
  while (!regions.IsEmpty()) delete regions.Tail();
  while (!portals.IsEmpty()) delete portals.Tail();
  while (!levels.IsEmpty()) delete levels.Tail();

  // Delete scene and mesh
  if (scene) delete scene;
  if (mesh) delete mesh;

  // Delete label and name
  if (label) free(label);
  if (name) free(name);
}



void MPHouse::
InsertImage(MPImage *image)
{
  // Insert image
  image->house = this;
  image->house_index = images.NEntries();
  images.Insert(image);

  // Insert rgbd image
  rgbd.InsertImage(&(image->rgbd));

  // Update bounding box
  bbox.Union(image->position);
}



void MPHouse::
RemoveImage(MPImage *image)
{
  // Remove image
  MPImage *tail = images.Tail();
  tail->house_index = image->house_index;
  images[image->house_index] = tail;
  images.RemoveTail();
  image->house = NULL;
  image->house_index = -1;

  // Remove rgbd image
  rgbd.RemoveImage(&(image->rgbd));
}



void MPHouse::
InsertPanorama(MPPanorama *panorama)
{
  // Insert panorama
  panorama->house = this;
  panorama->house_index = panoramas.NEntries();
  panoramas.Insert(panorama);

  // Update bounding box
  bbox.Union(panorama->position);
}



void MPHouse::
RemovePanorama(MPPanorama *panorama)
{
  // Remove panorama
  MPPanorama *tail = panoramas.Tail();
  tail->house_index = panorama->house_index;
  panoramas[panorama->house_index] = tail;
  panoramas.RemoveTail();
  panorama->house = NULL;
  panorama->house_index = -1;
}



void MPHouse::
InsertCategory(MPCategory *category)
{
  // Insert category
  category->house = this;
  category->house_index = categories.NEntries();
  categories.Insert(category);
}



void MPHouse::
RemoveCategory(MPCategory *category)
{
  // Remove category
  MPCategory *tail = categories.Tail();
  tail->house_index = category->house_index;
  categories[category->house_index] = categories.Tail();
  categories.RemoveTail();
  category->house = NULL;
  category->house_index = -1;
}



void MPHouse::
InsertSegment(MPSegment *segment)
{
  // Insert segment
  segment->house = this;
  segment->house_index = segments.NEntries();
  segments.Insert(segment);
}



void MPHouse::
RemoveSegment(MPSegment *segment)
{
  // Remove segment
  MPSegment *tail = segments.Tail();
  tail->house_index = segment->house_index;
  segments[segment->house_index] = segments.Tail();
  segments.RemoveTail();
  segment->house = NULL;
  segment->house_index = -1;
}



void MPHouse::
InsertObject(MPObject *object)
{
  // Insert object
  object->house = this;
  object->house_index = objects.NEntries();
  objects.Insert(object);

  // Update bounding box
  bbox.Union(object->obb.BBox());
}



void MPHouse::
RemoveObject(MPObject *object)
{
  // Remove object
  MPObject *tail = objects.Tail();
  tail->house_index = object->house_index;
  objects[object->house_index] = objects.Tail();
  objects.RemoveTail();
  object->house = NULL;
  object->house_index = -1;
}



void MPHouse::
InsertVertex(MPVertex *vertex)
{
  // Insert vertex
  vertex->house = this;
  vertex->house_index = vertices.NEntries();
  vertices.Insert(vertex);

  // Update bounding box
  bbox.Union(vertex->position);
}



void MPHouse::
RemoveVertex(MPVertex *vertex)
{
  // Remove vertex
  MPVertex *tail = vertices.Tail();
  tail->house_index = vertex->house_index;
  vertices[vertex->house_index] = vertices.Tail();
  vertices.RemoveTail();
  vertex->house = NULL;
  vertex->house_index = -1;
}



void MPHouse::
InsertSurface(MPSurface *surface)
{
  // Insert surface
  surface->house = this;
  surface->house_index = surfaces.NEntries();
  surfaces.Insert(surface);

  // Update bounding box
  bbox.Union(surface->bbox);
}



void MPHouse::
RemoveSurface(MPSurface *surface)
{
  // Remove surface
  MPSurface *tail = surfaces.Tail();
  tail->house_index = surface->house_index;
  surfaces[surface->house_index] = surfaces.Tail();
  surfaces.RemoveTail();
  surface->house = NULL;
  surface->house_index = -1;
}



void MPHouse::
InsertRegion(MPRegion *region)
{
  // Insert region
  region->house = this;
  region->house_index = regions.NEntries();
  regions.Insert(region);

  // Update bounding box
  bbox.Union(region->bbox);
}



void MPHouse::
RemoveRegion(MPRegion *region)
{
  // Remove region
  MPRegion *tail = regions.Tail();
  tail->house_index = region->house_index;
  regions[region->house_index] = regions.Tail();
  regions.RemoveTail();
  region->house = NULL;
  region->house_index = -1;
}



void MPHouse::
InsertPortal(MPPortal *portal)
{
  // Insert portal
  portal->house = this;
  portal->house_index = portals.NEntries();
  portals.Insert(portal);

  // Update bounding box
  bbox.Union(portal->span.BBox());
}



void MPHouse::
RemovePortal(MPPortal *portal)
{
  // Remove portal
  MPPortal *tail = portals.Tail();
  tail->house_index = portal->house_index;
  portals[portal->house_index] = portals.Tail();
  portals.RemoveTail();
  portal->house = NULL;
  portal->house_index = -1;
}



void MPHouse::
InsertLevel(MPLevel *level)
{
  // Insert level
  level->house = this;
  level->house_index = levels.NEntries();
  levels.Insert(level);

  // Update bounding box
  bbox.Union(level->bbox);
}



void MPHouse::
RemoveLevel(MPLevel *level)
{
  // Remove level
  MPLevel *tail = levels.Tail();
  tail->house_index = level->house_index;
  levels[level->house_index] = levels.Tail();
  levels.RemoveTail();
  level->house = NULL;
  level->house_index = -1;
}



MPSegment *MPHouse::
FindSegment(int id) const
{
  // Check each segment
  for (int i = 0; i < segments.NEntries(); i++) {
    MPSegment *segment = segments.Kth(i);
    if (segment->id == id) return segment;
  }

  // Did not find segment
  return NULL;
}



MPCategory *MPHouse::
FindCategory(int label_id) const
{
  // Check each category
  for (int i = 0; i < categories.NEntries(); i++) {
    MPCategory *category = categories.Kth(i);
    if (category->label_id == label_id) return category;
  }

  // Did not find category
  return NULL;
}



MPCategory *MPHouse::
FindCategory(const char *label_name) const
{
  // Check each category
  for (int i = 0; i < categories.NEntries(); i++) {
    MPCategory *category = categories.Kth(i);
    if (!strcmp(category->label_name, label_name)) return category;
  }

  // Did not find category
  return NULL;
}



MPRegion *MPHouse::
FindRegion(const R3Point& position, RNLength max_distance) const
{
  // Check each region
  MPRegion *best_region = NULL;
  RNLength best_distance = max_distance;
  for (int i = 0; i < regions.NEntries(); i++) {
    MPRegion *region = regions.Kth(i);
    RNLength bbox_distance = R3Distance(region->bbox, position);
    if (bbox_distance > max_distance) continue;
    R2Polygon polygon = region->FloorPolygon();
    RNLength polygon_distance = R2Distance(polygon, R2Point(position.X(), position.Y()));
    if (polygon_distance < best_distance) {
      best_distance = polygon_distance;
      best_region = region;
    }
  }

  // Return best region
  return best_region;
}



MPLevel *MPHouse::
FindLevel(const R3Point& position, RNLength max_dz) const
{
  // Check every level
  MPLevel *best_level = NULL;
  for (int i = 0; i < levels.NEntries(); i++) {
    MPLevel *level = levels.Kth(i);
    RNScalar dz = fabs(level->position.Z() - position.Z());
    if (dz < max_dz) { best_level = level; max_dz = dz; }
  }

  // Return best level
  return best_level;
}



MPImage *MPHouse::
FindClosestImage(const R3Point& query_position, const R3Vector& query_normal,
  RNLength max_distance, RNBoolean check_normal, RNBoolean check_visibility) const
{
  // Check all images
  MPImage *closest_image = NULL;
  RNLength closest_d = max_distance;
  R3Halfspace query_halfspace(query_position, query_normal);
  for (int i = 0; i < images.NEntries(); i++) {
    MPImage *image = images.Kth(i);
    RNLength d = R3Distance(query_position, image->position);
    if (d >= closest_d) continue;
    if (check_normal) {
      if (!R3Intersects(query_halfspace, image->position)) continue;
    }
    closest_image = image;
    closest_d = d;
  }
  
  // Return closest image
  return closest_image;
}



MPVertex *MPHouse::
FindClosestVertex(const R3Point& query_position, const R3Vector& query_normal,
  RNLength max_distance, RNBoolean check_normal, RNBoolean check_visibility) const
{
  // Check all vertices
  MPVertex *closest_vertex = NULL;
  RNLength closest_d = max_distance;
  for (int i = 0; i < vertices.NEntries(); i++) {
    MPVertex *vertex = vertices.Kth(i);
    RNLength d = R3Distance(query_position, vertex->position);
    if (d >= closest_d) continue;
    if (check_normal) {
      if (query_normal.Dot(vertex->normal) < 0) continue;
    }
    closest_vertex = vertex;
    closest_d = d;
  }
  
  // Return closest vertex
  return closest_vertex;
}



MPRegion *MPHouse::
FindClosestRegion(const R3Point& query_position, const R3Vector& query_normal,
  RNLength max_distance, RNBoolean check_normal, RNBoolean check_visibility) const
{
  // Check all regions
  MPRegion *closest_region = NULL;
  RNLength closest_d = max_distance;
  R3Halfspace query_halfspace(query_position, query_normal);
  for (int i = 0; i < regions.NEntries(); i++) {
    MPRegion *region = regions.Kth(i);
    RNLength d = R3Distance(query_position, region->position);
    if (d >= closest_d) continue;
    if (check_normal) {
      if (!R3Intersects(query_halfspace, region->position)) continue;
    }
    closest_region = region;
    closest_d = d;
  }
  
  // Return closest region
  return closest_region;
}



void MPHouse::
Draw(RNFlags draw_flags) const
{
  // Draw contents
  DrawImages(draw_flags);
  DrawPanoramas(draw_flags);
  DrawObjects(draw_flags);
  DrawRegions(draw_flags);
  DrawPortals(draw_flags);
  DrawLevels(draw_flags);
  DrawMesh(draw_flags);
  DrawScene(draw_flags);
}



void MPHouse::
DrawLevels(RNFlags draw_flags) const
{
  // Draw levels
  for (int i = 0; i < levels.NEntries(); i++) {
    MPLevel *level = levels.Kth(i);
    level->Draw(draw_flags);
  }
}



void MPHouse::
DrawPortals(RNFlags draw_flags) const
{
  // Draw portals
  for (int i = 0; i < portals.NEntries(); i++) {
    MPPortal *portal = portals.Kth(i);
    portal->Draw(draw_flags);
  }
}



void MPHouse::
DrawRegions(RNFlags draw_flags) const
{
  // Draw regions
  for (int i = 0; i < regions.NEntries(); i++) {
    MPRegion *region = regions.Kth(i);
    region->Draw(draw_flags);
  }
}



void MPHouse::
DrawSurfaces(RNFlags draw_flags) const
{
  // Draw surfaces
  for (int i = 0; i < surfaces.NEntries(); i++) {
    MPSurface *surface = surfaces.Kth(i);
    surface->Draw(draw_flags);
  }
}



void MPHouse::
DrawVertices(RNFlags draw_flags) const
{
  // Draw vertices
  for (int i = 0; i < vertices.NEntries(); i++) {
    MPVertex *vertex = vertices.Kth(i);
    vertex->Draw(draw_flags);
  }
}



void MPHouse::
DrawObjects(RNFlags draw_flags) const
{
  // Draw objects
  for (int i = 0; i < objects.NEntries(); i++) {
    MPObject *object = objects.Kth(i);
    object->Draw(draw_flags);
  }
}



void MPHouse::
DrawSegments(RNFlags draw_flags) const
{
  // Draw segments
  for (int i = 0; i < segments.NEntries(); i++) {
    MPSegment *segment = segments.Kth(i);
    segment->Draw(draw_flags);
  }
}



void MPHouse::
DrawPanoramas(RNFlags draw_flags) const
{
  // Draw panoramas
  for (int i = 0; i < panoramas.NEntries(); i++) {
    MPPanorama *panorama = panoramas.Kth(i);
    panorama->Draw(draw_flags);
  }
}



void MPHouse::
DrawImages(RNFlags draw_flags) const
{
  // Draw images
  for (int i = 0; i < images.NEntries(); i++) {
    MPImage *image = images.Kth(i);
    image->Draw(draw_flags);
  }
}



void MPHouse::
DrawMesh(RNFlags draw_flags) const
{
  // Check mesh
  if (!mesh) return;

  // Check draw flags
  if (!draw_flags[MP_SHOW_MESH]) return;

  // Draw faces
  if (draw_flags[MP_DRAW_FACES]) {
    if (draw_flags[MP_COLOR_FOR_PICK]) {
      glDisable(GL_LIGHTING);
      for (int i = 0; i < mesh->NFaces(); i++) {
        R3MeshFace *face = mesh->Face(i);
        LoadIndex(i, MP_MESH_TAG);
        mesh->DrawFace(face);
      }
    }
    else if (draw_flags[MP_COLOR_BY_LABEL]) {
      if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
      else glEnable(GL_LIGHTING);
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < mesh->NFaces(); i++) {
        R3MeshFace *face = mesh->Face(i);
        int category = mesh->FaceCategory(face);
        LoadColor(category);
        R3LoadNormal(mesh->FaceNormal(face));
        for (int j = 0; j < 3; j++) {
          R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
          R3LoadPoint(mesh->VertexPosition(vertex));
        }
      }
      glEnd();
    }
    else if (draw_flags[MP_COLOR_BY_INDEX]) {
      if (draw_flags & MP_COLOR_FOR_PICK) glDisable(GL_LIGHTING);
      else glEnable(GL_LIGHTING);
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < mesh->NFaces(); i++) {
        R3MeshFace *face = mesh->Face(i);
        int segment = mesh->FaceSegment(face);
        LoadColor(segment);
        R3LoadNormal(mesh->FaceNormal(face));
        for (int j = 0; j < 3; j++) {
          R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
          R3LoadPoint(mesh->VertexPosition(vertex));
        }
      }
      glEnd();
    }
    else {
      glDisable(GL_LIGHTING);
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < mesh->NFaces(); i++) {
        R3MeshFace *face = mesh->Face(i);
        for (int j = 0; j < 3; j++) {
          R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
          R3LoadRgb(mesh->VertexColor(vertex));
          R3LoadPoint(mesh->VertexPosition(vertex));
        }
      }
      glEnd();
    }
  }

  // Draw edges
  if (draw_flags[MP_DRAW_EDGES]) {
    if (!draw_flags[MP_COLOR_FOR_PICK]) {
      glDisable(GL_LIGHTING);
      glColor3d(0.0, 1.0, 0.0);
      mesh->DrawEdges();
    }
  }
    
  // Draw vertices
  if (draw_flags[MP_DRAW_VERTICES]) {
    if (draw_flags[MP_COLOR_FOR_PICK]) {
      glDisable(GL_LIGHTING);
      glBegin(GL_POINTS);
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        R3MeshFace *face = mesh->FaceOnVertex(vertex);
        if (!face) continue;
        LoadIndex(mesh->FaceID(face), MP_MESH_TAG);
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
      glEnd();
    }
    else if (draw_flags[MP_COLOR_BY_LABEL]) {
      glDisable(GL_LIGHTING);
      glBegin(GL_POINTS);
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        R3MeshFace *face = mesh->FaceOnVertex(vertex);
        if (!face) continue;
        int category = mesh->FaceCategory(face);
        LoadColor(category);
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
      glEnd();
    }
    else if (draw_flags[MP_COLOR_BY_INDEX]) {
      glDisable(GL_LIGHTING);
      glBegin(GL_POINTS);
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        R3MeshFace *face = mesh->FaceOnVertex(vertex);
        if (!face) continue;
        int segment = mesh->FaceSegment(face);
        LoadColor(segment);
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
      glEnd();
    }
    else {
      glDisable(GL_LIGHTING);
      glBegin(GL_POINTS);
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        R3LoadRgb(mesh->VertexColor(vertex));
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
      glEnd();
    }
  }
}



void MPHouse::
DrawScene(RNFlags draw_flags) const
{
  // Check scene
  if (!scene) return;

  // Check draw flags
  if (!draw_flags[MP_SHOW_SCENE]) return;

  // Draw faces
  if (draw_flags[MP_DRAW_FACES]) {
    glEnable(GL_LIGHTING);
    glColor3d(1.0, 1.0, 1.0); 
    scene->Draw(R3_DEFAULT_DRAW_FLAGS);
  }

  // Draw edges
  if (draw_flags[MP_DRAW_EDGES]) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    scene->Draw(R3_EDGES_DRAW_FLAG);
  }
}



int MPHouse::
ReadFile(const char *filename)
{
  // Check file type
  char file_type[1024];
  FILE *fp = fopen(filename, "r");
  if (!fp) { fprintf(stderr, "Unable to open %s\n", filename); return 0; }
  fscanf(fp, "%s", file_type);
  fclose(fp);

  // Read file
  if (!strcmp(file_type, "ASCII")) return ReadAsciiFile(filename);
  else fprintf(stderr, "Unable to read file %s, unrecognized type: %s\n", filename, file_type); 

  // Unrecognized file type
  return 0;
}



int MPHouse::
WriteFile(const char *filename) const
{
  // Write file of appropriate type
  return WriteAsciiFile(filename);
}



////////////////////////////////////////////////////////////////////////
// Ascii file parsing
////////////////////////////////////////////////////////////////////////

int MPHouse::
ReadAsciiFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open house file %s\n", filename);
    return 0;
  }

  // Useful variables
  char cmd[1024], version[1024], name_buffer[1024], label_buffer[1024];
  int nimages, npanoramas, nvertices, nsurfaces, nsegments, nobjects, ncategories, nregions, nportals, nlevels;
  int house_index, level_index, region_index, surface_index, category_index, object_index, panorama_index, id, dummy;
  RNScalar height, area;
  R3Point position;
  R3Vector normal;
  R3Box box;

  // Read file type and version
  fscanf(fp, "%s%s", cmd, version);
  if (strcmp(cmd, "ASCII")) {
    fprintf(stderr, "Unable to read ascii file %s, wrong type: %s\n", filename, cmd);
    return 0;
  }

  // Read header
  if (!strcmp(version, "1.0")) {
    nsegments = 0;
    nobjects = 0;
    ncategories = 0;
    nportals = 0;
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%s", name_buffer);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%d", &nimages);
    fscanf(fp, "%d", &npanoramas);
    fscanf(fp, "%d", &nvertices);
    fscanf(fp, "%d", &nsurfaces);
    fscanf(fp, "%d", &nregions);
    fscanf(fp, "%d", &nlevels);
    fscanf(fp, "%lf%lf%lf%lf%lf%lf", &box[0][0], &box[0][1], &box[0][2], &box[1][0], &box[1][1], &box[1][2]);
    for (int i = 0; i < 8; i++) fscanf(fp, "%d", &dummy);
  }
  else {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%s", name_buffer);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%d", &nimages);
    fscanf(fp, "%d", &npanoramas);
    fscanf(fp, "%d", &nvertices);
    fscanf(fp, "%d", &nsurfaces);
    fscanf(fp, "%d", &nsegments);
    fscanf(fp, "%d", &nobjects);
    fscanf(fp, "%d", &ncategories);
    fscanf(fp, "%d", &nregions);
    fscanf(fp, "%d", &nportals);
    fscanf(fp, "%d", &nlevels);
    for (int i = 0; i < 5; i++) fscanf(fp, "%d", &dummy);
    fscanf(fp, "%lf%lf%lf%lf%lf%lf", &box[0][0], &box[0][1], &box[0][2], &box[1][0], &box[1][1], &box[1][2]);
    for (int i = 0; i < 5; i++) fscanf(fp, "%d", &dummy);
  }

  // Fill in house info
  this->name = (strcmp(name_buffer, "-")) ? strdup(name_buffer) : NULL;
  this->label = (strcmp(label_buffer, "-")) ? strdup(label_buffer) : NULL;
  this->bbox = box;

  // Read levels
  for (int i = 0; i < nlevels; i++) {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &dummy);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    fscanf(fp, "%lf%lf%lf%lf%lf%lf", &box[0][0], &box[0][1], &box[0][2], &box[1][0], &box[1][1], &box[1][2]);
    for (int j = 0; j < 5; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "L")) { fprintf(stderr, "Error reading level %d\n", i); return 0; }
    MPLevel *level = new MPLevel();
    level->position = position;
    level->label = (strcmp(label_buffer, "-")) ? strdup(label_buffer) : NULL;
    level->bbox = box;
    InsertLevel(level);
  }
    
  // Read regions
  for (int i = 0; i < nregions; i++) {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &level_index);
    fscanf(fp, "%d%d", &dummy, &dummy);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    fscanf(fp, "%lf%lf%lf%lf%lf%lf", &box[0][0], &box[0][1], &box[0][2], &box[1][0], &box[1][1], &box[1][2]);
    fscanf(fp, "%lf", &height);
    for (int j = 0; j < 4; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "R")) { fprintf(stderr, "Error reading region %d\n", i); return 0; }
    MPRegion *region = new MPRegion();
    region->position = position;
    region->label = (strcmp(label_buffer, "-")) ? strdup(label_buffer) : NULL;
    region->bbox = box;
    region->height = (height > 0) ? height : bbox.ZMax() - position.Z();
    InsertRegion(region);
    if (level_index >= 0) {
      MPLevel *level = levels.Kth(level_index);
      level->InsertRegion(region);
    }
  }
    
  // Read portals
  for (int i = 0; i < nportals; i++) {
    int region0_index, region1_index;
    R3Point p0, p1;
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &region0_index);
    fscanf(fp, "%d", &region1_index);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%lf%lf%lf", &p0[0], &p0[1], &p0[2]);
    fscanf(fp, "%lf%lf%lf", &p1[0], &p1[1], &p1[2]);
    for (int j = 0; j < 4; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "P")) { fprintf(stderr, "Error reading portal %d\n", i); return 0; }
    MPPortal *portal = new MPPortal();
    portal->span.Reset(p0, p1);
    portal->label = (strcmp(label_buffer, "-")) ? strdup(label_buffer) : NULL;
    InsertPortal(portal);
    if (region0_index >= 0) {
      MPRegion *region0 = regions.Kth(region0_index);
      region0->InsertPortal(portal, 0);
    }
    if (region1_index >= 0) {
      MPRegion *region1 = regions.Kth(region1_index);
      region1->InsertPortal(portal, 1);
    }
  }
    
  // Read surfaces
  for (int i = 0; i < nsurfaces; i++) {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &region_index);
    fscanf(fp, "%d", &dummy);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    fscanf(fp, "%lf%lf%lf", &normal[0], &normal[1], &normal[2]);
    fscanf(fp, "%lf%lf%lf%lf%lf%lf", &box[0][0], &box[0][1], &box[0][2], &box[1][0], &box[1][1], &box[1][2]);
    for (int j = 0; j < 5; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "S")) { fprintf(stderr, "Error reading surface %d\n", i); return 0; }
    MPSurface *surface = new MPSurface();
    surface->position = position;
    surface->normal = normal;
    surface->label = (strcmp(label_buffer, "-")) ? strdup(label_buffer) : NULL;
    surface->bbox = box;
    InsertSurface(surface);
    if (region_index >= 0) {
      MPRegion *region = regions.Kth(region_index);
      region->InsertSurface(surface);
    }
  }
    
  // Read vertices
  for (int i = 0; i < nvertices; i++) {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &surface_index);
    fscanf(fp, "%s", label_buffer);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    fscanf(fp, "%lf%lf%lf", &normal[0], &normal[1], &normal[2]);
    for (int j = 0; j < 3; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "V")) { fprintf(stderr, "Error reading vertex %d\n", i); return 0; }
    MPVertex *vertex = new MPVertex();
    vertex->position = position;
    vertex->normal = normal;
    vertex->label = (strcmp(label_buffer, "-")) ? strdup(label_buffer) : NULL;
    InsertVertex(vertex);
    if (surface_index >= 0) {
      MPSurface *surface = surfaces.Kth(surface_index);
      surface->InsertVertex(vertex);
    }
  }

  // Read panoramas
  for (int i = 0; i < npanoramas; i++) {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%s", name_buffer);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &region_index);
    fscanf(fp, "%d", &dummy);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    for (int j = 0; j < 5; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "P")) { fprintf(stderr, "Error reading panorama %d\n", i); return 0; }
    MPPanorama *panorama = new MPPanorama();
    panorama->position = position;
    panorama->name = (strcmp(name_buffer, "-")) ? strdup(name_buffer) : NULL;
    InsertPanorama(panorama);
    if (region_index >= 0) {
      MPRegion *region = regions.Kth(region_index);
      region->InsertPanorama(panorama);
    }
  }

  // Read images
  for (int i = 0; i < nimages; i++) {
    double intrinsics[9];
    double extrinsics[16];
    int camera_index, yaw_index, width, height;
    char depth_filename[1024], color_filename[1024];
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &panorama_index);
    fscanf(fp, "%s", name_buffer);
    fscanf(fp, "%d", &camera_index);
    fscanf(fp, "%d", &yaw_index);
    for (int j = 0; j < 16; j++) fscanf(fp, "%lf", &extrinsics[j]);
    for (int j = 0; j < 9; j++) fscanf(fp, "%lf", &intrinsics[j]);
    fscanf(fp, "%d%d", &width, &height);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    for (int j = 0; j < 5; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "I")) { fprintf(stderr, "Error reading image %d\n", i); return 0; }
    sprintf(depth_filename, "%s_d%d_%d.png", name_buffer, camera_index, yaw_index);
    sprintf(color_filename, "%s_i%d_%d.jpg", name_buffer, camera_index, yaw_index);
    MPImage *image = new MPImage();
    image->name = (strcmp(name_buffer, "-")) ? strdup(name_buffer) : NULL;
    image->camera_index = camera_index;
    image->yaw_index = yaw_index;
    image->rgbd.SetNPixels(width, height);
    image->rgbd.SetExtrinsics(R4Matrix(extrinsics));
    image->rgbd.SetIntrinsics(R3Matrix(intrinsics));
    image->rgbd.SetDepthFilename(depth_filename);
    image->rgbd.SetColorFilename(color_filename);
    image->rgbd.SetName(name_buffer);
    image->extrinsics = R4Matrix(extrinsics);
    image->intrinsics = R3Matrix(intrinsics);
    image->width = width;
    image->height = height;
    image->position = position;
    InsertImage(image);
    if (panorama_index >= 0) {
      MPPanorama *panorama = panoramas.Kth(panorama_index);
      panorama->InsertImage(image);
    }
  }

  // Read categories
  for (int i = 0; i < ncategories; i++) {
    int label_id, mpcat40_id;
    char label_name[1024], mpcat40_name[1024];
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d %s", &label_id, label_name);
    fscanf(fp, "%d %s", &mpcat40_id, mpcat40_name);
    for (int j = 0; j < 5; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "C")) { fprintf(stderr, "Error reading category %d\n", i); return 0; }
    char *label_namep = label_name; while (*label_namep) { if (*label_namep == '#') *label_namep = ' '; label_namep++; }
    char *mpcat40_namep = mpcat40_name; while (*mpcat40_namep) { if (*mpcat40_namep == '#') *mpcat40_namep = ' '; mpcat40_namep++; }
    MPCategory *category = new MPCategory();
    category->label_id = label_id;
    category->mpcat40_id = mpcat40_id;
    if (strcmp(label_name, "-")) category->label_name = strdup(label_name);
    if (strcmp(mpcat40_name, "-")) category->mpcat40_name = strdup(mpcat40_name);
    InsertCategory(category);
  }
    
  // Read objects
  for (int i = 0; i < nobjects; i++) {
    R3Vector axis0, axis1, radius;
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &region_index);
    fscanf(fp, "%d", &category_index);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    fscanf(fp, "%lf%lf%lf", &axis0[0], &axis0[1], &axis0[2]);
    fscanf(fp, "%lf%lf%lf", &axis1[0], &axis1[1], &axis1[2]);
    fscanf(fp, "%lf%lf%lf", &radius[0], &radius[1], &radius[2]);
    for (int j = 0; j < 8; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "O")) { fprintf(stderr, "Error reading object %d\n", i); return 0; }
    MPObject *object = new MPObject();
    object->position = position;
    object->obb = R3OrientedBox(position, axis0, axis1, radius[0], radius[1], radius[2]);
    InsertObject(object);
    if (region_index >= 0) {
      MPRegion *region = regions.Kth(region_index);
      region->InsertObject(object);
    }
    if (category_index >= 0) {
      MPCategory *category = categories.Kth(category_index);
      category->InsertObject(object);
    }
  }
    
  // Read segments
  for (int i = 0; i < nsegments; i++) {
    fscanf(fp, "%s", cmd);
    fscanf(fp, "%d", &house_index);
    fscanf(fp, "%d", &object_index);
    fscanf(fp, "%d", &id);
    fscanf(fp, "%lf", &area);
    fscanf(fp, "%lf%lf%lf", &position[0], &position[1], &position[2]);
    fscanf(fp, "%lf%lf%lf%lf%lf%lf", &box[0][0], &box[0][1], &box[0][2], &box[1][0], &box[1][1], &box[1][2]);
    for (int j = 0; j < 5; j++) fscanf(fp, "%d", &dummy);
    if (strcmp(cmd, "E")) { fprintf(stderr, "Error reading segment %d\n", i); return 0; }
    MPSegment *segment = new MPSegment();
    segment->id = id;
    segment->area = area;
    segment->position = position;
    segment->bbox = box;
    InsertSegment(segment);
    if (object_index >= 0) {
      MPObject *object = objects.Kth(object_index);
      object->InsertSegment(segment);
    }
  }
    
  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int MPHouse::
WriteAsciiFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open house file %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "ASCII 1.1\n");
  fprintf(fp, "H  ");
  fprintf(fp, "%s  ", (name) ? name : "-");
  fprintf(fp, "%s  ", (label) ? label : "-");
  fprintf(fp, "%d ", images.NEntries());
  fprintf(fp, "%d ", panoramas.NEntries());
  fprintf(fp, "%d ", vertices.NEntries());
  fprintf(fp, "%d ", surfaces.NEntries());
  fprintf(fp, "%d ", segments.NEntries());
  fprintf(fp, "%d ", objects.NEntries());
  fprintf(fp, "%d ", categories.NEntries());
  fprintf(fp, "%d ", regions.NEntries());
  fprintf(fp, "%d  ", portals.NEntries());
  fprintf(fp, "%d  ", levels.NEntries());
  for (int i = 0; i < 5; i++) fprintf(fp, "0 ");
  fprintf(fp, "%g %g %g  %g %g %g  ", bbox[0][0], bbox[0][1], bbox[0][2], bbox[1][0], bbox[1][1], bbox[1][2]);
  for (int i = 0; i < 5; i++) fprintf(fp, "0 ");
  fprintf(fp, "\n");

  // Write levels
  for (int i = 0; i < levels.NEntries(); i++) {
    MPLevel *level = levels.Kth(i);
    fprintf(fp, "L  ");
    fprintf(fp, "%d  ", level->house_index);
    fprintf(fp, "0  ");
    fprintf(fp, "%s  ", (level->label) ? (level->label) : "-");
    fprintf(fp, "%g %g %g  ", level->position.X(), level->position.Y(), level->position.Z());
    fprintf(fp, "%g %g %g  %g %g %g  ", level->bbox[0][0], level->bbox[0][1], level->bbox[0][2], level->bbox[1][0], level->bbox[1][1], level->bbox[1][2]);
    for (int j = 0; j < 5; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Write regions
  for (int i = 0; i < regions.NEntries(); i++) {
    MPRegion *region = regions.Kth(i);
    fprintf(fp, "R  ");
    fprintf(fp, "%d ", region->house_index);
    fprintf(fp, "%d  ", (region->level) ? region->level->house_index : -1);
    fprintf(fp, "0 ");
    fprintf(fp, "0  ");
    fprintf(fp, "%s  ", (region->label) ? region->label : "-");
    fprintf(fp, "%g %g %g  ", region->position.X(), region->position.Y(), region->position.Z());
    fprintf(fp, "%g %g %g  %g %g %g  ", region->bbox[0][0], region->bbox[0][1], region->bbox[0][2], region->bbox[1][0], region->bbox[1][1], region->bbox[1][2]);
    fprintf(fp, "%g  ", region->height);
    for (int j = 0; j < 4; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Write portals
  for (int i = 0; i < portals.NEntries(); i++) {
    MPPortal *portal = portals.Kth(i);
    fprintf(fp, "R  ");
    fprintf(fp, "%d ", portal->house_index);
    fprintf(fp, "%d  ", (portal->regions[0]) ? portal->regions[0]->house_index : -1);
    fprintf(fp, "%d  ", (portal->regions[1]) ? portal->regions[1]->house_index : -1);
    fprintf(fp, "%s  ", (portal->label) ? portal->label : "-");
    fprintf(fp, "%g %g %g  ", portal->span[0][0], portal->span[0][1], portal->span[0][2]);
    fprintf(fp, "%g %g %g  ", portal->span[1][0], portal->span[1][1], portal->span[1][2]);
    for (int j = 0; j < 4; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Write surfaces
  for (int i = 0; i < surfaces.NEntries(); i++) {
    MPSurface *surface = surfaces.Kth(i);
    fprintf(fp, "S  ");
    fprintf(fp, "%d %d  ", surface->house_index, (surface->region) ? surface->region->house_index : -1);
    fprintf(fp, "0  ");
    fprintf(fp, "%s  ", (surface->label) ? surface->label : "-");
    fprintf(fp, "%g %g %g  ", surface->position.X(), surface->position.Y(), surface->position.Z());
    fprintf(fp, "%g %g %g  ", surface->normal.X(), surface->normal.Y(), surface->normal.Z());
    fprintf(fp, "%g %g %g  %g %g %g  ", surface->bbox[0][0], surface->bbox[0][1], surface->bbox[0][2], surface->bbox[1][0], surface->bbox[1][1], surface->bbox[1][2]);
    for (int j = 0; j < 5; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Write vertices
  for (int i = 0; i < surfaces.NEntries(); i++) {
    MPSurface *surface = surfaces.Kth(i);
    for (int j = 0; j < surface->vertices.NEntries(); j++) {
      MPVertex *vertex = surface->vertices.Kth(j);
      fprintf(fp, "V  ");
      fprintf(fp, "%d %d  ", vertex->house_index, (vertex->surface) ? vertex->surface->house_index : -1);
      fprintf(fp, "%s  ", (vertex->label) ? vertex->label : "-");
      fprintf(fp, "%g %g %g  ", vertex->position.X(), vertex->position.Y(), vertex->position.Z());
      fprintf(fp, "%g %g %g  ", vertex->normal.X(), vertex->normal.Y(), vertex->normal.Z());
      for (int j = 0; j < 3; j++) fprintf(fp, "0 ");
      fprintf(fp, "\n");
    }
  }

  // Write panoramas
  for (int i = 0; i < panoramas.NEntries(); i++) {
    MPPanorama *panorama = panoramas.Kth(i);
    fprintf(fp, "P  ");
    fprintf(fp, "%s  ", (panorama->name) ? panorama->name : "-");
    fprintf(fp, "%d %d  ", panorama->house_index, (panorama->region) ? panorama->region->house_index : -1);
    fprintf(fp, "0  ");
    fprintf(fp, "%g %g %g  ", panorama->position.X(), panorama->position.Y(), panorama->position.Z());
    for (int j = 0; j < 5; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }

  // Write images
  for (int i = 0; i < images.NEntries(); i++) {
    MPImage *image = images.Kth(i);
    fprintf(fp, "I  ");
    fprintf(fp, "%d ", image->house_index);
    fprintf(fp, "%d   ", (image->panorama) ? image->panorama->house_index : -1);
    fprintf(fp, "%s  ", (image->name) ? image->name : "-");
    fprintf(fp, "%d ", image->camera_index);
    fprintf(fp, "%d  ", image->yaw_index);
    for (int j = 0; j < 16; j++) fprintf(fp, "%g ", image->extrinsics[j/4][j%4]);
    fprintf(fp, " ");
    for (int j = 0; j < 9; j++) fprintf(fp, "%g ", image->intrinsics[j/3][j%3]);
    fprintf(fp, " ");
    fprintf(fp, "%d %d  ", image->width, image->height);
    fprintf(fp, "%g %g %g  ", image->position.X(), image->position.Y(), image->position.Z());
    for (int j = 0; j < 5; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }

  // Write categories
  for (int i = 0; i < categories.NEntries(); i++) {
    MPCategory *category = categories.Kth(i);
    char label_name[1024], mpcat40_name[1024];
    strncpy(label_name, (category->label_name) ? category->label_name : "-", 1024);
    strncpy(mpcat40_name, (category->mpcat40_name) ? category->mpcat40_name : "-", 1024);
    char *label_namep = label_name; while (*label_namep) { if (*label_namep == ' ') *label_namep = '#'; label_namep++; }
    char *mpcat40_namep = mpcat40_name; while (*mpcat40_namep) { if (*mpcat40_namep == ' ') *mpcat40_namep = '#'; mpcat40_namep++; }
    fprintf(fp, "C  ");
    fprintf(fp, "%d  ", category->house_index);
    fprintf(fp, "%d %s  ", category->label_id, label_name);
    fprintf(fp, "%d %s  ", category->mpcat40_id, mpcat40_name);
    for (int j = 0; j < 5; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Write objects
  for (int i = 0; i < objects.NEntries(); i++) {
    MPObject *object = objects.Kth(i);
    fprintf(fp, "O  ");
    fprintf(fp, "%d ", object->house_index);
    fprintf(fp, "%d ", (object->region) ? object->region->house_index : -1);
    fprintf(fp, "%d  ", (object->category) ? object->category->house_index : -1);
    fprintf(fp, "%g %g %g  ", object->obb.Center().X(), object->obb.Center().Y(), object->obb.Center().Z());
    fprintf(fp, "%g %g %g  ", object->obb.Axis(0).X(), object->obb.Axis(0).Y(), object->obb.Axis(0).Z());
    fprintf(fp, "%g %g %g  ", object->obb.Axis(1).X(), object->obb.Axis(1).Y(), object->obb.Axis(1).Z());
    fprintf(fp, "%g %g %g  ", object->obb.Radius(0), object->obb.Radius(1), object->obb.Radius(2));
    for (int j = 0; j < 8; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Write segments
  for (int i = 0; i < segments.NEntries(); i++) {
    MPSegment *segment = segments.Kth(i);
    fprintf(fp, "E  ");
    fprintf(fp, "%d ", segment->house_index);
    fprintf(fp, "%d ", (segment->object) ? segment->object->house_index : -1);
    fprintf(fp, "%d  ", segment->id);
    fprintf(fp, "%g  ", segment->area);
    fprintf(fp, "%g %g %g  ", segment->position.X(), segment->position.Y(), segment->position.Z());
    fprintf(fp, "%g %g %g  %g %g %g  ", segment->bbox[0][0], segment->bbox[0][1], segment->bbox[0][2], segment->bbox[1][0], segment->bbox[1][1], segment->bbox[1][2]);
    for (int j = 0; j < 5; j++) fprintf(fp, "0 ");
    fprintf(fp, "\n");
  }
    
  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Other file parsing
////////////////////////////////////////////////////////////////////////

int MPHouse::
ReadMeshFile(const char *filename)
{
  // Allocate mesh
  mesh = new R3Mesh();
  if (!mesh) {
    fprintf(stderr, "Unable to allocate mesh for %s\n", filename);
    return 0;
  }

  // Read mesh from file
  if (!mesh->ReadFile(filename)) {
    delete mesh;
    return 0;
  }

  // Update segments
  if (segments.NEntries() > 0) {
    RNMap<int, MPSegment *> id_to_segment;
    for (int i = 0; i < segments.NEntries(); i++) {
      MPSegment *segment = segments.Kth(i);
      id_to_segment.Insert(segment->id, segment);
    }
    for (int i = 0; i < mesh->NFaces(); i++) {
      R3MeshFace *face = mesh->Face(i);
      int segment_id = mesh->FaceMaterial(face);
      if (segment_id < 0) continue;
      MPSegment *segment = NULL;
      if (id_to_segment.Find(segment_id, &segment)) {
        segment->bbox.Union(mesh->FaceBBox(face));
        segment->faces.Insert(face);
        segment->mesh = mesh;
      }
    }
  }

  // Update bbox
  bbox.Union(mesh->BBox());

  // Return success
  return 1;
}



int MPHouse::
ReadSceneFile(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate scene
  scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene for %s\n", filename);
    return 0;
  }

  // Read scene from file
  if (!scene->ReadFile(filename)) {
    delete scene;
    return 0;
  }

  // Process scene
  scene->RemoveReferences();
  scene->RemoveTransformations();
  scene->RemoveHierarchy();

  // Update bbox
  bbox.Union(scene->BBox());

  // Return success
  return 1;
}



int MPHouse::
ReadCategoryFile(const char *filename)
{
  // Open file
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open category file %s\n", filename);
    return 0;
  }

  // Read keys from first line
  RNArray<char *> keys;
  char key_buffer[4096];
  if (fgets(key_buffer, 4096, fp)) {
    char *token = key_buffer;
    while (*token &&  (*token != '\r') && (*token != '\n')) {
      while (*token == ' ') token++;
      keys.Insert(token);
      while (*token && (*token != '\t') && (*token != '\r') && (*token != '\n')) token++;
      *(token++) = '\0';
    }
  }
  
  // Extract indices of interesting info
  int label_id_k = -1;
  int label_name_k = -1;
  int mpcat40_id_k = -1;
  int mpcat40_name_k = -1;
  for (int i = 0; i < keys.NEntries(); i++) {
    if (!strcmp(keys[i], "index")) label_id_k = i;
    else if (!strcmp(keys[i], "raw_category")) label_name_k = i; 
    else if (!strcmp(keys[i], "mpcat40index")) mpcat40_id_k = i; 
    else if (!strcmp(keys[i], "mpcat40")) mpcat40_name_k = i; 
  }

  // Check if found key fields in header
  if ((label_id_k < 0) || (label_name_k < 0) || (mpcat40_id_k < 0)) {
    fprintf(stderr, "Did not find index, raw_category, and mpcat40index in header of %s\n", filename);
    return 0;
  }

  // Read subsequent lines of file
  char value_buffer[4096];
  while (fgets(value_buffer, 4096, fp)) {
    // Read values
    RNArray<char *> values;
    char *token = value_buffer;
    while (*token &&  (*token != '\r') && (*token != '\n')) {
      while (*token == ' ') token++;
      values.Insert(token);
      while (*token && (*token != '\t') && (*token != '\r') && (*token != '\n')) token++;
      *(token++) = '\0';
    }

    // Create category
    MPCategory *category = new MPCategory();
    category->label_id = atoi(values[label_id_k]);
    category->label_name = strdup(values[label_name_k]);
    category->mpcat40_id = atoi(values[mpcat40_id_k]);
    category->mpcat40_name = strdup(values[mpcat40_name_k]);

    // Insert category
    InsertCategory(category);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int MPHouse::
ReadSegmentFile(const char *filename)
{
  // Check mesh
  if (!mesh) return 0;

  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE* fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open json file %s\n", filename);
    return 0;
  }

  // Read file 
  std::string text;
  fseek(fp, 0, SEEK_END);
  long const size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char* buffer = new char[size + 1];
  unsigned long const usize = static_cast<unsigned long const>(size);
  if (fread(buffer, 1, usize, fp) != usize) { fprintf(stderr, "Unable to read %s\n", filename); return 0; }
  else { buffer[size] = 0; text = buffer; }
  delete[] buffer;

  // Close file
  fclose(fp);

  // Parse file
  Json::Value json_root;
  Json::Reader json_reader;
  if (!json_reader.parse(text, json_root, false)) {
    fprintf(stderr, "Unable to parse %s\n", filename);
    return 0;
  }

  // Check file
  if (!json_root.isMember("segIndices")) {
    fprintf(stderr, "Segment file %s has no segIndices\n", filename);
    return 0;
  }

  // Parse segment identifiers
  Json::Value json_segments = json_root["segIndices"];
  RNMap<int, MPSegment *> id_to_segment;
  for (Json::ArrayIndex i = 0; i < json_segments.size(); i++) {
    Json::Value json_segment = json_segments[i];
    int face_index = (int) i;
    if ((face_index >= 0) && (face_index < mesh->NFaces())) {
      MPSegment *segment = NULL;
      int segment_id = atoi(json_segment.asString().c_str());
      if (!id_to_segment.Find(segment_id, &segment)) {
        segment = new MPSegment();
        segment->id = segment_id;
        segment->mesh = mesh;
        InsertSegment(segment);
        id_to_segment.Insert(segment_id, segment);
      }
      R3MeshFace *face = mesh->Face(face_index);
      mesh->SetFaceMaterial(face, segment_id);
      segment->faces.Insert(face);
      segment->area += mesh->FaceArea(face);
      segment->bbox.Union(mesh->FaceBBox(face));
      segment->position = segment->bbox.Centroid();
    }
  }
  
  // Return success
  return 1;
}



int MPHouse::
ReadObjectFile(const char *filename)
{
  // Open file
  FILE* fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open json file %s\n", filename);
    return 0;
  }

  // Read file 
  std::string text;
  fseek(fp, 0, SEEK_END);
  long const size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char* buffer = new char[size + 1];
  unsigned long const usize = static_cast<unsigned long const>(size);
  if (fread(buffer, 1, usize, fp) != usize) { fprintf(stderr, "Unable to read %s\n", filename); return 0; }
  else { buffer[size] = 0; text = buffer; }
  delete[] buffer;

  // Close file
  fclose(fp);

  // Create maps to find categories
  RNMap<int, MPCategory *> id_to_category;
  RNMap<std::string, MPCategory *> name_to_category;
  for (int i = 0; i < categories.NEntries(); i++) {
    MPCategory *category = categories.Kth(i);
    id_to_category.Insert(category->label_id, category);
    name_to_category.Insert(category->label_name, category);
  }
  
  // Create a map to find segments
  RNMap<int, MPSegment *> id_to_segment;
  for (int i = 0; i < segments.NEntries(); i++) {
    MPSegment *segment = segments.Kth(i);
    id_to_segment.Insert(segment->id, segment);
  }
  
  // Parse file
  Json::Value json_root;
  Json::Reader json_reader;
  if (!json_reader.parse(text, json_root, false)) {
    fprintf(stderr, "Unable to parse %s\n", filename);
    return 0;
  }

  // Check file
  if (!json_root.isMember("segGroups")) {
    fprintf(stderr, "Annotation file %s has no groups\n", filename);
    return 0;
  }

  // Parse objects
  Json::Value json_groups = json_root["segGroups"];
  for (Json::ArrayIndex object_index = 0; object_index < json_groups.size(); object_index++) {
    Json::Value json_group = json_groups[object_index];

    // Create object
    MPObject *object = new MPObject();

    // Associate object with category
    if (json_group.isMember("label")) {
      Json::Value json_label = json_group["label"];
      std::string label_name = json_label.asString();
      MPCategory *category = NULL;
      if (name_to_category.Find(label_name, &category)) {
        category->InsertObject(object);
      }
    }
    else if (json_group.isMember("label_index")) {
      Json::Value json_label_index = json_group["label_index"];
      int label_id = atoi(json_label_index.asString().c_str());
      MPCategory *category = NULL;
      if (id_to_category.Find(label_id, &category)) {
        category->InsertObject(object);
      }
    }

    // Associate segments with object
    if (json_group.isMember("segments")) {
      Json::Value json_segments = json_group["segments"];
      for (Json::ArrayIndex i = 0; i < json_segments.size(); i++) {
        Json::Value json_segment = json_segments[i];
        int segment_id = atoi(json_segment.asString().c_str());
        MPSegment *segment = NULL;
        if (id_to_segment.Find(segment_id, &segment)) {
          object->InsertSegment(segment);
        }
      }
    }
  
    // Parse obb
    if (json_group.isMember("obb")) {
      Json::Value json_obb = json_group["obb"];

      // Get obb centroid
      R3Point centroid(0, 0, 0);
      if (json_obb.isMember("centroid")) {
        Json::Value json_items = json_obb["centroid"];
        if (json_items.size() >= 3) {
          for (int i = 0; i < 3; i++) {
            Json::Value json_item = json_items[i];
            centroid[i] = json_item.asDouble();
          }
        }
      }

      // Get obb axes lengths
      R3Vector axesLengths(0, 0, 0);
      if (json_obb.isMember("axesLengths")) {
        Json::Value json_items = json_obb["axesLengths"];
        if (json_items.size() >= 3) {
          for (int i = 0; i < 3; i++) {
            Json::Value json_item = json_items[i];
            axesLengths[i] = json_item.asDouble();
          }
        }
      }

      // Get obb axes
      R3Vector normalizedAxes[3] = { R3Vector(0,0,0), R3Vector(0,0,0), R3Vector(0,0,0) };
      if (json_obb.isMember("normalizedAxes")) {
        Json::Value json_items = json_obb["normalizedAxes"];
        if (json_items.size() >= 9) {
          for (int i = 0; i < 9; i++) {
            Json::Value json_item = json_items[i];
            normalizedAxes[i/3][i%3] = json_item.asDouble();
          }
        }
      }

      // Set obb
      object->position = centroid;
      object->obb.Reset(centroid, normalizedAxes[0], normalizedAxes[1], axesLengths[0], axesLengths[1], axesLengths[2]);
    }

    // Associate object with region
    MPRegion *region = FindRegion(object->obb.Centroid());
    if (region) region->InsertObject(object);

    // Insert object into house
    InsertObject(object);
  }

  // Update stuff
  if (mesh) {
    // Update mesh
    for (int i = 0; i < objects.NEntries(); i++) {
      MPObject *object = objects.Kth(i);
      MPCategory *category = object->category;
      int object_index = i;
      int category_index = (category) ? category->label_id : -1;
      for (int j = 0; j < object->segments.NEntries(); j++) {
        MPSegment *segment = object->segments.Kth(j);
        for (int k = 0; k < segment->faces.NEntries(); k++) {
          R3MeshFace *face = segment->faces.Kth(k);
          mesh->SetFaceMaterial(face, segment->id);
          mesh->SetFaceSegment(face, object_index);
          mesh->SetFaceCategory(face, category_index);
        }
      }
    }

    // Compute object oriented bounding boxes
    for (int i = 0; i < objects.NEntries(); i++) {
      MPObject *object = objects.Kth(i);
      if (object->obb.MaxRadius() <= 0) {
        RNArray<R3Point *> points;
        for (int j = 0; j < object->segments.NEntries(); j++) {
          MPSegment *segment = object->segments.Kth(j);
          for (int k = 0; k < segment->faces.NEntries(); k++) {
            R3MeshFace *face = segment->faces.Kth(k);
            for (int m = 0; m < 3; m++) {
              R3MeshVertex *vertex = mesh->VertexOnFace(face, m);
              points.Insert((R3Point *) &(mesh->VertexPosition(vertex)));
            }
          }
        }
        object->obb = R3OrientedBox(points);
      }
    }
  }

  // Return success
  return 1;
}



int MPHouse::
ReadConfigurationFile(const char *filename) 
{
  // Read file
  RGBDConfiguration configuration;
  if (!configuration.ReadFile(filename)) {
    fprintf(stderr, "Unable to read configuration from %s\n", filename);
    return 0;
  }

  // Create images
  for (int i = 0; i < configuration.NImages(); i++) {
    RGBDImage *image = configuration.Image(i);
    R3Point viewpoint = image->WorldViewpoint();
    int width = (image->NPixels(RN_X) > 0) ? image->NPixels(RN_X) : 1280;
    int height = (image->NPixels(RN_Y) > 0) ? image->NPixels(RN_Y) : 1024;
    
    // Get image name
    char image_name_buffer[4096];
    strncpy(image_name_buffer, image->DepthFilename(), 4096);
    char *name = (strrchr(image_name_buffer, '/')) ? strrchr(image_name_buffer, '/')+1 : image_name_buffer;
    char *camera_index = strstr(name, "_d");
    if (camera_index) { *camera_index = '\0'; camera_index++; *camera_index = '\0'; camera_index++; }
    char *yaw_index = strchr(camera_index, '_');
    if (yaw_index) { *yaw_index = '\0'; yaw_index++; }
    char *ext_index = strchr(yaw_index, '.');
    if (ext_index) { *ext_index = '\0'; ext_index++; }
    
    // Create image
    MPImage *img = new MPImage();
    img->name = (name && strcmp(name, "-")) ? strdup(name) : NULL;
    img->camera_index = (camera_index) ? atoi(camera_index) : -1;
    img->yaw_index = (yaw_index) ? atoi(yaw_index) : -1;
    img->extrinsics = image->CameraToWorld().Inverse().Matrix();
    img->intrinsics = image->Intrinsics();
    img->width = width;
    img->height = height;
    img->position = image->WorldViewpoint();
    InsertImage(img);

    // Fill in rgbd info
    char depth_filename[1024], color_filename[1024];
    sprintf(depth_filename, "%s_d%d_%d.png", img->name, img->camera_index, img->yaw_index);
    sprintf(color_filename, "%s_i%d_%d.jpg", img->name, img->camera_index, img->yaw_index);
    img->rgbd.SetNPixels(img->width, img->height);
    img->rgbd.SetExtrinsics(R4Matrix(img->extrinsics));
    img->rgbd.SetIntrinsics(R3Matrix(img->intrinsics));
    img->rgbd.SetDepthFilename(depth_filename);
    img->rgbd.SetColorFilename(color_filename);
    img->rgbd.SetName(img->name);

    // Find panorama
    MPPanorama *panorama = NULL;
    for (int j = 0; j < panoramas.NEntries(); j++) {
      MPPanorama *p = panoramas.Kth(j);
      if (p->name && name && (!strcmp(p->name, name))) {
        if (R3SquaredDistance(p->position, viewpoint) < 0.5) {
          panorama = p;
          break;
        }
      }
    }

    // Insert image into panorama
    if (panorama) {
      // Just insert image
      panorama->InsertImage(img);
    }
    else {
      // Create panorama and insert image
      panorama = new MPPanorama();
      panorama->name = (name && strcmp(name, "-")) ? strdup(name) : NULL;
      panorama->position = viewpoint;
      panorama->InsertImage(img);
      InsertPanorama(panorama);

      // Find region for panorama
      MPRegion *region = FindRegion(viewpoint);
      if (region) region->InsertPanorama(panorama);
    }
  }
  
  // Update bbox
  bbox.Union(configuration.WorldBBox());

  // Return success
  return 1;
}



