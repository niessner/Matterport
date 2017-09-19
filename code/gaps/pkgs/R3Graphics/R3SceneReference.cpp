/* Source file for the R3 scene reference class */



/* Include files */

#include "R3Graphics.h"



/* Member functions */

R3SceneReference::
R3SceneReference(R3Scene *referenced_scene)
  : referenced_scene(referenced_scene),
    materials(),
    info(),
    name(NULL),
    data(NULL)
{
}



R3SceneReference::
R3SceneReference(R3Scene *referenced_scene, const RNArray<R3Material *>& materials)
  : referenced_scene(referenced_scene),
    materials(),
    info(),
    name(NULL),
    data(NULL)
{
  // Copy pointers to materials
  this->materials = materials;
}



R3SceneReference::
~R3SceneReference(void)
{
  // Delete name
  if (name) free(name);
}



const R3Box& R3SceneReference::
BBox(void) const
{
  // Return bounding box
  if (!referenced_scene) return R3null_box;
  return referenced_scene->BBox();
}                        



const R3Point R3SceneReference::
Centroid(void) const
{
  // Return centroid
  if (!referenced_scene) return R3zero_point;
  return referenced_scene->Centroid();
}                        



const RNInterval R3SceneReference::
NFacets(void) const
{
  // Return nfacets
  if (!referenced_scene) return RNzero_interval;
  return referenced_scene->NFacets();
}                        



const RNLength R3SceneReference::
Length(void) const
{
  // Return length
  if (!referenced_scene) return 0;
  return referenced_scene->Length();
}                        



const RNArea R3SceneReference::
Area(void) const
{
  // Return area
  if (!referenced_scene) return 0;
  return referenced_scene->Area();
}                        



const RNVolume R3SceneReference::
Volume(void) const
{
  // Return volume
  if (!referenced_scene) return 0;
  return referenced_scene->Volume();
}                        



const R3Point R3SceneReference::
ClosestPoint(const R3Point& point) const
{
  // Return closest point
  if (!referenced_scene) return R3zero_point;
  return referenced_scene->ClosestPoint(point);
}                        



void R3SceneReference::
InsertInfo(const char *key, const char *value) 
{
  // Insert key-value pair
  info.Insert(key, value);
}



void R3SceneReference::
ReplaceInfo(const char *key, const char *value) 
{
  // Replace key-value pair
  info.Replace(key, value);
}



void R3SceneReference::
RemoveInfo(const char *key) 
{
  // Insert key-value pair
  info.Remove(key);
}



void R3SceneReference::
Draw(const R3DrawFlags draw_flags, const RNArray<R3Material *> *mats) const
{
  // Check scene
  if (!referenced_scene) return;
  
  // Draw scene with materials
  referenced_scene->Root()->Draw(draw_flags, &materials);
}



