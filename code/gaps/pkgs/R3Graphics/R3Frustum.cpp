// Source file for view frustum class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"



////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////

R3Frustum::
R3Frustum(void)
{
}



R3Frustum::
R3Frustum(const R3Camera& camera)
{
  // Set camera
  SetCamera(camera);
}



R3Frustum::
R3Frustum(const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up, 
  RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist)
{
  // Set camera
  SetCamera(R3Camera(viewpoint, towards, up, xfov, yfov, neardist, fardist));
}



void R3Frustum::
SetCamera(const R3Camera& camera)
{
  // Set camera
  this->camera = camera;

  // Compute halfspaces
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      halfspaces[dir][dim] = camera.Halfspace(dir, dim);
    }
  }
}



int R3Frustum::
Intersects(const R3Point& point) const
{
  // Check each halfspace
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      if (!R3Intersects(halfspaces[dir][dim], point)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



int R3Frustum::
Intersects(const R3Span& span) const
{
  // Check each halfspace
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      if (!R3Intersects(halfspaces[dir][dim], span)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



int R3Frustum::
Intersects(const R3Box& box) const
{
  // Check each halfspace
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      if (!R3Intersects(halfspaces[dir][dim], box)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



int R3Frustum::
Contains(const R3Point& point) const
{
  // Check each halfspace
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      if (!R3Contains(halfspaces[dir][dim], point)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



int R3Frustum::
Contains(const R3Span& span) const
{
  // Check each halfspace
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      if (!R3Contains(halfspaces[dir][dim], span)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



int R3Frustum::
Contains(const R3Box& box) const
{
  // Check each halfspace
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      if (!R3Contains(halfspaces[dir][dim], box)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



void R3Frustum::
Draw(void) const
{
  // Draw camera
  camera.Draw();
}








