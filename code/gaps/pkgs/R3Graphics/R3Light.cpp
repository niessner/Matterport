/* Source file for the R3 light class */



/* Include files */

#include "R3Graphics.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Light);



/* Public variables */

RNScalar R3ambient_light_intensity = 0.2;
RNRgb R3ambient_light_color(1.0, 1.0, 1.0);



/* Public functions */

int 
R3InitLight()
{
    /* Return success */
    return TRUE;
}



void 
R3StopLight()
{
}



R3Light::
R3Light(const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      active(1),
      intensity(1),
      color(0.5, 0.5, 0.5),
      id(0)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Light::
R3Light(const R3Light& light, const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      active(light.active),
      intensity(light.intensity),
      color(light.color),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Light::
R3Light(const RNRgb& color,
	RNScalar intensity,
	RNBoolean active,
        const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      active(active),
      intensity(intensity),
      color(color),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Light::
~R3Light(void)
{
    // Unload
    if (id > 0) { fprintf(stderr, "Will not get here because id is not used yet.\n"); }

    // Remove light from scene
    if (scene) scene->RemoveLight(this);

    // Free name
    if (name) free(name);
}



R3Light *R3Light::
Copy(void) const
{
  // Return copy of appropriate type
  if (ClassID() == R3DirectionalLight::CLASS_ID()) return new R3DirectionalLight(*((R3DirectionalLight *) this));
  else if (ClassID() == R3PointLight::CLASS_ID()) return new R3PointLight(*((R3PointLight *) this));
  else if (ClassID() == R3SpotLight::CLASS_ID()) return new R3SpotLight(*((R3SpotLight *) this));
  else if (ClassID() == R3AreaLight::CLASS_ID()) return new R3AreaLight(*((R3AreaLight *) this));

  // Should never get here
  RNAbort("Unrecognized light class id: %d\n", ClassID());
  return NULL;
}



void R3Light::
SetName(const char *name)
{
    // Set name
    if (this->name) free(this->name);
    if (name) this->name = strdup(name);
    else this->name = NULL;
}



void R3Light::
SetActive(RNBoolean active)
{
    // Set active
    this->active = active;
}



void R3Light::
SetIntensity(RNScalar intensity)
{
    // Set intensity
    this->intensity = intensity;
}



void R3Light::
SetColor(const RNRgb& color)
{
    // Set color
    this->color = color;
}



void R3Light::
Transform(const R3Transformation& transformation)
{
  // Transform light
}



RNRgb R3Light::
IrradianceAtPoint(const R3Point& point) const
{
    // Return irradiance at point
    if (!IsActive()) return RNblack_rgb;
    return IntensityAtPoint(point) * color;
}



