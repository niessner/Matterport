/* Include file for R3 graphics module */

#ifndef __R3__GRAPHICS__H__
#define __R3__GRAPHICS__H__



/* Class declarations */

class R3Scene;
class R3SceneNode;
class R3SceneElement;



/* Dependency include files */

#include "R3Shapes/R3Shapes.h"



/* Material include files */

#include "R3Graphics/R3Brdf.h"
#include "R3Graphics/R2Texture.h"
#include "R3Graphics/R3Material.h"



/* Light include files */

#include "R3Graphics/R3Light.h"
#include "R3Graphics/R3DirectionalLight.h"
#include "R3Graphics/R3PointLight.h"
#include "R3Graphics/R3SpotLight.h"
#include "R3Graphics/R3AreaLight.h"



/* Viewing include files */

#include "R3Graphics/R2Viewport.h"
#include "R3Graphics/R3Camera.h"
#include "R3Graphics/R3Frustum.h"
#include "R3Graphics/R3Viewer.h"



/* Scene include files */

#include "R3Graphics/R3SceneReference.h"
#include "R3Graphics/R3SceneElement.h"
#include "R3Graphics/R3SceneNode.h"
#include "R3Graphics/R3Scene.h"



/* Initialization functions */

int R3InitGraphics(void);
void R3StopGraphics(void);



#endif











