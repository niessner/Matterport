/* Include file for RN graphics module */

#ifndef __RN__GRFX__H__
#define __RN__GRFX__H__



/* Turn off warnings about deprecated functions (for GLU and GLUT) */

#if (RN_OS == RN_MAC) 
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif



/* Graphics library include files */

#if (RN_2D_GRFX == RN_XLIB)
#   ifndef __RN_XLIB__
#       include "X11/Xlib.h"
#       include "X11/X.h"
#       define __RN_XLIB__
#   endif
#endif

#if (RN_3D_GRFX == RN_IRISGL)
#   ifndef __RN_IRISGL__
#       include <gl/gl.h>
#       define __RN_IRISGL__
#   endif
#   ifndef __RN_SPHEREGL__
#       include <gl/sphere.h>
#       define __RN_SPHEREGL__
#   endif
#elif (RN_3D_GRFX == RN_OPENGL)
#   ifndef __RN_OPENGL__
#       if (RN_OS == RN_MAC) 
#           include <OpenGL/gl.h>
#           include <OpenGL/glu.h>
#       else
#           include <GL/gl.h>
#           include <GL/glu.h>
#       endif
#       define __RN_OPENGL__
#   endif
#elif (RN_3D_GRFX == RN_3DR)
#   ifndef __RN_3DR__
#       include "3dr.h"
#       include "3dg.h"
        extern R3dHandle_t R3dr_rc;
        extern G3dHandle_t R3dr_gc;
#       define __RN_3DR__
#   endif
#endif



/* Compatibility definitions -- not sure about this */

#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif



/* Initialization functions */

int RNInitGrfx(void);
void RNStopGrfx(void);



/* 2D viewing functions */

void R2SetViewport(int xmin, int ymin, int xmax, int ymax);
void R2SetWindow(float xmin, float ymin, float xmax, float ymax);
void R2ScaleWindow(float xscale, float yscale);
void R2TranslateWindow(float xtranslate, float ytranslate);
void R2WindowToViewport(float wx, float wy, int *vx, int *vy);
void R2WindowToViewport(double wx, double wy, int *vx, int *vy);
void R2WindowToViewport(float wx, float wy, short *vx, short *vy);
void R2WindowToViewport(double wx, double wy, short *vx, short *vy);
void R2ViewportToWindow(int vx, int vy, float *wx, float *wy);
void R2ViewportToWindow(int vx, int vy, double *wx, double *wy);



/* 2D Polygon/Line drawing functions */

void R2BeginPolygon(void);
void R2EndPolygon(void);
void R2BeginLine(void);
void R2EndLine(void);
void R2BeginLoop(void);
void R2EndLoop(void);

void R2LoadPoint(float x, float y);
void R2LoadPoint(double x, double y);
void R2LoadPoint(const float point[2]);
void R2LoadPoint(const double point[2]);



/* 2D image drawing functions */

void R2DrawImage(int x, int y, int width, int height, int depth, void *data);



/* 2D primitive drawing functions */

void R2DrawText(float x, float y, const char *str);
void R2DrawText(double x, double y, const char *str);



/* 3D Polygon/Line drawing functions */

void R3BeginPolygon(void);
void R3EndPolygon(void);
void R3BeginLine(void);
void R3EndLine(void);
void R3BeginLoop(void);
void R3EndLoop(void);

void R3LoadPoint(float x, float y, float z);
void R3LoadPoint(double x, double y, float z);
void R3LoadPoint(const float point[3]);
void R3LoadPoint(const double point[3]);

void R3LoadNormal(float x, float y, float z);
void R3LoadNormal(double x, double y, float z);
void R3LoadNormal(const float normal[3]);
void R3LoadNormal(const double normal[3]);

void R3LoadTextureCoords(float x, float y);
void R3LoadTextureCoords(double x, double y);
void R3LoadTextureCoords(const float texcoords[2]);
void R3LoadTextureCoords(const double texcoords[2]);
void R3LoadTextureCoords(float x, float y, float z);
void R3LoadTextureCoords(double x, double y, double z);



/* 3D primitive drawing functions */

void R3DrawText(float x, float y, float z, const char *str);
void R3DrawText(double x, double y, double z, const char *str);



/* Color drawing functions */

void RNLoadRgb(float r, float g, float b);
void RNLoadRgb(double r, double g, double b);
void RNLoadRgb(const float rgb[3]);
void RNLoadRgb(const double rgb[3]);



/* Error handling functions */

void RNGrfxError(const char *message);



/* Graphics library specific functions */

#if (RN_2D_GRFX == RN_XLIB)
    int RNInitGrfx(Display *display, Window window, GC gc);
#endif



/* Inline functions */

#include "RNGrfx.I"



#endif







