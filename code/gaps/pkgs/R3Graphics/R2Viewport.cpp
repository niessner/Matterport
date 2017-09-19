/* Source file for the viewport class */



/* Include files */

#include "R3Graphics.h"



/* Public variables */

R2Viewport R2default_viewport(0, 0, 128, 128);



int R3InitViewport() 
{
    /* Return OK status */
    return TRUE;
}



void R3StopViewport()
{
}



R2Viewport::
R2Viewport(void)
{
}



R2Viewport::
R2Viewport (int xmin, int ymin, int width, int height)
    : xmin(xmin), 
      ymin(ymin),
      width(width),
      height(height)
{
}



void R2Viewport::
Move(int xmin, int ymin)
{
    // Set viewport variables
    this->xmin = xmin;
    this->ymin = ymin;
}



void R2Viewport::
Resize(int width, int height)
{
    // Set viewport variables
    this->width = width;
    this->height = height;
}



void R2Viewport::
Resize(int xmin, int ymin, int width, int height)
{
    // Set viewport variables
    this->xmin = xmin;
    this->ymin = ymin;
    this->width = width;
    this->height = height;
}





