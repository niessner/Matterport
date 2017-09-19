/* Source file for GAPS graphics module */



/* Include files */

#include "RNBasics/RNBasics.h"



/* Xlib variables */

#if (RN_2D_GRFX == RN_XLIB)
    Display *RNgrfx_xdisplay = NULL;
    Window RNgrfx_xwindow = NULL;
    GC RNgrfx_xgc;
    XPoint RNgrfx_xpoints[256];
    int RNgrfx_xnpoints = 0;
#endif



/* 2D viewing variables */

float RNgrfx_window_xscale = 1.0;
float RNgrfx_window_yscale = 1.0;
float RNgrfx_window_width = 1.0;
float RNgrfx_window_height = 1.0;
float RNgrfx_window_xcenter = 0.5;
float RNgrfx_window_ycenter = 0.5;
float RNgrfx_window_xmin = 0.0;
float RNgrfx_window_ymin = 0.0;
float RNgrfx_window_xmax = 1.0;
float RNgrfx_window_ymax = 1.0;
int RNgrfx_viewport_width = 100;
int RNgrfx_viewport_height = 100;
int RNgrfx_viewport_xcenter = 50;
int RNgrfx_viewport_ycenter = 50;
int RNgrfx_viewport_xmin = 0;
int RNgrfx_viewport_ymin = 0;
int RNgrfx_viewport_xmax = 100;
int RNgrfx_viewport_ymax = 100;



/* Private bookkeeping variables */

static int RNgrfx_active_count = 0;



int 
RNInitGrfx(void)
{
    // Check whether are already initialized 
    if ((RNgrfx_active_count++) > 0) return 1;

    // Initialize submodules 
    // ???

    // return OK status 
    return 1;
}



#if (RN_2D_GRFX == RN_XLIB)

int 
RNInitGrfx(Display *display, Window window, GC gc)
{
    // Initialize Xlib variables 
    RNgrfx_xdisplay = display;
    RNgrfx_xwindow = window;
    RNgrfx_xgc = gc;

    // Initialize package
    return RNInitGrfx();
}

#endif



void 
RNStopGrfx(void)
{
    // Check whether have been initialized 
    if ((--RNgrfx_active_count) > 0) return;

    // Stop submodules 
    // ???
}



/* 2D viewing functions */

void 
R2SetViewport(int xmin, int ymin, int xmax, int ymax)
{
    /* Update viewport */
    RNgrfx_viewport_xmin = xmin;
    RNgrfx_viewport_ymin = ymin;
    RNgrfx_viewport_xmax = xmax;
    RNgrfx_viewport_ymax = ymax;
    RNgrfx_viewport_xcenter = (RNgrfx_viewport_xmax + RNgrfx_viewport_xmin) / 2;
    RNgrfx_viewport_ycenter = (RNgrfx_viewport_ymax + RNgrfx_viewport_ymin) / 2;
    RNgrfx_viewport_width = RNgrfx_viewport_xmax - RNgrfx_viewport_xmin;
    RNgrfx_viewport_height = RNgrfx_viewport_ymax - RNgrfx_viewport_ymin;

    /* Update scale */
    RNgrfx_window_xscale = RNgrfx_viewport_width / RNgrfx_window_width;
    RNgrfx_window_yscale = RNgrfx_viewport_height / RNgrfx_window_height;
}



void 
R2SetWindow(float xmin, float ymin, float xmax, float ymax)
{
    /* Update window */
    RNgrfx_window_xmin = xmin;
    RNgrfx_window_ymin = ymin;
    RNgrfx_window_xmax = xmax;
    RNgrfx_window_ymax = ymax;
    RNgrfx_window_xcenter = (RNgrfx_window_xmax + RNgrfx_window_xmin) / 2.0;
    RNgrfx_window_ycenter = (RNgrfx_window_ymax + RNgrfx_window_ymin) / 2.0;
    RNgrfx_window_width = RNgrfx_window_xmax - RNgrfx_window_xmin;
    RNgrfx_window_height = RNgrfx_window_ymax - RNgrfx_window_ymin;

    /* Update scale */
    RNgrfx_window_xscale = RNgrfx_viewport_width / RNgrfx_window_width;
    RNgrfx_window_yscale = RNgrfx_viewport_height / RNgrfx_window_height;
}



void 
R2ScaleWindow(float xscale, float yscale)
{
    /* Update window */
    RNgrfx_window_xmin = RNgrfx_window_xcenter - 0.5 * xscale * RNgrfx_window_width;
    RNgrfx_window_ymin = RNgrfx_window_ycenter - 0.5 * yscale * RNgrfx_window_height;
    RNgrfx_window_xmax = RNgrfx_window_xcenter + 0.5 * xscale * RNgrfx_window_width;
    RNgrfx_window_ymax = RNgrfx_window_ycenter + 0.5 * yscale * RNgrfx_window_height;
    RNgrfx_window_width = RNgrfx_window_xmax - RNgrfx_window_xmin;
    RNgrfx_window_height = RNgrfx_window_ymax - RNgrfx_window_ymin;

    /* Update scale */
    RNgrfx_window_xscale = RNgrfx_viewport_width / RNgrfx_window_width;
    RNgrfx_window_yscale = RNgrfx_viewport_height / RNgrfx_window_height;
}



void 
R2TranslateWindow(float xtranslate, float ytranslate)
{
    /* Update window */
    RNgrfx_window_xmin += xtranslate;
    RNgrfx_window_ymin += ytranslate;
    RNgrfx_window_xmax += xtranslate;
    RNgrfx_window_ymax += ytranslate;
    RNgrfx_window_xcenter = (RNgrfx_window_xmax + RNgrfx_window_xmin) / 2.0;
    RNgrfx_window_ycenter = (RNgrfx_window_ymax + RNgrfx_window_ymin) / 2.0;
}



void 
RNGrfxError(const char *message)
{
    // Fatal error
    fprintf(stderr, "%s\n", message);
    abort();
}


