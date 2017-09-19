/* Global Xlib variables */

#if (RN_2D_GRFX == RN_XLIB)
    extern Display *RNgrfx_xdisplay;
    extern Window RNgrfx_xwindow;
    extern GC RNgrfx_xgc;
    extern XPoint RNgrfx_xpoints[];
    extern int RNgrfx_xnpoints;
#endif



/* Global 2D viewing variables */

extern float RNgrfx_window_xscale;
extern float RNgrfx_window_yscale;
extern float RNgrfx_window_width;
extern float RNgrfx_window_height;
extern float RNgrfx_window_xcenter;
extern float RNgrfx_window_ycenter;
extern float RNgrfx_window_xmin;
extern float RNgrfx_window_ymin;
extern float RNgrfx_window_xmax;
extern float RNgrfx_window_ymax;
extern int RNgrfx_viewport_width;
extern int RNgrfx_viewport_height;
extern int RNgrfx_viewport_xcenter;
extern int RNgrfx_viewport_ycenter;
extern int RNgrfx_viewport_xmin;
extern int RNgrfx_viewport_ymin;
extern int RNgrfx_viewport_xmax;
extern int RNgrfx_viewport_ymax;



/* Inline functions */

inline void
R2WindowToViewport(float wx, float wy, int *vx, int *vy)
{
    // Transform x coordinate by window/viewport
    *vx = RNgrfx_viewport_xmin + (int) ((wx - RNgrfx_window_xmin) * RNgrfx_window_xscale);
    *vy = RNgrfx_viewport_ymin + (int) ((wy - RNgrfx_window_ymin) * RNgrfx_window_yscale);
    // *vy = RNgrfx_viewport_ymax - *vy;
}
	


inline void
R2WindowToViewport(double wx, double wy, int *vx, int *vy)
{
    // Transform x coordinate by window/viewport
    *vx = RNgrfx_viewport_xmin + (int) ((wx - RNgrfx_window_xmin) * RNgrfx_window_xscale);
    *vy = RNgrfx_viewport_ymin + (int) ((wy - RNgrfx_window_ymin) * RNgrfx_window_yscale);
    // *vy = RNgrfx_viewport_ymax - *vy;
}
	


inline void
R2WindowToViewport(float wx, float wy, short *vx, short *vy)
{
    // Transform x coordinate by window/viewport
    *vx = (short) (RNgrfx_viewport_xmin + (wx - RNgrfx_window_xmin) * RNgrfx_window_xscale);
    *vy = (short) (RNgrfx_viewport_ymin + (wy - RNgrfx_window_ymin) * RNgrfx_window_yscale);
    // *vy = RNgrfx_viewport_ymax - *vy;
}
	


inline void
R2WindowToViewport(double wx, double wy, short *vx, short *vy)
{
    // Transform x coordinate by window/viewport
    *vx = (short) (RNgrfx_viewport_xmin + (wx - RNgrfx_window_xmin) * RNgrfx_window_xscale);
    *vy = (short) (RNgrfx_viewport_ymin + (wy - RNgrfx_window_ymin) * RNgrfx_window_yscale);
    // *vy = RNgrfx_viewport_ymax - *vy;
}
	


inline void
R2ViewportToWindow(int vx, int vy, float *wx, float *wy)
{
    // Transform x coordinate by window/viewport
    *wx = RNgrfx_window_xmin + (float) ((vx - RNgrfx_viewport_xmin) / RNgrfx_window_xscale);
    *wy = RNgrfx_window_ymin + (float) ((vy - RNgrfx_viewport_ymin) / RNgrfx_window_yscale);
    // *wy = RNgrfx_window_ymax - (*wy - RNgrfx_window_ymin);
}
	


inline void
R2ViewportToWindow(int vx, int vy, double *wx, double *wy)
{
    // Transform x coordinate by window/viewport
    *wx = RNgrfx_window_xmin + (double) ((vx - RNgrfx_viewport_xmin) / RNgrfx_window_xscale);
    *wy = RNgrfx_window_ymin + (double) ((vy - RNgrfx_viewport_ymin) / RNgrfx_window_yscale);
    // *wy = RNgrfx_window_ymax - (*wy - RNgrfx_window_ymin);
}
	


inline void 
R2BeginPolygon(void)
{
    // Begin drawing polygon (follow with R2DrawPoint* and R2EndPolygon)
#if (RN_2D_GRFX == RN_IRISGL)
    bgnpolygon();
#elif (RN_2D_GRFX == RN_OPENGL)
    glBegin(GL_POLYGON);
#elif (RN_2D_GRFX == RN_XLIB)
    RNgrfx_xnpoints = 0;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2EndPolygon(void)
{
    // End drawing polygon
#if (RN_2D_GRFX == RN_IRISGL)
    endpolygon();
#elif (RN_2D_GRFX == RN_OPENGL)
    glEnd();
#elif (RN_2D_GRFX == RN_XLIB)
    XFillPolygon(RNgrfx_xdisplay, RNgrfx_xwindow, RNgrfx_xgc,
	RNgrfx_xpoints, RNgrfx_xnpoints, Convex, CoordModeOrigin);
#else
    RNGrfxError("Not Implemented");
#endif
}
	


inline void 
R2BeginLine(void)
{
    // Begin drawing line (follow with R2DrawPoint* and R2EndLine)
#if (RN_2D_GRFX == RN_IRISGL)
    bgnline();
#elif (RN_2D_GRFX == RN_OPENGL)
    glBegin(GL_LINE_STRIP);
#elif (RN_2D_GRFX == RN_XLIB)
    RNgrfx_xnpoints = 0;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2EndLine(void)
{
    // End drawing line
#if (RN_2D_GRFX == RN_IRISGL)
    endline();
#elif (RN_2D_GRFX == RN_OPENGL)
    glEnd();
#elif (RN_2D_GRFX == RN_XLIB)
    XDrawLines(RNgrfx_xdisplay, RNgrfx_xwindow, RNgrfx_xgc,
	RNgrfx_xpoints, RNgrfx_xnpoints, CoordModeOrigin);
#else
    RNGrfxError("Not Implemented");
#endif
}
	


inline void 
R2BeginLoop(void)
{
    // Begin drawing closed line (follow with R2DrawPoint* and R2EndLoop)
#if (RN_2D_GRFX == RN_IRISGL)
    bgnclosedline();
#elif (RN_2D_GRFX == RN_OPENGL)
    glBegin(GL_LINE_LOOP);
#elif (RN_2D_GRFX == RN_XLIB)
    RNgrfx_xnpoints = 0;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2EndLoop(void)
{
    // End drawing closed line
#if (RN_2D_GRFX == RN_IRISGL)
    endclosedline();
#elif (RN_2D_GRFX == RN_OPENGL)
    glEnd();
#elif (RN_2D_GRFX == RN_XLIB)
    XDrawLine(RNgrfx_xdisplay, RNgrfx_xwindow, RNgrfx_xgc,
	RNgrfx_xpoints[RNgrfx_xnpoints-1].x, RNgrfx_xpoints[RNgrfx_xnpoints-1].y,
	RNgrfx_xpoints[0].x, RNgrfx_xpoints[0].y);
    XDrawLines(RNgrfx_xdisplay, RNgrfx_xwindow, RNgrfx_xgc,
	RNgrfx_xpoints, RNgrfx_xnpoints, CoordModeOrigin);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2LoadPoint(float x, float y)
{
    // Load vertex (within R2BeginXXX and R2EndXXX)
#if (RN_2D_GRFX == RN_IRISGL)
    float point[2];
    point[0] = x; 
    point[1] = y;
    v2f(point);
#elif (RN_2D_GRFX == RN_OPENGL)
    glVertex2f(x, y);
#elif (RN_2D_GRFX == RN_XLIB)
    R2WindowToViewport(x, y, &RNgrfx_xpoints[RNgrfx_xnpoints].x, &RNgrfx_xpoints[RNgrfx_xnpoints].y);
    RNgrfx_xnpoints++;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2LoadPoint(double x, double y)
{
    // Load vertex (within R2BeginXXX and R2EndXXX)
#if (RN_2D_GRFX == RN_IRISGL)
    double point[2];
    point[0] = x; 
    point[1] = y;
    v2d(point);
#elif (RN_2D_GRFX == RN_OPENGL)
    glVertex2d(x, y);
#elif (RN_2D_GRFX == RN_XLIB)
    R2WindowToViewport(x, y, &RNgrfx_xpoints[RNgrfx_xnpoints].x, &RNgrfx_xpoints[RNgrfx_xnpoints].y);
    RNgrfx_xnpoints++;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2LoadPoint(const float point[2])
{
    // Load vertex (within R2BeginXXX and R2EndXXX)
#if (RN_2D_GRFX == RN_IRISGL)
    v2f(point);
#elif (RN_2D_GRFX == RN_OPENGL)
    glVertex2fv(point);
#elif (RN_2D_GRFX == RN_XLIB)
    R2WindowToViewport(point[0], point[1], &RNgrfx_xpoints[RNgrfx_xnpoints].x, &RNgrfx_xpoints[RNgrfx_xnpoints].y);
    RNgrfx_xnpoints++;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R2LoadPoint(const double point[2])
{
    // Load vertex (within R2BeginXXX and R2EndXXX)
#if (RN_2D_GRFX == RN_IRISGL)
    v2d(point);
#elif (RN_2D_GRFX == RN_OPENGL)
    glVertex2dv(point);
#elif (RN_2D_GRFX == RN_XLIB)
    R2WindowToViewport(point[0], point[1], &RNgrfx_xpoints[RNgrfx_xnpoints].x, &RNgrfx_xpoints[RNgrfx_xnpoints].y);
    RNgrfx_xnpoints++;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void
R2DrawImage(int x, int y, int width, int height, int depth, const unsigned char *data)
{
#if (RN_3D_GRFX == RN_OPENGL)
    // Set projection matrix
    glMatrixMode(GL_PROJECTION);  
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height);

    // Set model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Set position for image
    glRasterPos2i(x, y);

    // Determine image format
    GLenum format = GL_LUMINANCE;
    if (depth == 1) format = GL_LUMINANCE;
    else if (depth == 2) format = GL_LUMINANCE_ALPHA;
    else if (depth == 3) format = GL_RGB;
    else if (depth == 4) format = GL_RGBA;
    else RNAbort("Illegal image image");

    // Draw pixels
    glDrawPixels(width, height, format, GL_UNSIGNED_BYTE, data);
  
    // Reset model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // Reset projection matrix
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void
R2DrawText(float x, float y, const char *str)
{
#if (RN_3D_GRFX == RN_OPENGL)
    glRasterPos2f(x, y);
    glCallLists(strlen(str), GL_UNSIGNED_BYTE, (GLubyte *) str);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void
R2DrawText(double x, double y, const char *str)
{
#if (RN_3D_GRFX == RN_OPENGL)
    glRasterPos2d(x, y);
    glCallLists(strlen(str), GL_UNSIGNED_BYTE, (GLubyte *) str);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3BeginPolygon(void)
{
    // Begin drawing polygon (follow with R3DrawPoint* and R3EndPolygon)
#if (RN_3D_GRFX == RN_IRISGL)
    bgnpolygon();
#elif (RN_3D_GRFX == RN_OPENGL)
    glBegin(GL_POLYGON);
#elif (RN_3D_GRFX == RN_3DR)
    G3dBeginPrim(G3D_PRM_POLYGON, 256);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3EndPolygon(void)
{
    // End drawing polygon
#if (RN_3D_GRFX == RN_IRISGL)
    endpolygon();
#elif (RN_3D_GRFX == RN_OPENGL)
    glEnd();
#elif (RN_3D_GRFX == RN_3DR)
    G3dEndPrim(R3dr_gc);
#else
    RNGrfxError("Not Implemented");
#endif
}
	


inline void 
R3BeginLine(void)
{
    // Begin drawing line (follow with R3DrawPoint* and R3EndLine)
#if (RN_3D_GRFX == RN_IRISGL)
    bgnline();
#elif (RN_3D_GRFX == RN_OPENGL)
    glBegin(GL_LINE_STRIP);
#elif (RN_3D_GRFX == RN_3DR)
    G3dBeginPrim(G3D_PRM_POLYLINE, 256);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3EndLine(void)
{
    // End drawing line
#if (RN_3D_GRFX == RN_IRISGL)
    endline();
#elif (RN_3D_GRFX == RN_OPENGL)
    glEnd();
#elif (RN_3D_GRFX == RN_3DR)
    G3dEndPrim(R3dr_gc);
#else
    RNGrfxError("Not Implemented");
#endif
}
	


inline void 
R3BeginLoop(void)
{
    // Begin drawing closed line (follow with R3DrawPoint* and R3EndLoop)
#if (RN_3D_GRFX == RN_IRISGL)
    bgnclosedline();
#elif (RN_3D_GRFX == RN_OPENGL)
    glBegin(GL_LINE_LOOP);
#elif (RN_3D_GRFX == RN_3DR)
    R3BeginLine(); // ???
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3EndLoop(void)
{
    // End drawing closed line
#if (RN_3D_GRFX == RN_IRISGL)
    endclosedline();
#elif (RN_3D_GRFX == RN_OPENGL)
    glEnd();
#elif (RN_3D_GRFX == RN_3DR)
    R3EndLine(); // ???
#else
    RNGrfxError("Not Implemented");
#endif
}
	


inline void 
R3LoadPoint(float x, float y, float z)
{
    // Load vertex (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    float point[3];
    point[0] = x; 
    point[1] = y; 
    point[2] = z;
    v3f(point);
#elif (RN_3D_GRFX == RN_OPENGL)
    glVertex3f(x, y, z);
#elif (RN_3D_GRFX == RN_3DR)
    // Normal, texcoords, color ???
    float point[3];
    point[0] = x; 
    point[1] = y; 
    point[2] = z;
    G3dAddPrimVtxF((PointF_t *) point, NULL, NULL, NULL);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadPoint(double x, double y, double z)
{
    // Load vertex (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    double point[3];
    point[0] = x; 
    point[1] = y; 
    point[2] = z;
    v3f(point);
#elif (RN_3D_GRFX == RN_OPENGL)
    glVertex3d(x, y, z);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadPoint(const float point[3])
{
    // Load vertex (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    v3f(point);
#elif (RN_3D_GRFX == RN_OPENGL)
    glVertex3fv(point);
#elif (RN_3D_GRFX == RN_3DR)
    // Normal, texcoords, color ???
    G3dAddPrimVtxF((PointF_t *) point, NULL, NULL, NULL);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadPoint(const double point[3])
{
    // Load vertex (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    v3d(point);
#elif (RN_3D_GRFX == RN_OPENGL)
    glVertex3dv(point);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadNormal(float x, float y, float z)
{
    // Load normal vector 
#if (RN_3D_GRFX == RN_IRISGL)
    float normal[3];
    normal[0] = x; 
    normal[1] = y; 
    normal[2] = z;
    n3f(normal);
#elif (RN_3D_GRFX == RN_OPENGL)
    glNormal3f(x, y, z);
#elif (RN_3D_GRFX == RN_3DR)
    R3dr_vertex_normal.x = x;
    R3dr_vertex_normal.y = y;
    R3dr_vertex_normal.z = z;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadNormal(double x, double y, double z)
{
    // Load normal vector 
#if (RN_3D_GRFX == RN_IRISGL)
    double normal[3];
    normal[0] = x; 
    normal[1] = y; 
    normal[2] = z;
    n3d(normal);
#elif (RN_3D_GRFX == RN_OPENGL)
    glNormal3d(x, y, z);
#elif (RN_3D_GRFX == RN_3DR)
    R3dr_vertex_normal.x = x;
    R3dr_vertex_normal.y = y;
    R3dr_vertex_normal.z = z;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadNormal(const float normal[3])
{
    // Load normal vector 
#if (RN_3D_GRFX == RN_IRISGL)
    n3f(normal);
#elif (RN_3D_GRFX == RN_OPENGL)
    glNormal3fv(normal);
#elif (RN_3D_GRFX == RN_3DR)
    R3dr_vertex_normal.x = normal[0];
    R3dr_vertex_normal.y = normal[1];
    R3dr_vertex_normal.z = normal[2];
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadNormal(const double normal[3])
{
    // Load normal vector 
#if (RN_3D_GRFX == RN_IRISGL)
    n3d(normal);
#elif (RN_3D_GRFX == RN_OPENGL)
    glNormal3dv(normal);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadTextureCoords(float x, float y)
{
    // Load texture coordinate (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    float texcoords[2];
    texcoords[0] = x; 
    texcoords[1] = y; 
    t2f(texcoords);
#elif (RN_3D_GRFX == RN_OPENGL)
    glTexCoord2f(x, y);
#elif (RN_3D_GRFX == RN_3DR)
    R3dr_vertex_texcoords.x = x;
    R3dr_vertex_texcoords.y = y;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadTextureCoords(double x, double y)
{
    // Load texture coordinate (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    double texcoords[2];
    texcoords[0] = x; 
    texcoords[1] = y; 
    t2d(texcoords);
#elif (RN_3D_GRFX == RN_OPENGL)
    glTexCoord2d(x, y);
#elif (RN_3D_GRFX == RN_3DR)
    R3dr_vertex_texcoords.x = x;
    R3dr_vertex_texcoords.y = y;
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadTextureCoords(const float texcoords[2])
{
    // Load texture coordinate (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    t2f(texcoords);
#elif (RN_3D_GRFX == RN_OPENGL)
    glTexCoord2fv(texcoords);
#elif (RN_3D_GRFX == RN_3DR)
    R3dr_vertex_texcoords.x = texcoords[0];
    R3dr_vertex_texcoords.y = texcoords[1];
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadTextureCoords(const double texcoords[2])
{
    // Set texture coordinate (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    t2d(texcoords);
#elif (RN_3D_GRFX == RN_OPENGL)
    glTexCoord2dv(texcoords);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadTextureCoords(float x, float y, float z)
{
    // Load texture coordinate (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    float texcoords[3];
    texcoords[0] = x; 
    texcoords[1] = y; 
    texcoords[2] = z; 
    t3f(texcoords);
#elif (RN_3D_GRFX == RN_OPENGL)
    glTexCoord3f(x, y, z);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
R3LoadTextureCoords(double x, double y, double z)
{
    // Load texture coordinate (within R3BeginXXX and R3EndXXX)
#if (RN_3D_GRFX == RN_IRISGL)
    double texcoords[2];
    texcoords[0] = x; 
    texcoords[1] = y; 
    texcoords[2] = z; 
    t2d(texcoords);
#elif (RN_3D_GRFX == RN_OPENGL)
    glTexCoord3d(x, y, z);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void
R3DrawText(float x, float y, float z, const char *str)
{
#if (RN_3D_GRFX == RN_OPENGL)
    glRasterPos3f(x, y, z);
    glCallLists(strlen(str), GL_UNSIGNED_BYTE, (GLubyte *) str);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void
R3DrawText(double x, double y, double z, const char *str)
{
#if (RN_3D_GRFX == RN_OPENGL)
    glRasterPos3d(x, y, z);
    glCallLists(strlen(str), GL_UNSIGNED_BYTE, (GLubyte *) str);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
RNLoadRgb(const float rgb[3])
{
    // Load rgb
#if (RN_3D_GRFX == RN_IRISGL)
    c3f(rgb);
#elif (RN_3D_GRFX == RN_OPENGL)
    glColor3fv(rgb);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
RNLoadRgb(double red, double green, double blue)
{
    // Load rgb
#if (RN_3D_GRFX == RN_OPENGL)
    glColor3d(red, green, blue);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
RNLoadRgb(float red, float green, float blue)
{
    // Load rgb
#if (RN_3D_GRFX == RN_OPENGL)
    glColor3f(red, green, blue);
#else
    RNGrfxError("Not Implemented");
#endif
}



inline void 
RNLoadRgb(const double rgb[3])
{
    // Load rgb
#if (RN_3D_GRFX == RN_IRISGL)
    // Assume rgb values are [0,1]
    assert((rgb[0] >= 0.0) && rgb[0] <= 1.0);
    unsigned long r = (unsigned long) (rgb[0] * 255.0);
    assert((rgb[1] >= 0.0) && rgb[1] <= 1.0);
    unsigned long g = (unsigned long) (rgb[1] * 255.0);
    assert((rgb[2] >= 0.0) && rgb[2] <= 1.0);
    unsigned long b = (unsigned long) (rgb[2] * 255.0);
    cpack(0xFF000000 | (b << 24) && 0x00FF0000 | (g << 16) && 0x0000FF00 | (r << 24) && 0x000000FF);
#elif (RN_3D_GRFX == RN_OPENGL)
    glColor3dv(rgb);
#else
    RNGrfxError("Not Implemented");
#endif
}




