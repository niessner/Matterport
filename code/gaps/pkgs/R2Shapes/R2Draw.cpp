/* Source file for the draw utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



void R3Matrix::
Load(void) const
{
    // Load matrix onto stack replacing top of stack
#if ((RN_3D_GRFX == RN_IRISGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    // R3Matrix matrix(Transpose());
    // loadmatrix((Matrix) matrix.m);
#elif (RN_3D_GRFX == RN_OPENGL) 
    double t[4][4];
    t[0][0] = m[0][0]; t[0][1] = m[1][0]; t[0][2] = 0.0;     t[0][3] = 0.0;
    t[1][0] = m[0][1]; t[1][1] = m[1][1]; t[1][2] = 0.0;     t[1][3] = 0.0;
    t[2][0] = 0.0;     t[2][1] = 0.0;     t[2][2] = 1.0;     t[2][3] = 0.0;
    t[3][0] = m[0][2]; t[3][1] = m[1][2]; t[3][2] = 0.0;     t[3][3] = 1.0;
    glLoadMatrixd((const GLdouble *) t);
#else
    RNAbort("Not Implemented");
#endif
}



void R3Matrix::
Draw(void) const
{
    // Multiply top of stack by matrix
#if ((RN_3D_GRFX == RN_IRISGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    // R3Matrix matrix(Transpose());
    // multmatrix((Matrix) matrix.m);
#elif ((RN_3D_GRFX == RN_OPENGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    GLdouble t[4][4];
    t[0][0] = m[0][0]; t[0][1] = m[1][0]; t[0][2] = 0.0;     t[0][3] = 0.0;
    t[1][0] = m[0][1]; t[1][1] = m[1][1]; t[1][2] = 0.0;     t[1][3] = 0.0;
    t[2][0] = 0.0;     t[2][1] = 0.0;     t[2][2] = 1.0;     t[2][3] = 0.0;
    t[3][0] = m[0][2]; t[3][1] = m[1][2]; t[3][2] = 0.0;     t[3][3] = 1.0;
    glMultMatrixd((const GLdouble *) t);
#elif ((RN_3D_GRFX == RN_OPENGL) && (RN_MATH_PRECISION == RN_DOUBLE_PRECISION))
    GLdouble t[4][4];
    t[0][0] = m[0][0]; t[0][1] = m[1][0]; t[0][2] = 0.0;     t[0][3] = 0.0;
    t[1][0] = m[0][1]; t[1][1] = m[1][1]; t[1][2] = 0.0;     t[1][3] = 0.0;
    t[2][0] = 0.0;     t[2][1] = 0.0;     t[2][2] = 1.0;     t[2][3] = 0.0;
    t[3][0] = m[0][2]; t[3][1] = m[1][2]; t[3][2] = 0.0;     t[3][3] = 1.0;
    glMultMatrixd((const GLdouble *) t);
#else
    RNAbort("Not Implemented");
#endif
}



void R3Matrix::
Push(void) const
{
    // Push top of stack
#if (RN_3D_GRFX == RN_IRISGL)
    pushmatrix(); 
#elif (RN_3D_GRFX == RN_OPENGL)
    glPushMatrix();
#else
    RNAbort("Not Implemented");
#endif

    // Multiply top of stack by matrix
    Draw();
}



void R3Matrix::
Pop(void) const
{
    // Pop top of stack
#if (RN_3D_GRFX == RN_IRISGL)
    popmatrix();
#elif (RN_3D_GRFX == RN_OPENGL)
    glPopMatrix();
#else
    RNAbort("Not Implemented");
#endif
}



void R2Affine::
Load(void) const
{
    // Load matrix onto stack
    Matrix().Load();

    // Set mirror flag
    if (IsMirrored()) {}
}



void R2Affine::
Draw(void) const
{
    // Multiply top of stack by matrix
    Matrix().Draw();

    // Set mirror flag
    if (IsMirrored()) {}
}



void R2Affine::
Push(void) const
{
    // Push matrix onto stack
    Matrix().Push();

    // Set mirror flag
    if (IsMirrored()) {}
}



void R2Affine::
Pop(void) const
{
    // Pop matrix off stack
    Matrix().Pop();

    // Restore mirror flag
    if (IsMirrored()) {}
}



void R2Vector::
Draw(void) const
{
    // Draw vector from (0,0,0)
    R2BeginLine();
    R2LoadPoint(R2zero_point);
    R2LoadPoint((R2zero_point + *this));
    R2EndLine();
}



void R2Point::
Draw(void) const
{
    // Draw point
    R2BeginLine();
    R2LoadPoint(v);
    R2LoadPoint(v);
    R2EndLine();
}



void R2Line::
Draw(void) const
{
    // Draw line
    R2BeginLine();
    R2LoadPoint(Point());
    R2LoadPoint(Point() + Vector());
    R2EndLine();
}



void R2Ray::
Draw(void) const
{
    // Draw ray
    R2BeginLine();
    R2LoadPoint(Start());
    R2LoadPoint((Start() + Vector()));
    R2EndLine();
}



void R2Span::
Draw(void) const
{
    // Draw span
    R2BeginLine();
    R2LoadPoint(Start());
    R2LoadPoint(End());
    R2EndLine();
}



void R2Halfspace::
Draw(void) const
{
    // Draw halfspace line
    Line().Draw();
}



void R2Arc::
Draw(const R2DrawFlags draw_flags) const
{
    // Draw arc
    if (draw_flags[R2_EDGES_DRAW_FLAG]) {
	R2BeginLine();
	R2LoadPoint(StartPoint());
	int start_index = (int) (1 + R2circle_npoints * StartAngle() / RN_TWO_PI);
	int stop_index = (int) (R2circle_npoints * StopAngle() / RN_TWO_PI);
	for (int i = start_index; i <= stop_index; i++) 
	    R2LoadPoint(Center() + R2circle_points[i % R2circle_npoints] * Radius());
	R2LoadPoint(StopPoint());
	R2EndLine();
    }
}



void R2Polyline::
Draw(const R2DrawFlags draw_flags) const
{
    // Draw polyline
    R2BeginLine();
    for (int i = 0; i < npoints; i++) R2LoadPoint(points[i]);
    R2EndLine();
}




void R2Box::
Draw(const R2DrawFlags draw_flags) const
{
    // Draw surface
    if (draw_flags[R2_SURFACES_DRAW_FLAG]) {
        R2BeginPolygon();
        R2LoadPoint(Corner(RN_LO, RN_LO));
        R2LoadPoint(Corner(RN_HI, RN_LO));
        R2LoadPoint(Corner(RN_HI, RN_HI));
        R2LoadPoint(Corner(RN_LO, RN_HI));
	R2EndPolygon();
    }

    // Draw edges
    if (draw_flags[R2_EDGES_DRAW_FLAG]) {
	R2BeginLoop();
        R2LoadPoint(Corner(RN_LO, RN_LO));
        R2LoadPoint(Corner(RN_HI, RN_LO));
        R2LoadPoint(Corner(RN_HI, RN_HI));
        R2LoadPoint(Corner(RN_LO, RN_HI));
	R2EndLoop();
    }    
}



void R2Circle::
Draw(const R2DrawFlags draw_flags) const
{
    // Draw surface
    if (draw_flags[R2_SURFACES_DRAW_FLAG]) {
	R2BeginPolygon();
	for (int i = 0; i < R2circle_npoints; i++) 
	    R2LoadPoint(center + R2circle_points[i] * radius);
	R2EndPolygon();
    }

    // Draw edges
    if (draw_flags[R2_EDGES_DRAW_FLAG]) {
	R2BeginLoop();
	for (int i = 0; i < R2circle_npoints; i++) 
	    R2LoadPoint(center + R2circle_points[i] * radius);
	R2EndLoop();
    }
}




void R2Polygon::
Draw(const R2DrawFlags draw_flags) const
{
    // Draw surface
    if (draw_flags[R2_SURFACES_DRAW_FLAG]) {
	R2BeginPolygon();
	for (int i = 0; i < npoints; i++) R2LoadPoint(points[i]);
	R2EndPolygon();
    }

    // Draw edges
    if (draw_flags[R2_EDGES_DRAW_FLAG]) {
	R2BeginLoop();
	for (int i = 0; i < npoints; i++) R2LoadPoint(points[i]);
	R2EndLoop();
    }
}




