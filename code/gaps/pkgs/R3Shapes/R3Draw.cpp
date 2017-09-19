/* Source file for the draw utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Private variables */

int R3draw_mirrored = FALSE;



/* Private variables */

#if (RN_3D_GRFX == RN_3DR)
    static PointF_t R3dr_vertex_normal { 0.0, 0.0, 0.0 };
    static PointF_t R3dr_vertex_texcoords = { 0.0, 0.0, 0.0 };
    static ColorFA_t R3dr_vertex_color = { 0.0, 0.0, 0.0, 0.0 };
#endif




void R4Matrix::
Load(void) const
{
    // Load matrix onto stack replacing top of stack
#if ((RN_3D_GRFX == RN_IRISGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    R4Matrix matrix(Transpose());
    loadmatrix((Matrix) matrix.m);
#elif ((RN_3D_GRFX == RN_OPENGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    R4Matrix matrix(Transpose());
    glLoadMatrixf((const GLfloat *) matrix.m);
#elif ((RN_3D_GRFX == RN_OPENGL) && (RN_MATH_PRECISION == RN_DOUBLE_PRECISION))
    R4Matrix matrix(Transpose());
    glLoadMatrixd((const GLdouble *) matrix.m);
#elif ((RN_3D_GRFX == RN_3DR) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    G3dSetTransform(R3dr_gc, m);
#else
    RNAbort("Not Implemented");
#endif
}



void R4Matrix::
Draw(void) const
{
    // Multiply top of stack by matrix
#if ((RN_3D_GRFX == RN_IRISGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    R4Matrix matrix(Transpose());
    multmatrix((Matrix) matrix.m);
#elif ((RN_3D_GRFX == RN_OPENGL) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    R4Matrix matrix(Transpose());
    glMultMatrixf((const GLfloat *) matrix.m);
#elif ((RN_3D_GRFX == RN_OPENGL) && (RN_MATH_PRECISION == RN_DOUBLE_PRECISION))
    R4Matrix matrix(Transpose());
    glMultMatrixd((const GLdouble *) matrix.m);
#elif ((RN_3D_GRFX == RN_3DR) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    G3dPreMultTransform(R3dr_gc, m);
#else
    RNAbort("Not Implemented");
#endif
}



void R4Matrix::
Push(void) const
{
    // Push top of stack
#if (RN_3D_GRFX == RN_IRISGL)
    pushmatrix(); 
#elif (RN_3D_GRFX == RN_OPENGL)
    glPushMatrix();
#elif (RN_3D_GRFX == RN_3DR)
    G3dPushTransform(R3dr_gc);
#else
    RNAbort("Not Implemented");
#endif

    // Multiply top of stack by matrix
    Draw();
}



void R4Matrix::
Pop(void) const
{
    // Pop top of stack
#if (RN_3D_GRFX == RN_IRISGL)
    popmatrix();
#elif (RN_3D_GRFX == RN_OPENGL)
    glPopMatrix();
#elif (RN_3D_GRFX == RN_3DR)
    G3dPopTransform(R3dr_gc);
#else
    RNAbort("Not Implemented");
#endif
}



void R3Triad::
Draw(void) const
{
    // Draw three rays
    R3BeginLine();
    R3LoadPoint(R3zero_point);
    R3LoadPoint((R3zero_point + Axis(RN_X)));
    R3EndLine();
    R3BeginLine();
    R3LoadPoint(R3zero_point);
    R3LoadPoint((R3zero_point + Axis(RN_Y)));
    R3EndLine();
    R3BeginLine();
    R3LoadPoint(R3zero_point);
    R3LoadPoint((R3zero_point + Axis(RN_Z)));
    R3EndLine();
}



void R3CoordSystem::
Draw(void) const
{
    // Draw three rays
    R3BeginLine();
    R3LoadPoint(Origin());
    R3LoadPoint((Origin() + Axes().Axis(RN_X)));
    R3EndLine();
    R3BeginLine();
    R3LoadPoint(Origin());
    R3LoadPoint((Origin() + Axes().Axis(RN_Y)));
    R3EndLine();
    R3BeginLine();
    R3LoadPoint(Origin());
    R3LoadPoint((Origin() + Axes().Axis(RN_Z)));
    R3EndLine();
}




void R3Affine::
Load(void) const
{
    // Load matrix onto stack
    Matrix().Load();

    // Set mirror flag
    if (IsMirrored()) {
        R3draw_mirrored = !R3draw_mirrored;
#if (RN_3D_GRFX == RN_IRISGL)
        frontface(R3draw_mirrored);
        backface(!R3draw_mirrored);
#elif (RN_3D_GRFX == RN_OPENGL)
        glFrontFace((R3draw_mirrored) ? GL_CW : GL_CCW);
#elif (RN_3D_GRFX == RN_3DR)
	G3dSetState(R3dr_gc, G3DL_FRONT_CCW, (R3draw_mirrored) ? 0 : 1);
#else
	RNAbort("Not Implemented");
#endif
    }
}



void R3Affine::
Draw(void) const
{
    // Multiply top of stack by matrix
    Matrix().Draw();

    // Set mirror flag - this screws up Push/Pop ???
    if (IsMirrored()) {
        R3draw_mirrored = !R3draw_mirrored;
#if (RN_3D_GRFX == RN_IRISGL)
        frontface(R3draw_mirrored);
        backface(!R3draw_mirrored);
#elif (RN_3D_GRFX == RN_OPENGL)
        glFrontFace((R3draw_mirrored) ? GL_CW : GL_CCW);
#elif (RN_3D_GRFX == RN_3DR)
	G3dSetState(R3dr_gc, G3DL_FRONT_CCW, (R3draw_mirrored) ? 0 : 1);
#else
	RNAbort("Not Implemented");
#endif
    }
}



void R3Affine::
Push(void) const
{
    // Push matrix onto stack
    Matrix().Push();

    // Set mirror flag
    if (IsMirrored()) {
        R3draw_mirrored = !R3draw_mirrored;
#if (RN_3D_GRFX == RN_IRISGL)
        frontface(R3draw_mirrored);
        backface(!R3draw_mirrored);
#elif (RN_3D_GRFX == RN_OPENGL)
        glFrontFace((R3draw_mirrored) ? GL_CW : GL_CCW);
#elif (RN_3D_GRFX == RN_3DR)
	G3dSetState(R3dr_gc, G3DL_FRONT_CCW, (R3draw_mirrored) ? 0 : 1);
#else
	RNAbort("Not Implemented");
#endif
    }
}



void R3Affine::
Pop(void) const
{
    // Pop matrix off stack
    Matrix().Pop();

    // Restore mirror flag
    if (IsMirrored()) {
        R3draw_mirrored = !R3draw_mirrored;
#if (RN_3D_GRFX == RN_IRISGL)
        frontface(R3draw_mirrored);
        backface(!R3draw_mirrored);
#elif (RN_3D_GRFX == RN_OPENGL)
        glFrontFace((R3draw_mirrored) ? GL_CW : GL_CCW);
#elif (RN_3D_GRFX == RN_3DR)
	G3dSetState(R3dr_gc, G3DL_FRONT_CCW, (R3draw_mirrored) ? 0 : 1);
#else
	RNAbort("Not Implemented");
#endif
    }
}



void R3Vector::
Draw(void) const
{
    // Draw vector from (0,0,0)
    R3BeginLine();
    R3LoadPoint(R3zero_point);
    R3LoadPoint((R3zero_point + *this));
    R3EndLine();
}



void R3Point::
Draw(void) const
{
    // Draw point
    R3BeginLine();
    R3LoadPoint(v);
    R3LoadPoint(v);
    R3EndLine();
}



void R3Line::
Draw(void) const
{
    // Draw line
    R3BeginLine();
    R3LoadPoint(Point());
    R3LoadPoint(Point() + Vector());
    R3EndLine();
}



void R3Ray::
Draw(void) const
{
    // Draw ray
    R3BeginLine();
    R3LoadPoint(Start());
    R3LoadPoint((Start() + Vector()));
    R3EndLine();
}



void R3Span::
Draw(void) const
{
    // Draw span
    R3BeginLine();
    R3LoadPoint(Start());
    R3LoadPoint(End());
    R3EndLine();
}



void R3Plane::
Draw(void) const
{
    // Draw plane
    RNAbort("Not implemented");
}



void R3Halfspace::
Draw(void) const
{
    // Draw halfspace plane
    Plane().Draw();
}



void R3Triangle::
Draw(const R3DrawFlags draw_flags) const
{
    // Ammend draw flags
    R3DrawFlags flags(Flags() & draw_flags);
    RNDimension dim = 0, dim1 = 0, dim2 = 0;

    // Unroll flags/loops for efficiency 
    switch (flags) {
    case R3_SURFACES_DRAW_FLAG:
      // No shading
      R3BeginPolygon();
      R3LoadPoint(v[0]->Position());
      R3LoadPoint(v[1]->Position());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_SURFACE_NORMALS_DRAW_FLAG:
      // Flat shading
      R3BeginPolygon();
      R3LoadNormal(Normal());
      R3LoadPoint(v[0]->Position());
      R3LoadPoint(v[1]->Position());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_VERTEX_COLORS_DRAW_FLAG:
      // Color interpolation
      R3BeginPolygon();
      RNLoadRgb(v[0]->Color());
      R3LoadPoint(v[0]->Position());
      RNLoadRgb(v[1]->Color());
      R3LoadPoint(v[1]->Position());
      RNLoadRgb(v[2]->Color());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_SURFACE_NORMALS_DRAW_FLAG | R3_VERTEX_NORMALS_DRAW_FLAG:
      // Gouraud shading
      R3BeginPolygon();
      R3LoadNormal(v[0]->Normal());
      R3LoadPoint(v[0]->Position());
      R3LoadNormal(v[1]->Normal());
      R3LoadPoint(v[1]->Position());
      R3LoadNormal(v[2]->Normal());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_SURFACE_NORMALS_DRAW_FLAG | 
         R3_SURFACE_TEXTURE_DRAW_FLAG | R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG: 
      // Flat shading with texture coordinates
      R3BeginPolygon();
      R3LoadNormal(Normal());
      R3LoadTextureCoords(v[0]->TextureCoords());
      R3LoadPoint(v[0]->Position());
      R3LoadTextureCoords(v[1]->TextureCoords());
      R3LoadPoint(v[1]->Position());
      R3LoadTextureCoords(v[2]->TextureCoords());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_SURFACE_NORMALS_DRAW_FLAG | 
         R3_SURFACE_TEXTURE_DRAW_FLAG: 
      // Flat shading with texture coordinate generation
      dim = Normal().MaxDimension();
      dim1 = (dim + 1) % 3;
      dim2 = (dim + 2) % 3;
      R3BeginPolygon();
      R3LoadNormal(Normal());
      R3LoadTextureCoords(v[0]->Position()[dim1], v[0]->Position()[dim2]);
      R3LoadPoint(v[0]->Position());
      R3LoadTextureCoords(v[1]->Position()[dim1], v[1]->Position()[dim2]);
      R3LoadPoint(v[1]->Position());
      R3LoadTextureCoords(v[2]->Position()[dim1], v[2]->Position()[dim2]);
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_SURFACE_NORMALS_DRAW_FLAG | R3_VERTEX_NORMALS_DRAW_FLAG | 
         R3_SURFACE_TEXTURE_DRAW_FLAG | R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG: 
      // Gouraud shading with texture coordinates
      R3BeginPolygon();
      R3LoadTextureCoords(v[0]->TextureCoords());
      R3LoadNormal(v[0]->Normal());
      R3LoadPoint(v[0]->Position());
      R3LoadTextureCoords(v[1]->TextureCoords());
      R3LoadNormal(v[1]->Normal());
      R3LoadPoint(v[1]->Position());
      R3LoadTextureCoords(v[2]->TextureCoords());
      R3LoadNormal(v[2]->Normal());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    case R3_SURFACES_DRAW_FLAG | R3_SURFACE_NORMALS_DRAW_FLAG | R3_VERTEX_NORMALS_DRAW_FLAG | 
         R3_SURFACE_TEXTURE_DRAW_FLAG: 
      // Gouraud shading with texture coordinate generation
      dim = Normal().MaxDimension();
      dim1 = (dim + 1) % 3;
      dim2 = (dim + 2) % 3;
      R3BeginPolygon();
      R3LoadTextureCoords(v[0]->Position()[dim1], v[0]->Position()[dim2]);
      R3LoadNormal(v[0]->Normal());
      R3LoadPoint(v[0]->Position());
      R3LoadTextureCoords(v[1]->Position()[dim1], v[1]->Position()[dim2]);
      R3LoadNormal(v[1]->Normal());
      R3LoadPoint(v[1]->Position());
      R3LoadTextureCoords(v[2]->Position()[dim1], v[2]->Position()[dim2]);
      R3LoadNormal(v[2]->Normal());
      R3LoadPoint(v[2]->Position());
      R3EndPolygon();
      break;

    default:
        // Draw surface
        if (flags[R3_SURFACES_DRAW_FLAG]) {
	    // Begin polygon
	    R3BeginPolygon();

	    // Load polygon normal
	    if ((flags[R3_SURFACE_NORMALS_DRAW_FLAG]) &&
		(!flags[R3_VERTEX_NORMALS_DRAW_FLAG]))
	        R3LoadNormal(Normal());

	    // Compute stuff for texture coordinates
	    if ((draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) &&
		(!flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG])) {
	        dim = Normal().MaxDimension();
		dim1 = (dim + 1) % 3;
		dim2 = (dim + 2) % 3;
	    }

	    // Load triangle vertices
	    for (int i = 0; i < 3; i++) {
		// Load vertex color
  		if (flags[R3_VERTEX_COLORS_DRAW_FLAG]) {
                    RNLoadRgb(v[i]->Color());
 		}

		// Load vertex normal
		if (flags[R3_VERTEX_NORMALS_DRAW_FLAG])
		    R3LoadNormal(v[i]->normal);

		// Load vertex texture coordinates 
		if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) {
		    if (flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) R3LoadTextureCoords(v[i]->TextureCoords());
		    else R3LoadTextureCoords(v[i]->Position()[dim1], v[i]->Position()[dim2]);
		}

		// Load vertex
		R3LoadPoint(v[i]->Position());
	    }

	    // End polygon 
	    R3EndPolygon();
	}
	
	// Draw edges
	if (flags[R3_EDGES_DRAW_FLAG]) {
	    R3BeginLoop();
            R3LoadPoint(v[0]->Position());
            R3LoadPoint(v[1]->Position());
            R3LoadPoint(v[2]->Position());
	    R3EndLoop();
	}
	break;
    }
}



void R3TriangleArray::
Draw(const R3DrawFlags draw_flags) const
{
    // Draw all triangles
    for (int i = 0; i < triangles.NEntries(); i++)
        triangles.Kth(i)->Draw(draw_flags);
}



void R3Circle::
Draw(const R3DrawFlags draw_flags) const
{
    // Push matrix (Zaxis -> ScaledNormal + OffsetToCenter)
    R4Matrix matrix = R4identity_matrix;
    matrix.Translate(center.Vector());
    matrix.Rotate(R3posz_vector, Normal());
    matrix.Push();

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	// Load normal
	if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) 
	    R3LoadNormal(R3posz_vector);

	// Draw circle in XY plane
	R3BeginPolygon();
	for (int i = 0; i < R3circle_npoints; i++) {
	    if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) 
		R3LoadTextureCoords(R3circle_texcoords[i]);
	    R3LoadPoint(R3circle_points[i] * radius);
	}
	R3EndPolygon();
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	R3BeginLoop();
	for (int i = 0; i < R3circle_npoints; i++) 
	    R3LoadPoint(R3circle_points[i] * radius);
	R3EndLoop();
    }

    // Pop matrix
    matrix.Pop();
}



void R3Ellipse::
Draw(const R3DrawFlags draw_flags) const
{
    // Push matrix
    R4Matrix matrix = R4identity_matrix;
    matrix.Transform(cs.Matrix());
    matrix.Scale(R3Vector(radii.X(), radii.Y(), 1.0));
    matrix.Push();

    // Draw unit sphere
    R3unit_circle.Draw(draw_flags);

    // Pop matrix
    matrix.Pop();
}



void R3Rectangle::
Draw(const R3DrawFlags draw_flags) const
{
    // Check if box is empty
    if (IsEmpty()) return;

    /* Get texture coordinates */
    static R2Point texcoords[4] = {
	R2Point(0.0, 0.0),
	R2Point(1.0, 0.0),
	R2Point(1.0, 1.0),
	R2Point(0.0, 1.0)
    };

    /* Get box corner points */
    R3Point corners[4];
    corners[0] = Corner(RN_NN_QUADRANT);
    corners[1] = Corner(RN_PN_QUADRANT);
    corners[2] = Corner(RN_PP_QUADRANT);
    corners[3] = Corner(RN_NP_QUADRANT);

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
        if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) 
            R3LoadNormal(Normal());

        R3BeginPolygon();
        for (int j = 0; j < 4; j++) {
            if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG])
                R3LoadTextureCoords(texcoords[j]);
	    R3LoadPoint(corners[j]);
        }
        R3EndPolygon();
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	R3BeginLoop();
	for (int i = 0; i < 4; i++)
	    R3LoadPoint(corners[i]);
	R3EndLoop();
    }    
}



void R3Box::
Draw(const R3DrawFlags draw_flags) const
{
    static R3Vector normals[6] = {
	R3Vector(-1.0, 0.0, 0.0),
	R3Vector(1.0, 0.0, 0.0),
	R3Vector(0.0, -1.0, 0.0),
	R3Vector(0.0, 1.0, 0.0),
	R3Vector(0.0, 0.0, -1.0),
	R3Vector(0.0, 0.0, 1.0)
    };
    static R2Point texcoords[4] = {
	R2Point(0.0, 0.0),
	R2Point(1.0, 0.0),
	R2Point(1.0, 1.0),
	R2Point(0.0, 1.0)
    };
    static int surface_paths[6][4] = {
	{ 3, 0, 1, 2 },
	{ 4, 7, 6, 5 },
	{ 0, 4, 5, 1 },
	{ 7, 3, 2, 6 },
	{ 3, 7, 4, 0 },
        { 1, 5, 6, 2 }
    };
    static int outline_path[16] = {
	0, 1, 2, 3,
	0, 4, 5, 6,
	7, 4, 5, 1,
	2, 6, 7, 3
    };

    /* Get box corner points */
    R3Point corners[8];
    corners[0] = Corner(RN_NNN_OCTANT);
    corners[1] = Corner(RN_NNP_OCTANT);
    corners[2] = Corner(RN_NPP_OCTANT);
    corners[3] = Corner(RN_NPN_OCTANT);
    corners[4] = Corner(RN_PNN_OCTANT);
    corners[5] = Corner(RN_PNP_OCTANT);
    corners[6] = Corner(RN_PPP_OCTANT);
    corners[7] = Corner(RN_PPN_OCTANT);

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	for (int i = 0; i < 6; i++) {
	    R3BeginPolygon();
	    if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) 
		R3LoadNormal(normals[i]);
	    for (int j = 0; j < 4; j++) {
		if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) 
		    R3LoadTextureCoords(texcoords[j]);
		R3LoadPoint(corners[surface_paths[i][j]]);
	    }
	    R3EndPolygon();
	}
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	R3BeginLine();
	for (int i = 0; i < 16; i++)
	    R3LoadPoint(corners[outline_path[i]]);
	R3EndLine();
    }    
}



void R3OrientedBox::
Draw(const R3DrawFlags draw_flags) const
{
    // Check if box is empty
    if (IsEmpty()) return;

    static R2Point texcoords[4] = {
	R2Point(0.0, 0.0),
	R2Point(1.0, 0.0),
	R2Point(1.0, 1.0),
	R2Point(0.0, 1.0)
    };
    static int surface_paths[6][4] = {
	{ 3, 0, 1, 2 },
	{ 4, 7, 6, 5 },
	{ 0, 4, 5, 1 },
	{ 7, 3, 2, 6 },
	{ 3, 7, 4, 0 },
        { 1, 5, 6, 2 }
    };
    static int outline_path[16] = {
	0, 1, 2, 3,
	0, 4, 5, 6,
	7, 4, 5, 1,
	2, 6, 7, 3
    };

    /* Get box corner points */
    R3Point corners[8];
    corners[0] = Corner(RN_NNN_OCTANT);
    corners[1] = Corner(RN_NNP_OCTANT);
    corners[2] = Corner(RN_NPP_OCTANT);
    corners[3] = Corner(RN_NPN_OCTANT);
    corners[4] = Corner(RN_PNN_OCTANT);
    corners[5] = Corner(RN_PNP_OCTANT);
    corners[6] = Corner(RN_PPP_OCTANT);
    corners[7] = Corner(RN_PPN_OCTANT);

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	for (int i = 0; i < 6; i++) {
	    R3BeginPolygon();
	    if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) {
                RNScalar sign = (i < 3) ? 1.0 : -1.0;
                R3LoadNormal(sign * Axis(i%3));
            }
	    for (int j = 0; j < 4; j++) {
		if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) 
		    R3LoadTextureCoords(texcoords[j]);
		R3LoadPoint(corners[surface_paths[i][j]]);
	    }
	    R3EndPolygon();
	}
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	R3BeginLine();
	for (int i = 0; i < 16; i++)
	    R3LoadPoint(corners[outline_path[i]]);
	R3EndLine();
    }    
}



void R3Cylinder::
Draw(const R3DrawFlags draw_flags) const
{
#if (RN_3D_GRFX == RN_OPENGL)
    // Create GLU quadric
    static GLUquadricObj *cylinder = gluNewQuadric();

    // Push matrix
    R4Matrix matrix = R4identity_matrix;
    matrix.Translate(axis.Start().Vector());
    matrix.Rotate(R3posz_vector, axis.Vector());
    matrix.Push();

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	if (draw_flags[R3_VERTEX_NORMALS_DRAW_FLAG]) gluQuadricNormals(cylinder, (GLenum) GLU_SMOOTH);
	else if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) gluQuadricNormals(cylinder, (GLenum) GLU_FLAT);
	if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) gluQuadricTexture(cylinder, GL_TRUE);
	else gluQuadricTexture(cylinder, GL_FALSE);
	gluQuadricDrawStyle(cylinder, (GLenum) GLU_FILL);
	gluCylinder(cylinder, base.Radius(), top.Radius(), Height(), 16, 1);
	gluQuadricOrientation(cylinder, (GLenum) GLU_INSIDE);
	gluDisk(cylinder, 0.0, Radius(), 16, 1);
	gluQuadricOrientation(cylinder, (GLenum) GLU_OUTSIDE);
	glTranslated(0.0, 0.0, Height());
	gluDisk(cylinder, 0.0, Radius(), 16, 1);
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	gluQuadricNormals(cylinder, (GLenum) GLU_NONE);
	gluQuadricTexture(cylinder, GL_FALSE);
	gluQuadricDrawStyle(cylinder, (GLenum) GLU_SILHOUETTE);
	gluCylinder(cylinder, base.Radius(), top.Radius(), Height(), 16, 1);
    }

    // Pop matrix
    matrix.Pop();
#else
    RNAbort("Not Implemented");
#endif
}



void R3Cone::
Draw(const R3DrawFlags draw_flags) const
{
#if (RN_3D_GRFX == RN_OPENGL)
    // Create GLU quadric
    // Should create GLU quadric for each draw style ???
    static GLUquadricObj *cone = gluNewQuadric(); 

    // Push matrix
    R4Matrix matrix = R4identity_matrix;
    matrix.Translate(axis.Start().Vector());
    matrix.Rotate(R3posz_vector, axis.Vector());
    matrix.Push();

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	if (draw_flags[R3_VERTEX_NORMALS_DRAW_FLAG]) gluQuadricNormals(cone, (GLenum) GLU_SMOOTH);
	else if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) gluQuadricNormals(cone, (GLenum) GLU_FLAT);
	if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) gluQuadricTexture(cone, GL_TRUE);
	else gluQuadricTexture(cone, GL_FALSE);
	gluQuadricDrawStyle(cone, (GLenum) GLU_FILL);
	gluCylinder(cone, Radius(), 0.0, Height(), 16, 1);
	gluQuadricOrientation(cone, (GLenum) GLU_INSIDE);
	gluDisk(cone, 0.0, Radius(), 16, 1);
	gluQuadricOrientation(cone, (GLenum) GLU_OUTSIDE);
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	gluQuadricNormals(cone, (GLenum) GLU_NONE);
	gluQuadricTexture(cone, GL_FALSE);
	gluQuadricDrawStyle(cone, (GLenum) GLU_SILHOUETTE);
	gluCylinder(cone, Radius(), 0.0, Height(), 16, 1);
    }

    // Pop matrix
    matrix.Pop();
#else
    RNAbort("Not Implemented");
#endif
}



void R3Sphere::
Draw(const R3DrawFlags draw_flags) const
{
    // Draw sphere - sphdraw uses n3f only ???
#if (RN_3D_GRFX == RN_IRISGL)
    float sphparams[4];
    sphparams[0] = center.X();
    sphparams[1] = center.Y();
    sphparams[2] = center.Z();
    sphparams[3] = radius;
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	sphmode(SPH_PRIM, SPH_MESH);
	sphdraw(sphparams);
    }
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	sphmode(SPH_PRIM, SPH_LINE);
	sphdraw(sphparams);
    }

#elif (RN_3D_GRFX == RN_OPENGL)
    // Create GLU quadric
    static GLUquadricObj *sphere = gluNewQuadric(); // ???

    // Push matrix 
    R4Matrix matrix = R4identity_matrix;
    matrix.Translate(center.Vector());
    matrix.Push();

    // Draw surface
    if (draw_flags[R3_SURFACES_DRAW_FLAG]) {
	if (draw_flags[R3_VERTEX_NORMALS_DRAW_FLAG]) gluQuadricNormals(sphere, (GLenum) GLU_SMOOTH);
	else if (draw_flags[R3_SURFACE_NORMALS_DRAW_FLAG]) gluQuadricNormals(sphere, (GLenum) GLU_FLAT);
	if (draw_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) gluQuadricTexture(sphere, GL_TRUE);
	else gluQuadricTexture(sphere, GL_FALSE);
	gluQuadricDrawStyle(sphere, (GLenum) GLU_FILL);
	gluSphere(sphere, radius, 8, 8);
    }

    // Draw edges
    if (draw_flags[R3_EDGES_DRAW_FLAG]) {
	gluQuadricNormals(sphere, (GLenum) GLU_NONE);
	gluQuadricTexture(sphere, GL_FALSE);
	gluQuadricDrawStyle(sphere, (GLenum) GLU_SILHOUETTE);
	gluSphere(sphere, radius, 8, 8);
    }

    // Pop matrix
    matrix.Pop();
#else
    RNAbort("Not Implemented");
#endif
}



void R3Ellipsoid::
Draw(const R3DrawFlags draw_flags) const
{
    // Push matrix
    R4Matrix matrix = R4identity_matrix;
    matrix.Transform(cs.Matrix());
    matrix.Scale(radii);
    matrix.Push();

    // Draw unit sphere
    R3unit_sphere.Draw(draw_flags);

    // Pop matrix
    matrix.Pop();
}




