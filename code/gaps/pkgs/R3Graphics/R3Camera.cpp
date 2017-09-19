/* Source file for the R3 camera class */



/* Include files */

#include "R3Graphics.h"



/* Public variables */

RNLength R3default_camera_xfov(RN_PI / 6.0);
RNLength R3default_camera_yfov(RN_PI / 6.0);
R3Vector R3default_camera_right_vector(1.0, 0.0, 0.0);
R3Vector R3default_camera_towards_vector(0.0, 1.0, 0.0);
R3Vector R3default_camera_up_vector(0.0, 0.0, 1.0);
R3Camera R3default_camera(R3CoordSystem(R3Point(0,0,0),
    R3Triad(R3default_camera_towards_vector, R3default_camera_up_vector)),
    R3default_camera_xfov, R3default_camera_yfov, 0.001, 1000.0);



/* Public functions */

int 
R3InitCamera()
{
    /* Return success */
    return TRUE;
}



void 
R3StopCamera()
{
}



R3Camera::
R3Camera(void)
    : cs(R3Point(0,0,0), R3Triad(R3Vector(1,0,0), R3Vector(0,1,0), R3Vector(0,0,1))),
      xfov(0),
      yfov(0),
      neardist(0),
      fardist(0),
      value(0),
      data(NULL)
{
}



R3Camera::
R3Camera(const R3CoordSystem& cs,
	  RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist)
    : cs(cs),
      xfov(xfov),
      yfov(yfov),
      neardist(neardist),
      fardist(fardist),
      value(0),
      data(NULL)
{
}



R3Camera::
R3Camera(const R3Camera& camera)
    : cs(camera.cs),
      xfov(camera.xfov),
      yfov(camera.yfov),
      neardist(camera.neardist),
      fardist(camera.fardist),
      value(camera.value),
      data(NULL)
{
}



R3Camera::
R3Camera(const R3Point& origin, const R3Triad& triad,
		  RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist)
    : cs(origin, triad),
      xfov(xfov),
      yfov(yfov),
      neardist(neardist),
      fardist(fardist),
      value(0),
      data(NULL)
{
}



R3Camera::
R3Camera(const R3Point& origin, const R3Vector& towards, const R3Vector& up,
		  RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist)
    : cs(origin, R3Triad(towards, up)),
      xfov(xfov),
      yfov(yfov),
      neardist(neardist),
      fardist(fardist),
      value(0),
      data(NULL)
{
}



R3Camera::
R3Camera(const R3Point& origin, RNAngle pitch, RNAngle yaw, RNAngle roll,
		  RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist)
    : xfov(xfov),
      yfov(yfov),
      neardist(neardist),
      fardist(fardist),
      value(0),
      data(NULL)
{
    // Clamp yaw (roll is undefined at +-PI/2)
    const RNAngle max_yaw = 0.99 * RN_PI_OVER_TWO;
    if (yaw > max_yaw) yaw = max_yaw;
    if (yaw < -max_yaw) yaw = -max_yaw;

    // Derive  triad from pitch, yaw, and roll
    RNScalar cos_yaw = cos(yaw);
    R3Vector towards = R3Vector(cos_yaw * cos(pitch), cos_yaw * sin(pitch), sin(yaw));
    R3Vector up = R3posz_vector;
    up.Rotate(towards, roll);
    R3Triad triad(towards, up);

    // Set coordinate system
    cs = R3CoordSystem(origin, triad);
}



const RNAngle R3Camera::
Pitch(void) const
{
    // Return pitch angle
    R3Vector v1(Towards());
    RNAngle yaw = Yaw();
    if (RNIsNotZero(yaw)) v1.Rotate(Left(), yaw);
    RNAngle pitch = R3InteriorAngle(R3default_camera_towards_vector, v1);
    if (v1.Dot(R3default_camera_right_vector) > 0.0) pitch = -pitch;
    return pitch;
}



const RNAngle R3Camera::
Yaw(void) const
{
    // Return yaw angle
    RNAngle angle = R3InteriorAngle(R3default_camera_up_vector, Towards()); 
    return RN_PI_OVER_TWO - angle;
}



const RNAngle R3Camera::
Roll(void) const
{
#if FALSE
    // Return roll angle
    R3Vector v2(Right());
    RNAngle yaw = Yaw();
    if (RNIsNotZero(yaw)) v2.Rotate(Left(), yaw);
    RNAngle roll = R3InteriorAngle(R3default_camera_up_vector, v2) - RN_PI_OVER_TWO;
    return roll;
#else
    // Return roll angle
    R3Vector v2 = Towards() % R3default_camera_up_vector; // Right vector if no roll
    RNAngle roll = R3InteriorAngle(Right(), v2);
    if (v2.Dot(Up()) > 0.0) roll = -roll;
    return roll;
#endif
}



const R3Halfspace R3Camera::
Halfspace(RNDirection dir, RNDimension dim) const
{
    R3Vector normal;
    switch (dim) {
    case RN_X:
        normal = cs.Axes().Axis(RN_X); 
	if (dir == RN_LO) normal.Rotate(cs.Axes().Axis(RN_Y), xfov);
	else normal.Rotate(cs.Axes().Axis(RN_Y), RN_PI - xfov);
	return R3Halfspace(Origin(), normal);
        break;

    case RN_Y:
        normal = cs.Axes().Axis(RN_Y); 
	if (dir == RN_LO) normal.Rotate(cs.Axes().Axis(RN_X), -yfov);
	else normal.Rotate(cs.Axes().Axis(RN_X), yfov - RN_PI);
	return R3Halfspace(Origin(), normal);
        break;

    case RN_Z:
        if (dir == RN_LO) return R3Halfspace(Origin() - cs.Axes().Axis(RN_Z) * neardist, -cs.Axes().Axis(RN_Z));
	else return R3Halfspace(Origin() - cs.Axes().Axis(RN_Z) * fardist, cs.Axes().Axis(RN_Z));
        break;
    }

    // Should never get here
    RNAbort("Invalid dim in camera halfspace query");
    return R3null_halfspace;
}



const R4Matrix R3Camera::
Perspective(void) const
{
#if FALSE
    // Return 4x4 perspective matrix
    double m[4][4];
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPerspective(2.0 * RN_RAD2DEG(yfov), tan(xfov) / tan(yfov), neardist, fardist);
    glGetDoublev(GL_PROJECTION_MATRIX, m);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    return R4Matrix(m);
#else 
    // Return 4x4 perspective matrix -- derived from OpenGL Programming Guide
    RNScalar A = 1.0 / tan(xfov);
    RNScalar B = 1.0 / tan(yfov);
    RNScalar C = -(fardist + neardist) / (fardist - neardist);
    RNScalar D = -2.0 * neardist * fardist / (fardist - neardist);
    return R4Matrix(  A, 0.0,  0.0, 0.0,
		    0.0,   B,  0.0, 0.0,
		    0.0, 0.0,    C,   D,
		    0.0, 0.0, -1.0, 0.0);
#endif
}



const RNBoolean R3Camera::
operator==(const R3Camera& camera) const
{
    // Return whether camera is equal
    return (cs == camera.cs) && 
           (xfov == camera.xfov) &&
           (yfov == camera.yfov) &&
           (neardist == camera.neardist) &&
           (fardist == camera.fardist);
}



void R3Camera::
Reset(const R3Point& origin, const R3Vector& towards, const R3Vector& up)
{
    // Set coordinate system
    cs.SetOrigin(origin);
    cs.SetAxes(R3Triad(towards, up));
}



void R3Camera::
Reposition(const R3Point& origin)
{
    // Set position
    cs.SetOrigin(origin);
}



void R3Camera::
Reorient(const R3Vector& towards, const R3Vector& up)
{
    // Set orientation
    cs.SetAxes(R3Triad(towards, up));
}



void R3Camera::
Reorient(RNAngle pitch, RNAngle yaw, RNAngle roll)
{
    // Set orientation
    SetPitch(pitch);
    SetYaw(yaw);
    SetRoll(roll);
}



void R3Camera::
Mirror(const R3Plane& plane)
{
    // Mirror coordinate system over plane
    cs.Mirror(plane);
}



void R3Camera::
SetCoordSystem(const R3CoordSystem& cs)
{
    // Set coordinate system
    this->cs.SetOrigin(cs.Origin());;
    this->cs.SetAxes(cs.Axes());;
}



void R3Camera::
SetOrigin(const R3Point& origin)
{
    // Set origin 
    cs.SetOrigin(origin);
}



void R3Camera::
SetTowards(const R3Vector& towards)
{
    // Set towards vector
    cs.SetAxes(R3Triad(towards, cs.Axes().Axis(RN_Y)));
}



void R3Camera::
SetUp(const R3Vector& up)
{
    // Set up vector
    cs.SetAxes(R3Triad(Towards(), up));
}



void R3Camera::
SetRight(const R3Vector& right)
{
    // Set right vector
    cs.SetAxes(R3Triad(Towards(), right % Towards()));
}



void R3Camera::
SetPitch(RNAngle pitch)
{
    // Set pitch angle
    RotatePitch(pitch - Pitch());
}



void R3Camera::
SetYaw(RNAngle yaw)
{
    // Set yaw angle
    RotateYaw(yaw - Yaw());
}



void R3Camera::
SetRoll(RNAngle roll)
{
    // Set roll angle
    RotateRoll(roll - Roll());
}



void R3Camera::
SetFOV(RNAngle xfov, RNAngle yfov)
{
    // Set field of view
    this->xfov = xfov;
    this->yfov = yfov;
}



void R3Camera::
SetXFOV(RNAngle xfov)
{
    // Set x field of view
    this->xfov = xfov;
}



void R3Camera::
SetYFOV(RNAngle yfov)
{
    // Set y field of view
    this->yfov = yfov;
}



void R3Camera::
SetNear(RNLength neardist)
{
    // Set near plane distance
    this->neardist = neardist;
}



void R3Camera::
SetFar(RNLength fardist)
{
    // Set far plane distance
    this->fardist = fardist;
}



void R3Camera::
Translate(const R3Vector& translation)
{
    // Translate coordinate system
    cs.Translate(translation);
}



void R3Camera::
Rotate(const R3Vector& axis, RNAngle dtheta)
{
    // Rotate coordinate system
    cs.Rotate(axis, dtheta);
}



void R3Camera::
Transform(const R3Transformation& transformation)
{
    // Transform near and far distances
    R3Vector zaxis(cs.Axes().Axis(RN_Z));
    zaxis.Transform(transformation);
    RNLength scale = zaxis.Length();
    neardist *= scale;
    fardist *= scale;
    
    // Transform coordinate system
    cs.Transform(transformation);
}



R3Camera& R3Camera::
operator=(const R3Camera& camera)
{
    // Copy camera 
    this->cs = camera.cs;
    this->xfov = camera.xfov;
    this->yfov = camera.yfov;
    this->neardist = camera.neardist;
    this->fardist = camera.fardist;
    this->value = camera.value;
    return *this;
}



void R3Camera::
Load(RNBoolean select_mode) const
{
    // Set the perspective projection
#if (RN_3D_GRFX == RN_IRISGL)
    perspective(20.0 * RN_RAD2DEG(yfov), tan(xfov) / tan(yfov), neardist, fardist);
#elif (RN_3D_GRFX == RN_OPENGL)
    if (!select_mode) {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
    }
    gluPerspective(2.0 * RN_RAD2DEG(yfov), tan(xfov) / tan(yfov), neardist, fardist);
    glMatrixMode(GL_MODELVIEW);
#elif ((RN_3D_GRFX == RN_3DR) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
    PointFW_t loc;
    PointF_t dir, up;
    loc.x = Origin().X();   loc.y = Origin().Y();   loc.z = Origin().Z();   loc.w = 1.0;
    dir.x = Towards().X();  dir.y = Towards().Y();  dir.z = Towards().Z();
    up.x = Up().X();        up.y = Up().Y();        up.z = Up().Z();
    G3dSetCameraPosition(R3dr_gc, &loc, &dir, &up);
#else
    RNAbort("Not Implemented");
#endif

    // Load the rotation matrix
    cs.InverseMatrix().Load();
}



void R3Camera::
Draw(void) const
{
#if TRUE
    // Draw pyramid to near plane
    // static float pyramid_color[3] = { 1.0, 0.0, 0.0 };
    // R3LoadRgb(pyramid_color);
    RNScalar scale = 1.0;
    R3Point org = Origin() + Towards() * scale;
    R3Vector dx = Right() * scale * tan(xfov);
    R3Vector dy = Up() * scale * tan(yfov);
    R3Point ur = org + dx + dy;
    R3Point lr = org + dx - dy;
    R3Point ul = org - dx + dy;
    R3Point ll = org - dx - dy;
    R3BeginLine();
    R3LoadPoint(ur);
    R3LoadPoint(ul);
    R3LoadPoint(ll);
    R3LoadPoint(lr);
    R3LoadPoint(ur);
    R3LoadPoint(Origin());
    R3LoadPoint(lr);
    R3LoadPoint(ll);
    R3LoadPoint(Origin());
    R3LoadPoint(ul);
    R3EndLine();
#endif

#if FALSE
    // Draw coordinate system 
    // static float cs_color[3] = { 0.0, 1.0, 0.0 };
    // R3LoadRgb(cs_color);
    cs.Draw();
#endif
}


int R3CompareCameras(const void *data1, const void *data2)
{
  // Sort the cameras according to highest value first
  R3Camera *camera1 = *((R3Camera **) data1);
  R3Camera *camera2 = *((R3Camera **) data2);
  if (camera1->Value() > camera2->Value()) return -1;
  else if (camera1->Value() < camera2->Value()) return 1;
  return 0;
}
