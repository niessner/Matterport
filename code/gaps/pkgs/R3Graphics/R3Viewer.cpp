/* Source file for the R3 viewer class */



/* Include files */

#include "R3Graphics.h"



/* Public functions */

int 
R3InitViewer()
{
    /* Return success */
    return TRUE;
}



void 
R3StopViewer()
{
}



R3Viewer::
R3Viewer(void)
{
}



R3Viewer::
R3Viewer(const R3Viewer& viewer)
    : camera(viewer.camera),
      viewport(viewer.viewport)
{
}



R3Viewer::
R3Viewer(const R3Camera& camera, const R2Viewport& viewport)
    : camera(camera),
      viewport(viewport)
{
    // Adjust camera field of view to match viewport aspect ratio
    RNScalar aspect = ((RNScalar) viewport.Height()) / ((RNScalar) viewport.Width());
    this->camera.SetYFOV(atan(aspect * tan(camera.XFOV())));
}



R3Viewer::
~R3Viewer(void)
{
}



void R3Viewer::
SetCamera(const R3Camera& camera)
{
    // Set viewer's camera
    this->camera = camera;
    RNScalar aspect = ((RNScalar) viewport.Height()) / ((RNScalar) viewport.Width());
    this->camera.SetYFOV(atan(aspect * tan(camera.XFOV())));
}



void R3Viewer::
SetViewport(const R2Viewport& viewport)
{
    // Set viewer's viewport
    this->viewport = viewport;

    // Adjust camera field of view to match viewport aspect ratio
    RNScalar aspect = (RNScalar) viewport.Height() / (RNScalar) viewport.Width();
    camera.SetYFOV(atan(aspect * tan(camera.XFOV())));
}



void R3Viewer::
ResizeViewport (int width, int height) 
{
    // Resize viewport
    viewport.Resize(width, height);

    // Adjust camera field of view to match viewport aspect ratio
    RNScalar aspect = (RNScalar) height / (RNScalar) width;
    camera.SetYFOV(atan(aspect * tan(camera.XFOV())));
}



void R3Viewer::
ResizeViewport (int xmin, int ymin, int width, int height) 
{
    // Resize viewport
    viewport.Resize(xmin, ymin, width, height);

    // Adjust camera field of view to match viewport aspect ratio
    RNScalar aspect = (RNScalar) height / (RNScalar) width;
    camera.SetYFOV(atan(aspect * tan(camera.XFOV())));
}



void R3Viewer::
TranslateWorld (const R3Vector& translation)
{
    // Translate camera by -translation
    camera.Translate(-translation);
}



void R3Viewer::
ScaleWorld (const R3Point& origin, RNScalar factor)
{
    // Scale world by factor - translate camera
    assert(factor > 0.0);
    factor = (factor - 1.0) / factor;
    camera.Translate((origin - camera.Origin()) * factor);
}



void R3Viewer::
RotateWorld (const R3Point& origin, const R3Vector& axis, RNAngle dtheta)
{
    // Rotate camera around line (origin, axis)
    R3Line axis_line = R3Line(origin, axis);
    R3Point apex = camera.Origin();
    apex.Rotate(axis_line, dtheta);
    R3Vector towards = camera.Towards();
    towards.Rotate(axis, dtheta);
    R3Vector up = camera.Up();
    up.Rotate(axis, dtheta);
    camera.Reset(apex, towards, up);
}



void R3Viewer::
TranslateWorld(RNScalar factor, const R3Point& origin, int, int, int dx, int dy)
{
    // Translate world based on mouse (dx,dy)
    if ((dx == 0) && (dy == 0)) return;
    RNLength length = factor * R3Distance(origin, camera.Origin()) * tan(camera.XFOV());
    RNLength vx = length * (RNLength) dx / (RNLength) viewport.Width();
    RNLength vy = length * (RNLength) dy / (RNLength) viewport.Height();
    R3Vector vector = (camera.Right() * vx) + (camera.Up() * vy);
    TranslateWorld(vector);
}



void R3Viewer::
ScaleWorld(RNScalar factor, const R3Point& origin, int, int, int dx, int dy)
{
    // Scale world based on mouse dx
    if ((dx == 0) && (dy == 0)) return;
    RNScalar motion = (RNScalar) dx / (RNScalar) viewport.Width() + (RNScalar) dy / (RNScalar) viewport.Height();
    RNScalar scale = exp(2.0 * factor * motion);
    ScaleWorld(origin, scale);
}



void R3Viewer::
RotateWorld(RNScalar factor, const R3Point& origin, int, int, int dx, int dy)
{
    // Rotate world based on mouse (dx)
    if ((dx == 0) && (dy == 0)) return;
    RNLength vx = (RNLength) dx / (RNLength) viewport.Width();
    RNLength vy = (RNLength) dy / (RNLength) viewport.Height();
    RNAngle theta = factor * 4.0 * (fabs(vx) + fabs(vy));
    R3Vector vector = (camera.Right() * vx) + (camera.Up() * vy);
    R3Vector axis = camera.Towards() % vector;
    axis.Normalize();
    RotateWorld(origin, axis, theta);
}



void R3Viewer::
RotateWorldPitch(RNScalar factor, const R3Point& origin, int, int, int dx, int)
{
    // Rotate world based on mouse (dx)
    RNLength vx = factor * (RNLength) dx / (RNLength) viewport.Width();
    RotateWorldPitch(origin, 4.0 * vx);
}



void R3Viewer::
RotateWorldYaw(RNScalar factor, const R3Point& origin, int, int, int, int dy)
{
    // Rotate world based on mouse (dy)
    RNLength vy = factor * (RNLength) dy / (RNLength) viewport.Height();
    RotateWorldYaw(origin, 4.0 * vy);
}



void R3Viewer::
RotateWorldRoll(RNScalar factor, const R3Point& origin, int, int, int dx, int)
{
    // Rotate world based on mouse (dx)
    RNLength vx = factor * (RNLength) dx / (RNLength) viewport.Width();
    RotateWorldRoll(origin, 4.0 * vx);
}



void R3Viewer::
TranslateCamera(RNScalar factor, int, int, int, int)
{
    // Translate camera in direction of view
    TranslateCamera(camera.Towards() * factor);
}



void R3Viewer::
RotateCamera(RNScalar factor, int x, int y, int, int)
{
    // Rotate camera based on mouse (x,y)                                                                                           
    RNLength vx = (RNLength) ((2 * x) - viewport.Width()) / (RNLength) viewport.Width();
    RNLength vy = (RNLength) ((2 * y) - viewport.Height()) / (RNLength) viewport.Height();
    RNAngle theta = factor * (RN_PI/64.0) * (vx*vx + vy*vy);
    R3Vector vector = (camera.Right() * vx) + (camera.Up() * vy);
    R3Vector axis = camera.Towards() % vector;
    axis.Normalize();
    RotateCamera(axis, theta);
}



void R3Viewer::
RotateCameraPitch(RNScalar factor, int x, int, int, int)
{
    // Rotate camera based on mouse (x)
    RNLength vx = (RNLength) ((2 * x) - viewport.Width()) / (RNLength) viewport.Width();
    RNAngle dpitch = factor * (RN_PI/64.0) * vx * ((vx < 0.0) ? -vx : vx); // vx*vx*vx;
    camera.RotatePitch(-dpitch);
}



void R3Viewer::
RotateCameraYaw(RNScalar factor, int, int y, int, int)
{
    // Rotate camera based on mouse (x,y)
    RNLength vy = (RNLength) ((2 * y) - viewport.Height()) / (RNLength) viewport.Height();
    RNAngle dyaw = factor * (RN_PI/64.0) * vy * ((vy < 0.0) ? -vy : vy); // vy*vy*vy;
    camera.RotateYaw(dyaw);
}



void R3Viewer::
RotateCameraRoll(RNScalar factor, int x, int, int, int)
{
    // Rotate camera based on mouse (x,y)
    RNLength vx = (RNLength) ((2 * x) - viewport.Width()) / (RNLength) viewport.Width();
    RNAngle droll = factor * (RN_PI/64.0) * vx * ((vx < 0.0) ? -vx : vx); // vx*vx*vx;
    camera.RotateRoll(droll);
}



void R3Viewer::
FlyCamera(RNScalar pitch_factor, RNScalar yaw_factor, RNScalar roll_factor, RNScalar step_factor, 
    RNBoolean relax_yaw, RNBoolean relax_roll, int x, int y, int dx, int dy)
{
    // Update pitch
    if (RNIsNotZero(pitch_factor)) {
	RotateCameraPitch(pitch_factor, x, y, dx, dy);
    }

    // Update yaw
    if (RNIsNotZero(yaw_factor)) {
	RotateCameraYaw(yaw_factor, x, y, dx, dy);
    }
    else if (relax_yaw && (RNIsNotZero(camera.Yaw()))) {
	if (RNIsZero(camera.Yaw(), 0.005)) RotateCameraYaw(-1.0 * camera.Yaw());
	else RotateCameraYaw(-0.05 * camera.Yaw());
    }
	
    // Update roll
    if (RNIsNotZero(roll_factor)) {
	RotateCameraRoll(roll_factor, x, y, dx, dy);
    }
    else if (relax_roll && (RNIsNotZero(camera.Roll()))) {
	if (camera.Roll() < 0.01) RotateCameraRoll(-1.0 * camera.Roll());
	else RotateCameraRoll(-0.05 * camera.Roll());
    }
	
    // Translate camera in direction of view 
    if (RNIsNotZero(step_factor)) {
	RNLength vx = (RNLength) ((2 * x) - viewport.Width()) / (RNLength) viewport.Width();
	if (vx > 1.0) vx = 1.0;
	else if (vx < -1.0) vx = -1.0;
	step_factor *= 1.0 - vx*vx;
	TranslateCamera(camera.Towards() * step_factor);
    }
}



void R3Viewer::
WalkCamera(RNScalar pitch_factor, RNScalar yaw_factor, RNScalar roll_factor, RNScalar step_factor, 
    RNBoolean relax_yaw, RNBoolean relax_roll, int x, int y, int dx, int dy)
{
    // Update pitch
    if (RNIsNotZero(pitch_factor)) {
	RotateCameraPitch(pitch_factor, x, y, dx, dy);
    }

    // Update yaw
    if (RNIsNotZero(yaw_factor)) {
	RotateCameraYaw(yaw_factor, x, y, dx, dy);
    }
    else if (relax_yaw && (RNIsNotZero(camera.Yaw()))) {
	if (RNIsZero(camera.Yaw(), 0.005)) RotateCameraYaw(-1.0 * camera.Yaw());
	else RotateCameraYaw(-0.05 * camera.Yaw());
    }
	
    // Update roll
    if (RNIsNotZero(roll_factor)) {
	RotateCameraRoll(roll_factor, x, y, dx, dy);
    }
    else if (relax_roll && (RNIsNotZero(camera.Roll()))) {
	if (camera.Roll() < 0.01) RotateCameraRoll(-1.0 * camera.Roll());
	else RotateCameraRoll(-0.05 * camera.Roll());
    }
	
    // Translate camera in direction of view ALONG HORIZONTAL PLANE
    if (RNIsNotZero(step_factor)) {
	RNLength vx = (RNLength) ((2 * x) - viewport.Width()) / (RNLength) viewport.Width();
	if (vx > 1.0) vx = 1.0;
	else if (vx < -1.0) vx = -1.0;
	step_factor *= 1.0 - vx*vx;
	R3Vector v = camera.Towards();
	v.Project(R3Plane(camera.Origin(), R3default_camera_up_vector));
	TranslateCamera(v * step_factor);
    }
}



const R3Ray R3Viewer::
WorldRay(int x, int y) const
{
    // Return ray from camera origin to appropriate point on far plane
#if FALSE
    // This also works
    RNScalar dx = (RNScalar) (2 * (x - viewport.XCenter())) / (RNScalar) viewport.Width();
    RNScalar dy = (RNScalar) (2 * (y - viewport.YCenter())) / (RNScalar) viewport.Height();
    R3Vector v = camera.Towards();
    v.Rotate(camera.Up(), dx * camera.XFOV());
    v.Rotate(camera.Left(), dy * camera.YFOV());
    return R3Ray(camera.Origin(), v);
#else
    RNScalar dx = (RNScalar) (2 * (x - viewport.XCenter())) / (RNScalar) viewport.Width();
    RNScalar dy = (RNScalar) (2 * (y - viewport.YCenter())) / (RNScalar) viewport.Height();
    R3Point far_org = camera.Origin() + camera.Towards() * camera.Far();
    R3Vector far_right = camera.Right() * camera.Far() * tan(camera.XFOV());
    R3Vector far_up = camera.Up() * camera.Far() * tan(camera.YFOV());
    R3Point far_point = far_org + (far_right * dx) + (far_up * dy);
    return R3Ray(camera.Origin(), far_point);
#endif
}



const R2Point R3Viewer::
ViewportPoint(const R3Point& point) const 
{
    // Transform 3D point into canonical coordinate system
    R3Point p = Camera().CoordSystem().InverseMatrix() * point;
    if (RNIsPositiveOrZero(p.Z())) return R2infinite_point;

    // Compute 2D point projected onto viewport
    RNCoord x = viewport.XCenter() + 0.5 * viewport.Width() * p.X() / -(p.Z() * tan(Camera().XFOV()));
    RNCoord y = viewport.YCenter() + 0.5 * viewport.Height() * p.Y() / -(p.Z() * tan(Camera().YFOV()));

    // Return point projected onto viewport
    return R2Point(x, y);
}



const R2Box R3Viewer::
ViewportBBox(const R3Shape& shape) const 
{
    // ??? THIS FUNCTION DOES NOT HANDLE SOME INVISIBLE SHAPES PROPERLY ???
    // ??? IT RETURNS FULL SCREEN BOXES FOR INVISIBLE DIAGONAL BOXES    ???
    // ??? THAT SPAN VIEW ORIGIN IN EVERY DIMENSION                     ???

    // Compute bounding box of shape projection on viewport
    RNBoolean infront = FALSE;
    R3Box shape_box(shape.BBox());
    R2Box projection_box(R2null_box);
    R4Matrix m(Camera().CoordSystem().InverseMatrix());
    for (RNOctant octant = RN_NNN_OCTANT; octant <= RN_PPP_OCTANT; octant++) {
        // Transform corner into canonical coordinate system
        R3Point p = m * shape_box.Corner(octant);

	// Union point's projection into projection box
	if (RNIsPositiveOrZero(p.Z())) {
	    // Point is behind viewer - union extrema into projection box
	    if (RNIsNegative(p.X())) projection_box[RN_LO][RN_X] = -RN_INFINITY;
	    else if (RNIsPositive(p.X())) projection_box[RN_HI][RN_X] = RN_INFINITY;
	    if (RNIsNegative(p.Y())) projection_box[RN_LO][RN_Y] = -RN_INFINITY;
	    else if (RNIsPositive(p.Y())) projection_box[RN_HI][RN_Y] = RN_INFINITY;
	}
	else {
	    // Point is in front of viewer - Union projection point into projection box
	    RNCoord x = viewport.XCenter() + 0.5 * viewport.Width() * p.X() / -(p.Z() * tan(Camera().XFOV()));
	    RNCoord y = viewport.YCenter() + 0.5 * viewport.Height() * p.Y() / -(p.Z() * tan(Camera().YFOV()));
	    projection_box.Union(R2Point(x, y));
	    infront = TRUE;
	}
    }

    // Check if some part of shape box is in front of viewer
    if (!infront) return R2null_box;

    // Intersect projection box with viewport box
    R2Box viewport_box(Viewport().XMin(), Viewport().YMin(), Viewport().XMax(), Viewport().YMax());
    projection_box.Intersect(viewport_box);

    // Return projection box
    return projection_box;
}



const R2Circle R3Viewer::
ViewportBCircle(const R3Shape& shape) const 
{
    // Return bounding circle of shape projection on viewport
    R3Sphere bsphere(shape.BSphere());
    R3Point p = Camera().CoordSystem().InverseMatrix() * bsphere.Center();
    if (RNIsPositiveOrZero(p.Z())) return R2infinite_circle;

    // Compute sphere center projection
    RNCoord x = viewport.XCenter() + 0.5 * viewport.Width() * p.X() / -(p.Z() * tan(Camera().XFOV()));
    RNCoord y = viewport.YCenter() + 0.5 * viewport.Height() * p.Y() / -(p.Z() * tan(Camera().YFOV()));
    R2Point center(x, y);

    // Compute sphere radius projection
    RNScalar fov = (Camera().XFOV() > Camera().YFOV()) ? Camera().XFOV() : Camera().YFOV();
    RNScalar radius = 0.5 * viewport.Width() * bsphere.Radius() / -(p.Z() * fov);

    // Return bounding circle of projection
    return R2Circle(center, radius);
}



const RNArea R3Viewer::
ViewportArea(const R3Shape& shape) const 
{
    // Return area of shape projection on viewport
    R2Circle bcircle(ViewportBCircle(shape));
    if (bcircle.IsFinite()) return bcircle.Area();
    else return 0.0;
}















