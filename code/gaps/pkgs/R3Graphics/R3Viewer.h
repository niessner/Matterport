/* Include file for R3 viewer class */



/* Initialization functions */

int R3InitViewer();
void R3StopViewer();



/* Class definition */

class R3Viewer {
    public:
        // Constructor functions
	R3Viewer(void);
        R3Viewer(const R3Viewer& viewer);
        R3Viewer(const R3Camera& camera, const R2Viewport& viewport);
        virtual ~R3Viewer(void);

        // Viewer property functions/operators
        const R3Camera& Camera(void) const;
        const R2Viewport& Viewport(void) const;

        // Camera manipulation functions/operators
        void SetCamera(const R3Camera& camera);
        void RepositionCamera(const R3Point& eyepoint);
        void TranslateCamera(const R3Vector& translation);
        void ReorientCamera(const R3Vector& towards, const R3Vector& up);
        void RotateCamera(const R3Vector& axis, RNAngle dtheta);
        void RotateCameraPitch(RNAngle dpitch);
        void RotateCameraYaw(RNAngle dyaw);
        void RotateCameraRoll(RNAngle droll);
        void MirrorCamera(const R3Plane& plane);
	void ResetCamera(const R3Point& eyepoint, const R3Vector& towards, const R3Vector& up);

	// Viewport manipulation functions/operators
        void SetViewport(const R2Viewport& viewport);
	void MoveViewport(int xmin, int ymin);
	void ResizeViewport(int width, int height);
        void ResizeViewport(int xmin, int ymin, int width, int height);

        // World manipulation functions/operators
        void TranslateWorld(const R3Vector& translation);
        void ScaleWorld(const R3Point& origin, RNScalar factor);
        void RotateWorld(const R3Point& origin, const R3Vector& normal, RNAngle dtheta);
        void RotateWorldPitch(const R3Point& origin, RNAngle dpitch);
        void RotateWorldYaw(const R3Point& origin, RNAngle dyaw);
        void RotateWorldRoll(const R3Point& origin, RNAngle droll);

        // Camera-in-hand user interface utility functions
	virtual void TranslateCamera(RNScalar factor, int x, int y, int dx, int dy);
	virtual void RotateCamera(RNScalar factor, int x, int y, int dx, int dy);
	virtual void RotateCameraPitch(RNScalar factor, int x, int y, int dx, int dy);
	virtual void RotateCameraYaw(RNScalar factor, int x, int y, int dx, int dy);
	virtual void RotateCameraRoll(RNScalar factor, int x, int y, int dx, int dy);
	virtual void FlyCamera(RNScalar pitch_factor, RNScalar yaw_factor, RNScalar roll_factor, 
	    RNScalar step_factor, RNBoolean relax_yaw, RNBoolean relax_roll, int x, int y, int dx, int dy);
	virtual void WalkCamera(RNScalar pitch_factor, RNScalar yaw_factor, RNScalar roll_factor, 
	    RNScalar step_factor, RNBoolean relax_yaw, RNBoolean relax_roll, int x, int y, int dx, int dy);

        // World-in-hand user interface utility functions
	virtual void TranslateWorld(RNScalar factor, const R3Point& origin, int x, int y, int dx, int dy);
	virtual void ScaleWorld(RNScalar factor, const R3Point& origin, int x, int y, int dx, int dy);
	virtual void RotateWorld(RNScalar factor, const R3Point& origin, int x, int y, int dx, int dy);
	virtual void RotateWorldPitch(RNScalar factor, const R3Point& origin, int x, int y, int dx, int dy);
	virtual void RotateWorldYaw(RNScalar factor, const R3Point& origin, int x, int y, int dx, int dy);
	virtual void RotateWorldRoll(RNScalar factor, const R3Point& origin, int x, int y, int dx, int dy);

	// Camera/viewport relationship functions/operators
	const R3Ray WorldRay(int x, int y) const;
	const R2Point ViewportPoint(const R3Point& point) const;
	const R2Box ViewportBBox(const R3Shape& shape) const;
	const R2Circle ViewportBCircle(const R3Shape& shape) const;
	const RNArea ViewportArea(const R3Shape& shape) const;

        // Comparison operators
	const RNBoolean operator==(const R3Viewer& viewer) const;
	const RNBoolean operator!=(const R3Viewer& viewer) const;

        // Draw functions
        void Load(void) const;

    private:
        R3Camera camera;
        R2Viewport viewport;
};




/* Inline functions */


inline const R3Camera& R3Viewer::
Camera(void) const
{
    // Return viewer's camera
    return camera;
}



inline const R2Viewport& R3Viewer::
Viewport(void) const
{
    // Return viewer's viewport
    return viewport;
}



inline void R3Viewer::
RepositionCamera (const R3Point& eyepoint)
{
    // Reposition camera origin at point
    camera.Reposition(eyepoint);
}



inline void R3Viewer::
TranslateCamera (const R3Vector& translation)
{
    // Translate camera by translation
    camera.Translate(translation);
}



inline void R3Viewer::
ReorientCamera(const R3Vector& towards, const R3Vector& up)
{
    // Reorient camera
    camera.Reorient(towards, up);
}



inline void R3Viewer::
RotateCamera (const R3Vector& axis, RNAngle dtheta)
{
    // Rotate camera 
    camera.Rotate(axis, dtheta);
}



inline void R3Viewer::
RotateCameraPitch (RNAngle dpitch)
{
    // Rotate camera 
    camera.RotatePitch(dpitch);
}



inline void R3Viewer::
RotateCameraYaw (RNAngle dyaw)
{
    // Rotate camera 
    camera.RotateYaw(dyaw);
}



inline void R3Viewer::
RotateCameraRoll (RNAngle droll)
{
    // Rotate camera 
    camera.RotateRoll(droll);
}



inline void R3Viewer::
MirrorCamera (const R3Plane& plane)
{
    // Mirror camera 
    camera.Mirror(plane);
}



inline void R3Viewer::
ResetCamera(const R3Point& eyepoint, const R3Vector& towards, const R3Vector& up)
{
    // Rotate camera 
    camera.Reset(eyepoint, towards, up);
}



inline void R3Viewer::
RotateWorldPitch (const R3Point& origin, RNAngle dpitch)
{
    // Rotate around up axis
    RotateWorld(origin, camera.Up(), dpitch);
}



inline void R3Viewer::
RotateWorldYaw (const R3Point& origin, RNAngle dyaw)
{
    // Rotate around left axis
    RotateWorld(origin, camera.Left(), dyaw);
}



inline void R3Viewer::
RotateWorldRoll (const R3Point& origin, RNAngle droll)
{
    // Rotate around towards axis
    RotateWorld(origin, camera.Towards(), droll);
}



inline void R3Viewer::
MoveViewport (int xmin, int ymin) 
{
    // Move viewport
    viewport.Move(xmin, ymin);
}



inline const RNBoolean R3Viewer::
operator==(const R3Viewer& viewer) const
{
    // Return whether viewer is not equal
    return (camera == viewer.camera) && 
           (viewport == viewer.viewport); 
}



inline const RNBoolean R3Viewer::
operator!=(const R3Viewer& viewer) const
{
    // Return whether viewer is not equal
    return !(*this == viewer);
}



inline void R3Viewer::
Load(void) const
{
    // Load viewport and camera
    viewport.Load();
    camera.Load();
}



