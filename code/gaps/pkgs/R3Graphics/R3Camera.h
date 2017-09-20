/* Include file for R3 camera class */



/* Initialization functions */

int R3InitCamera();
void R3StopCamera();



/* Class definition */

class R3Camera {
    public:
        // Constructor functions
	R3Camera(void);
        R3Camera(const R3Camera& camera);
	R3Camera(const R3CoordSystem &cs, RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist);
	R3Camera(const R3Point& origin, const R3Triad& triad, RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist);
	R3Camera(const R3Point& origin, const R3Vector& towards, const R3Vector& up, RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist);
	R3Camera(const R3Point& origin, RNAngle pitch, RNAngle yaw, RNAngle roll, RNAngle xfov, RNAngle yfov, RNLength neardist, RNLength fardist);

	// Property functions/operators
        const R3Point& Origin(void) const;        
        const R3Vector Towards(void) const;        
        const R3Vector& Backwards(void) const;        
        const R3Vector& Up(void) const;        
        const R3Vector Down(void) const;        
        const R3Vector& Right(void) const;        
        const R3Vector Left(void) const;        
	const RNAngle Pitch(void) const;
	const RNAngle Yaw(void) const;
	const RNAngle Roll(void) const;
        const RNAngle XFOV(void) const;
        const RNAngle YFOV(void) const;
        const RNLength Near(void) const;
        const RNLength Far(void) const;
	const R3Triad& Triad(void) const;
	const R3CoordSystem& CoordSystem(void) const;
	const R3Halfspace Halfspace(RNDirection dir, RNDimension dim) const;
	const R4Matrix Perspective(void) const;
        const RNScalar Value(void) const;
        const void *Data(void) const;
	const RNBoolean operator==(const R3Camera& camera) const;
	const RNBoolean operator!=(const R3Camera& camera) const;

        // Manipulation functions/operators
        void Reset(const R3Point& origin, const R3Vector& towards, const R3Vector& up);
        void Reposition(const R3Point& origin);
        void Reorient(const R3Vector& towards, const R3Vector& up);
        void Reorient(RNAngle pitch, RNAngle yaw, RNAngle roll);
        void Mirror(const R3Plane& plane);
        void SetCoordSystem(const R3CoordSystem& cs);        
        void SetOrigin(const R3Point& origin);        
        void SetTowards(const R3Vector& towards);
        void SetBackwards(const R3Vector& backwards);
	void SetUp(const R3Vector& up);
        void SetDown(const R3Vector& down);
        void SetRight(const R3Vector& right);
        void SetLeft(const R3Vector& left);
        void SetPitch(RNAngle pitch);
        void SetYaw(RNAngle pitch);
        void SetRoll(RNAngle pitch);
        void SetFOV(RNAngle xfov, RNAngle yfov);
        void SetXFOV(RNAngle xfov);
        void SetYFOV(RNAngle yfov);
        void SetNear(RNLength neardist);
        void SetFar(RNLength fardist);
        void SetValue(RNScalar value);
        void SetData(void *data);
  
        // More manipulation functions/operators
        void Translate(const R3Vector& translation);
        void Rotate(const R3Vector& axis, RNAngle dtheta);
        void RotatePitch(RNAngle dpitch);
        void RotateYaw(RNAngle dyaw);
        void RotateRoll(RNAngle droll);
	void Transform(const R3Transformation& transformation);
	R3Camera& operator=(const R3Camera& camera);

        // Draw functions/operators
        void Load(RNBoolean select_mode = FALSE) const;
        void Draw(void) const;
        void Outline(void) const;

    private:
        R3CoordSystem cs;
        RNAngle xfov, yfov;
        RNScalar neardist, fardist;
        RNScalar value;
        void *data;
};



/* Public variables */

extern RNLength R3default_camera_xfov;
extern RNLength R3default_camera_yfov;
extern R3Vector R3default_camera_towards_vector;
extern R3Vector R3default_camera_up_vector;
extern R3Vector R3default_camera_right_vector;
extern R3Camera R3default_camera;



/* Public functions */

extern int R3CompareCameras(const void *data1, const void *data2);



/* Inline functions */

inline const R3Point& R3Camera::
Origin(void) const
{
    // Return origin
    return cs.Origin();
}



inline const R3Vector R3Camera::
Towards(void) const
{
    // Return towards vector
    return -(cs.Axes().Axis(RN_Z));
}



inline const R3Vector& R3Camera::
Backwards(void) const
{
    // Return backwards vector
    return cs.Axes().Axis(RN_Z);
}



inline const R3Vector& R3Camera::
Up(void) const
{
    // Return up vector
    return cs.Axes().Axis(RN_Y);
}



inline const R3Vector R3Camera::
Down(void) const
{
    // Return down vector
    return -Up();
}



inline const R3Vector& R3Camera::
Right(void) const
{
    // Return right vector
    return cs.Axes().Axis(RN_X);
}



inline const R3Vector R3Camera::
Left(void) const
{
    // Return left vector
    return -Right();
}



inline const RNAngle R3Camera::
XFOV(void) const
{
    // Return x field of view
    return xfov;
}



inline const RNAngle R3Camera::
YFOV(void) const
{
    // Return y field of view
    return yfov;
}



inline const RNLength R3Camera::
Near(void) const
{
    // Return near plane distance
    return neardist;
}



inline const RNLength R3Camera::
Far(void) const
{
    // Return far plane distance
    return fardist;
}



inline const R3Triad& R3Camera::
Triad(void) const
{
    // Return triad
    return cs.Axes();
}



inline const R3CoordSystem& R3Camera::
CoordSystem(void) const
{
    // Return coordinate system
    return cs;
}



inline const RNScalar R3Camera::
Value(void) const
{
    // Return user-defined value
    return value;
}



inline const void *R3Camera::
Data(void) const
{
    // Return user-defined data
    return data;
}



inline const RNBoolean R3Camera::
operator!=(const R3Camera& camera) const
{
    // Return whether camera is not equal
    return !(*this == camera);
}



inline void R3Camera::
SetBackwards(const R3Vector& backwards)
{
    // Set backwards vector
    SetTowards(-backwards);
}



inline void R3Camera::
SetDown(const R3Vector& down)
{
    // Set down vector
    SetUp(-down);
}



inline void R3Camera::
SetLeft(const R3Vector& left)
{
    // Set left vector
    SetRight(-left);
}



inline void R3Camera::
SetValue(RNScalar value)
{
    // Set user-defined value
    this->value = value;
}



inline void R3Camera::
SetData(void *data)
{
    // Set user-defined data
    this->data = data;
}



inline void R3Camera::
RotatePitch(RNAngle dpitch)
{
    // Rotate around up vector
    // Rotate(Up(), dpitch);
    Rotate(R3default_camera_up_vector, dpitch);
}



inline void R3Camera::
RotateYaw(RNAngle dyaw)
{
    // Rotate around left vector
    Rotate(Right(), dyaw);
}



inline void R3Camera::
RotateRoll(RNAngle droll)
{
    // Rotate around towards vector
    Rotate(Towards(), droll);
}



inline void R3Camera::
Outline(void) const
{
    // Draw camera
    Draw();
}



