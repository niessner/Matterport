/* Include file for the GAPS quaternion class */



/* Initialization functions */

int R3InitQuaternion();
void R3StopQuaternion();



/* Class definition */

class R3Quaternion /* : public R3Base */ {
    public:
        // Constructor functions
	R3Quaternion(void);
        R3Quaternion(const R3Quaternion& quaternion);
	R3Quaternion(RNScalar a, RNScalar b, RNScalar c, RNScalar d);
        R3Quaternion(const R3Quaternion& quaternion1, const R3Quaternion& quaternion2, RNScalar t);
        R3Quaternion(const R3Vector& axis, RNAngle theta);
        R3Quaternion(const R4Matrix& matrix, int dummy);
        R3Quaternion(int dimension, RNAngle theta);
        R3Quaternion(RNAngle pitch, RNAngle yaw, RNAngle roll);
	R3Quaternion(const RNScalar array[4]);
  
        // Property functions/operators
	const RNScalar A(void) const;
	const RNScalar B(void) const;
	const RNScalar C(void) const;
	const RNScalar D(void) const;
	const RNScalar operator[](int i) const;
	const RNBoolean IsZero(void) const;
        const R3Vector Axis(void) const;
        const RNAngle Angle(void) const;
        const R4Matrix Matrix(void) const;
        const R3Quaternion Inverse(void) const;

	// Relationship functions/operators
	const RNBoolean operator==(const R3Quaternion& quaternion) const;
	const RNBoolean operator!=(const R3Quaternion& quaternion) const;

        // Manipulation functions/operators
        void Flip(void);
        void Multiply(const R3Quaternion& quaternion);
	void XRotate(RNAngle radians);
	void YRotate(RNAngle radians);
	void ZRotate(RNAngle radians);
	void Rotate(const R3Vector& xyz_radians);
	void Rotate(int dimension, RNAngle radians);
	void Rotate(const R3Vector& axis, RNAngle radians);
        void Rotate(const R3Quaternion& quaternion);
        void Rotate(const R4Matrix& matrix);
	void Reset(RNScalar a, RNScalar b, RNScalar c, RNScalar d);

	// Assignment operators
	R3Quaternion& operator=(const R3Quaternion& quaternion);
	R3Quaternion& operator*=(const R3Quaternion& quaternion);

        // Arithmetic operators
	friend R3Quaternion operator*(const R3Quaternion& quaternion1, const R3Quaternion& quaternion2);
	friend R3Point operator*(const R3Quaternion& quaternion, const R3Point& point);
	friend R3Vector operator*(const R3Quaternion& quaternion, const R3Vector& vector);
	friend R3Quaternion R3QuaternionSlerp(const R3Quaternion& quaternion1, const R3Quaternion& quaternion2, RNScalar t);

	// Undocumented functions/operators
  	RNScalar& operator[](int i);

    private:
        // Internal functions
        void Normalize(void);
        void Conjugate(void);

    private:
	RNScalar v[4];
};



/* Public variables */

extern const R3Quaternion R3null_quaternion;
extern const R3Quaternion R3zero_quaternion;
extern const R3Quaternion R3identity_quaternion;



/* Inline functions */

inline const RNScalar R3Quaternion::
A(void) const
{
    return (v[0]);
}



inline const RNScalar R3Quaternion::
B(void) const
{
    return (v[1]);
}



inline const RNScalar R3Quaternion::
C(void) const
{
    return (v[2]);
}



inline const RNScalar R3Quaternion::
D(void) const
{
    return (v[3]);
}



inline const RNScalar R3Quaternion::
operator[](int k) const
{
    assert((k>=0) && (k<=3));
    return(v[k]);
}



inline const R3Vector R3Quaternion::
Axis(void) const
{
    // Return axis of rotation
    R3Vector axis(v[1], v[2], v[3]);
    axis.Normalize();
    return axis;
}



inline const RNAngle R3Quaternion::
Angle(void) const
{
    // Return amount of counter-clockwise rotation around axis
    return 2 * acos(v[0]);
}



inline const RNBoolean R3Quaternion::
IsZero (void) const
{
    // Return whether quaternion is zero
    return ((v[0] == 0.0) && (v[1] == 0.0) && (v[2] == 0.0) && (v[3] == 0.0));
}



inline void R3Quaternion::
Reset(RNScalar w, RNScalar x, RNScalar y, RNScalar z) 
{
    // Set all coords
    v[0] = w;
    v[1] = x;
    v[2] = y;
    v[3] = z;
}



inline RNScalar& R3Quaternion::
operator[] (RNDimension k) 
{
    assert((k>=0) && (k<=3));
    return(v[k]);
}



