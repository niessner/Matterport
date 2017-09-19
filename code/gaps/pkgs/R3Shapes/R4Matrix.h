/* Include file for the GAPS matrix class */



/* Initialization functions */

int R3InitMatrix();
void R3StopMatrix();



/* Class definition */

class R4Matrix /* : public R3Base */ {
    public:
        // Constructor functions
	R4Matrix(void);
        R4Matrix(const R4Matrix& matrix);
	R4Matrix(RNScalar a00, RNScalar a01, RNScalar a02, RNScalar a03,
   	         RNScalar a10, RNScalar a11, RNScalar a12, RNScalar a13,
                 RNScalar a20, RNScalar a21, RNScalar a22, RNScalar a23,
                 RNScalar a30, RNScalar a31, RNScalar a32, RNScalar a33);
	R4Matrix(const RNScalar* array);

        // Property functions/operators
	const RNBoolean IsZero(void) const;
	const RNBoolean IsIdentity(void) const;
        const RNBoolean IsIsotropic(void) const;
        const RNBoolean HasTranslation(void) const;
        const RNBoolean HasScale(void) const;
        const RNBoolean HasRotation(void) const;
        const RNBoolean HasMirror(void) const;
  	const RNScalar *operator[](int i) const;
        const RNScalar Determinant(void) const;
	const R4Matrix Transpose(void) const;
	const R4Matrix Inverse(void) const;
	const RNBoolean operator==(const R4Matrix& matrix) const;
	const RNBoolean operator!=(const R4Matrix& matrix) const;

        // Draw functions/operators
        void Load(void) const;
        void Draw(void) const;
        void Push(void) const;
        void Pop(void) const;

        // Manipulation functions/operators
	void Flip(void);
	void Invert(void);
	void XTranslate(RNScalar offset);
	void YTranslate(RNScalar offset);
	void ZTranslate(RNScalar offset);
	void Translate(RNScalar offset);
	void Translate(RNAxis axis, RNScalar offset);
	void Translate(const R3Vector& offset);
	void XScale(RNScalar scale);
	void YScale(RNScalar scale);
	void ZScale(RNScalar scale);
	void Scale(RNScalar scale);
	void Scale(RNAxis axis, RNScalar scale);
	void Scale(const R3Vector& scale);
	void XRotate(RNAngle radians);
	void YRotate(RNAngle radians);
	void ZRotate(RNAngle radians);
	void Rotate(const R3Vector& xyz_radians);
	void Rotate(RNAxis axis, RNAngle radians);
	void Rotate(const R3Vector& vector, RNAngle radians);
	void Rotate(const R3Vector& from, const R3Vector& to);
	void Rotate(const R3Quaternion& quaternion);
	void Transform(const R4Matrix& matrix);
	void Multiply(const R4Matrix& matrix);
	void Add(const R4Matrix& matrix);
	void Subtract(const R4Matrix& matrix);

	// Assignment operators
	R4Matrix& operator=(const R4Matrix& matrix);
	R4Matrix& operator+=(const R4Matrix& matrix);
	R4Matrix& operator-=(const R4Matrix& matrix);
	R4Matrix& operator*=(RNScalar a);
	R4Matrix& operator*=(const R4Matrix& matrix);
	R4Matrix& operator/=(RNScalar a);

        // Arithmetic operators
	friend R4Matrix operator-(const R4Matrix& matrix);
	friend R4Matrix operator+(const R4Matrix& matrix1, const R4Matrix& matrix2);
	friend R4Matrix operator-(const R4Matrix& matrix1, const R4Matrix& matrix2);
	friend R4Matrix operator*(RNScalar a, const R4Matrix& matrix);
	friend R4Matrix operator*(const R4Matrix& matrix, RNScalar a);
	friend R4Matrix operator*(const R4Matrix& matrix1, const R4Matrix& matrix2);
	friend R4Matrix operator/(const R4Matrix& matrix, RNScalar scale);
	friend R3Vector operator*(const R4Matrix& matrix, const R3Vector& vector);
	friend R3Point operator*(const R4Matrix& matrix, const R3Point& point);

        // Undocumented functions/operators
  	RNScalar *operator[](int i);

    private:
	RNScalar m[4][4];
};



/* Public variables */

extern const R4Matrix R4null_matrix;
extern const R4Matrix R4identity_matrix;



/* Utility functions */

RNScalar R4MatrixDet2 (
    RNScalar a, RNScalar b,
    RNScalar c, RNScalar d);

RNScalar R4MatrixDet3 (
    RNScalar a, RNScalar b, RNScalar c, 
    RNScalar d, RNScalar e, RNScalar f, 
    RNScalar g, RNScalar h, RNScalar i);

RNScalar R4MatrixDet4 (
    RNScalar a, RNScalar b, RNScalar c, RNScalar d, 
    RNScalar e, RNScalar f, RNScalar g, RNScalar h, 
    RNScalar i, RNScalar j, RNScalar k, RNScalar l, 
    RNScalar m, RNScalar n, RNScalar o, RNScalar p);



/* Inline functions */

inline const RNScalar *R4Matrix::
operator[] (int i) const
{
    assert ((i>=0)&&(i<=3));
    return m[i];
}



inline const RNBoolean R4Matrix::
operator==(const R4Matrix& matrix) const
{
    // Return whether or not another matrix is the same
    return (!RNCompare(this, &matrix, sizeof(R4Matrix)));
}



inline const RNBoolean R4Matrix::
operator!=(const R4Matrix& matrix) const
{
    // Return whether matrix is not equal
    return (!(*this == matrix));
}



inline void R4Matrix:: 
Transform(const R4Matrix& a)
{
    // Post-multiply transform
    *this = *this * a;
}



inline void R4Matrix:: 
Multiply(const R4Matrix& a)
{
    // Multiply matrix
    *this = *this * a;
}



inline R4Matrix& R4Matrix::
operator+=(const R4Matrix& a)
{
    // Add matrix entry-by-entry
    Add(a);
    return *this;
}



inline R4Matrix& R4Matrix::
operator-=(const R4Matrix& a)
{
    // Subtract matrix entry-by-entry
    Subtract(a);
    return *this;
}



inline R4Matrix& R4Matrix::
operator*=(const R4Matrix& a)
{
    // Multiply matrix
    Multiply(a);
    return *this;
}



inline R4Matrix 
operator*(RNScalar b, const R4Matrix& a)
{
    // Scale matrix
    return a * b;
}



inline RNScalar *R4Matrix::
operator[] (int i) 
{
    assert ((i>=0)&&(i<=3));
    return m[i];
}



