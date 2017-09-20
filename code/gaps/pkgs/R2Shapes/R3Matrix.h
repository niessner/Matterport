/* Include file for the GAPS matrix class */



/* Initialization functions */

int R3InitMatrix();
void R3StopMatrix();



/* Class definition */

class R3Matrix /* : public R3Base */ {
    public:
        // Constructor functions
	R3Matrix(void);
        R3Matrix(const R3Matrix& matrix);
	R3Matrix(RNScalar a00, RNScalar a01, RNScalar a02, 
   	         RNScalar a10, RNScalar a11, RNScalar a12, 
                 RNScalar a20, RNScalar a21, RNScalar a22);
	R3Matrix(const RNScalar* array);

        // Property functions/operators
	const int IsZero(void) const;
	const int IsIdentity(void) const;
        const RNBoolean IsIsotropic(void) const;
        const RNBoolean HasTranslation(void) const;
        const RNBoolean HasScale(void) const;
        const RNBoolean HasRotation(void) const;
        const RNBoolean HasMirror(void) const;
  	const RNScalar *operator[](int i) const;
        const RNScalar Determinant(void) const;
	const R3Matrix Transpose(void) const;
	const R3Matrix Inverse(void) const;
	const RNBoolean operator==(const R3Matrix& matrix) const;
	const RNBoolean operator!=(const R3Matrix& matrix) const;

        // Manipulation functions/operators
	void Flip(void);
	void Invert(void);
	void XTranslate(RNScalar offset);
	void YTranslate(RNScalar offset);
	void Translate(RNScalar offset);
	void Translate(RNAxis axis, RNScalar offset);
	void Translate(const R2Vector& offset);
	void XScale(RNScalar scale);
	void YScale(RNScalar scale);
	void Scale(RNScalar scale);
	void Scale(RNAxis axis, RNScalar scale);
	void Scale(const R2Vector& scale);
	void Rotate(RNAngle radians);
	void Transform(const R3Matrix& matrix);
	void Multiply(const R3Matrix& matrix);
	void Add(const R3Matrix& matrix);
	void Subtract(const R3Matrix& matrix);

	// Assignment operators
	R3Matrix& operator=(const R3Matrix& matrix);
	R3Matrix& operator+=(const R3Matrix& matrix);
	R3Matrix& operator-=(const R3Matrix& matrix);
	R3Matrix& operator*=(RNScalar a);
	R3Matrix& operator*=(const R3Matrix& matrix);
	R3Matrix& operator/=(RNScalar a);

        // Draw functions/operators
        void Load(void) const;
        void Draw(void) const;
        void Push(void) const;
        void Pop(void) const;

        // Arithmetic operators
	friend R3Matrix operator-(const R3Matrix& matrix);
	friend R3Matrix operator+(const R3Matrix& matrix1, const R3Matrix& matrix2);
	friend R3Matrix operator-(const R3Matrix& matrix1, const R3Matrix& matrix2);
	friend R3Matrix operator*(RNScalar a, const R3Matrix& matrix);
	friend R3Matrix operator*(const R3Matrix& matrix, RNScalar a);
	friend R3Matrix operator*(const R3Matrix& matrix1, const R3Matrix& matrix2);
	friend R3Matrix operator/(const R3Matrix& matrix, RNScalar scale);
	friend R2Vector operator*(const R3Matrix& matrix, const R2Vector& vector);
	friend R2Point operator*(const R3Matrix& matrix, const R2Point& point);

        // Undocumented functions/operators
  	RNScalar *operator[](int i);

    private:
	RNScalar m[3][3];
};



/* Public variables */

extern const R3Matrix R3null_matrix;
extern const R3Matrix R3identity_matrix;



/* Utility functions */

RNScalar R3MatrixDet2 (
    RNScalar a, RNScalar b,
    RNScalar c, RNScalar d);

RNScalar R3MatrixDet3 (
    RNScalar a, RNScalar b, RNScalar c, 
    RNScalar d, RNScalar e, RNScalar f, 
    RNScalar g, RNScalar h, RNScalar i);



/* Inline functions */

inline const RNScalar *R3Matrix::
operator[] (int i) const
{
    assert ((i>=0)&&(i<=2));
    return m[i];
}



inline const RNBoolean R3Matrix::
operator==(const R3Matrix& matrix) const
{
    // Return whether or not another matrix is the same
    return (!RNCompare(this, &matrix, sizeof(R3Matrix)));
}



inline const RNBoolean R3Matrix::
operator!=(const R3Matrix& matrix) const
{
    // Return whether matrix is not equal
    return (!(*this == matrix));
}



inline void R3Matrix:: 
Transform(const R3Matrix& a)
{
    // Post-multiply transform
    *this = *this * a;
}



inline void R3Matrix:: 
Multiply(const R3Matrix& a)
{
    // Multiply matrix
    *this = *this * a;
}



inline R3Matrix& R3Matrix::
operator+=(const R3Matrix& a)
{
    // Add matrix entry-by-entry
    Add(a);
    return *this;
}



inline R3Matrix& R3Matrix::
operator-=(const R3Matrix& a)
{
    // Subtract matrix entry-by-entry
    Subtract(a);
    return *this;
}



inline R3Matrix& R3Matrix::
operator*=(const R3Matrix& a)
{
    // Multiply matrix
    Multiply(a);
    return *this;
}



inline R3Matrix 
operator*(RNScalar b, const R3Matrix& a)
{
    // Scale matrix
    return a * b;
}



inline RNScalar *R3Matrix::
operator[] (int i) 
{
    assert ((i>=0)&&(i<=2));
    return m[i];
}



