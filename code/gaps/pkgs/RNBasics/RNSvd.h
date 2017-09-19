// Include file for SVD functions



// Function to perform singular value decomposition

extern void RNSvdDecompose(int m, int n, 
  const RNScalar *a,
  RNScalar *u, RNScalar *w, RNScalar *vt);



// Function to perform back substitution

extern void RNSvdBacksubstitute(int m, int n, 
  const RNScalar *u, const RNScalar *w, const RNScalar *vt, const RNScalar *b,
  RNScalar *x, RNScalar eps = RN_EPSILON);



// Function to solve system of equations by SVD and back substitution

extern void RNSvdSolve(int m, int n, 
  const RNScalar *a, const RNScalar *b,
  RNScalar *x, RNScalar eps = RN_EPSILON);

