// Include file for alignment utility functions



// Functions to align RNArray of points

R2Point R2Centroid(const RNArray<R2Point *>& points, const RNScalar *weights = NULL);
R2Diad R2PrincipleAxes(const R2Point& centroid, const RNArray<R2Point *>& points, const RNScalar *weights = NULL, RNScalar *variances = NULL);
RNLength R2AverageDistance(const R2Point& center, const RNArray<R2Point *>& points, const RNScalar *weights = NULL);
R2Affine R2NormalizationTransformation(const RNArray<R2Point *>& points,
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);
RNScalar R2AlignError(const RNArray<R2Point *>& points1, const RNArray<R2Point *>& points2, 
  const R3Matrix& matrix = R3identity_matrix, const RNScalar* weights = NULL);
R3Matrix R2AlignPoints(const RNArray<R2Point *>& points1, const RNArray<R2Point *>& poitns2, const RNScalar* weights = NULL, 
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);



// Functions to align C array of points

R2Point R2Centroid(int npoints, R2Point *points, const RNScalar *weights = NULL);
R2Diad R2PrincipleAxes(const R2Point& centroid, int npoints, R2Point *points, const RNScalar *weights = NULL, RNScalar *variances = NULL);
RNLength R2AverageDistance(const R2Point& center, int npoints, R2Point *points, const RNScalar *weights = NULL);
R2Affine R2NormalizationTransformation(int npoints, R2Point *points,
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);
RNScalar R2AlignError(int npoints, R2Point *points1, R2Point *points2,
  const R3Matrix& matrix = R3identity_matrix, RNScalar* weights = NULL);
R3Matrix R2AlignPoints(int npoints, R2Point *points1, R2Point *points2, const RNScalar* weights = NULL, 
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);



