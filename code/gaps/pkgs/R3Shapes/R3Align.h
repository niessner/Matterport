// Include file for alignment utility functions



// Functions to align RNArray of points

R3Box R3BoundingBox(const RNArray<R3Point *>& points);
R3Point R3Centroid(const RNArray<R3Point *>& points, const RNScalar *weights = NULL);
R3Triad R3PrincipleAxes(const R3Point& centroid, const RNArray<R3Point *>& points, const RNScalar *weights = NULL, RNScalar *variances = NULL);
RNLength R3AverageDistance(const R3Point& center, const RNArray<R3Point *>& points, const RNScalar *weights = NULL);
R3Affine R3NormalizationTransformation(const RNArray<R3Point *>& points,
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);
RNScalar R3AlignError(const RNArray<R3Point *>& points1, const RNArray<R3Point *>& points2, 
  const R4Matrix& matrix = R4identity_matrix, const RNScalar* weights = NULL);
R4Matrix R3AlignPoints(const RNArray<R3Point *>& points1, const RNArray<R3Point *>& points2, const RNScalar *weights = NULL, 
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);
R3Plane R3EstimatePlaneWithPCA(const RNArray<R3Point *>& points, const RNScalar *weights = NULL);
R3Plane R3EstimatePlaneWithRansac(const RNArray<R3Point *>& points, const RNScalar *weights = NULL,
  RNScalar tolerance = 0, int niterations = 16, RNScalar *max_inlier_fraction = NULL, RNScalar *avg_inlier_fraction = NULL);



// Functions to align C array of points

R3Box R3BoundingBox(int npoints, R3Point *points);
R3Point R3Centroid(int npoints, R3Point *points, const RNScalar *weights = NULL);
R3Triad R3PrincipleAxes(const R3Point& centroid, int npoints, R3Point *points, const RNScalar *weights = NULL, RNScalar *variances = NULL);
RNLength R3AverageDistance(const R3Point& center, int npoints, R3Point *points, const RNScalar *weights = NULL);
R3Affine R3NormalizationTransformation(int npoints, R3Point *points,
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);
RNScalar R3AlignError(int npoints, R3Point *points1, R3Point *points2, 
  const R4Matrix& matrix = R4identity_matrix, const RNScalar* weights = NULL);
R4Matrix R3AlignPoints(int npoints, R3Point *points1, R3Point *points2, const RNScalar* weights = NULL, 
  RNBoolean align_center = TRUE, RNBoolean align_rotation = TRUE, int align_scale = 1);
R3Plane R3EstimatePlaneWithRansac(int npoints, R3Point *points, const RNScalar *weights = NULL,
  RNScalar tolerance = 0, int niterations = 16, RNScalar *max_inlier_fraction = NULL, RNScalar *avg_inlier_fraction = NULL);
