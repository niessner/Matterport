// Source file for alignment utility functions



// Include files

#include "R3Shapes/R3Shapes.h"



R3Box
R3BoundingBox(const RNArray<R3Point *>& points)
{
  // Compute bounding box
  R3Box bbox = R3null_box;
  for (int i = 0; i < points.NEntries(); i++) {
    bbox.Union(*points[i]);
  }

  // Return bounding box
  return bbox;
}



R3Point 
R3Centroid(const RNArray<R3Point *>& points, const RNScalar *weights)
{
  // Compute center of mass
  R3Point centroid(0.0, 0.0, 0.0);
  RNScalar total_weight = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    R3Point *point = points[i];
    centroid += *point * weight;
    total_weight += weight;
  }

  // Compute average
  if (total_weight > 0) centroid /= total_weight;

  // Return center of mass
  return centroid;
}



R3Triad 
R3PrincipleAxes(const R3Point& centroid, const RNArray<R3Point *>& points, const RNScalar *weights, RNScalar *variances)
{
  // Compute covariance matrix
  RNScalar m[9] = { 0 };
  RNScalar total_weight = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    RNScalar x = points[i]->X() - centroid[0];
    RNScalar y = points[i]->Y() - centroid[1];
    RNScalar z = points[i]->Z() - centroid[2];
    m[0] += weight * x*x;
    m[1] += weight * x*y;
    m[2] += weight * x*z;
    m[3] += weight * y*x;
    m[4] += weight * y*y;
    m[5] += weight * y*z;
    m[6] += weight * z*x;
    m[7] += weight * z*y;
    m[8] += weight * z*z;
    total_weight += weight;
  }

  // Normalize covariance matrix
  if (total_weight == 0) return R3xyz_triad;
  for (int i = 0; i < 9; i++) {
    m[i] /= total_weight;
  }

  // Calculate SVD of second order moments
  RNScalar U[9];
  RNScalar W[3];
  RNScalar Vt[9];
  RNSvdDecompose(3, 3, m, U, W, Vt);  // m == U . DiagonalMatrix(W) . Vt

  // Principle axes are in Vt
  R3Vector axes[3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      axes[i][j] = Vt[3*i+j];
    }
  }
  
  // Normalize all axis vectors (probably not necessary)
  RNLength length0 = axes[0].Length();
  RNLength length1 = axes[1].Length();
  RNLength length2 = axes[2].Length();
  if (RNIsPositive(length0)) axes[0] /= length0;
  if (RNIsPositive(length1)) axes[1] /= length1;
  if (RNIsPositive(length2)) axes[2] /= length2;

  // Flip axes so that "heavier" on positive side
  int positive_count[3] = { 0, 0, 0 };
  int negative_count[3] = { 0, 0, 0 };
  for (int i = 0; i < points.NEntries(); i++) {
    R3Point *point = points[i];
    for (int j = 0; j < 3; j++) {
      RNScalar dot = axes[j].Dot(point->Vector());
      if (dot > 0.0) positive_count[j]++;
      else negative_count[j]++;
    }
  }
  for (int j =0; j < 3; j++) {
    if (positive_count[j] < negative_count[j]) {
      axes[j].Flip();
    }
  }

  // Compute orthonormal triad of axes
  axes[2] = axes[0] % axes[1];

  // Just checking
  assert(RNIsEqual(axes[0].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsEqual(axes[1].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsEqual(axes[2].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsZero(axes[0].Dot(axes[1]), RN_BIG_EPSILON));
  assert(RNIsZero(axes[1].Dot(axes[2]), RN_BIG_EPSILON));
  assert(RNIsZero(axes[0].Dot(axes[2]), RN_BIG_EPSILON));

  // Return variances (eigenvalues)
  if (variances) {
    variances[0] = W[0];
    variances[1] = W[1];
    variances[2] = W[2];
  }

  // Return triad of axes
  return R3Triad(axes[0], axes[1], axes[2]);
}



RNLength
R3AverageDistance(const R3Point& center, const RNArray<R3Point *>& points, const RNScalar *weights)
{
  // Compute sum of distances between a position on the surface and a center point
  RNScalar distance = 0.0;
  RNScalar total_weight = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    R3Point *point = points[i];
    distance += weight * R3Distance(*point, center);
    total_weight += weight;
  }

  // Compute average distance
  if (total_weight > 0) distance /= total_weight;

  // Return average distance
  return distance;
}



R3Affine 
R3NormalizationTransformation(const RNArray<R3Point *>& points,  RNBoolean translate, RNBoolean rotate, int scale) 
{
  // Initialize transformation
  R3Affine affine(R3identity_affine);

  // Compute center of mass
  R3Point centroid = R3Centroid(points);

  // Translate center of mass back to original (if not translating)
  if (!translate) {
    affine.Translate(centroid.Vector());
  }

  // Scale by inverse of radius
  if ((scale != 0) && (scale != 2)) {
    RNScalar radius = R3AverageDistance(centroid, points);
    if (RNIsPositive(radius)) affine.Scale(1.0 / radius);
  }

  // Rotate to align principal axes with XYZ
  if (rotate || (scale == 2)) {
    RNScalar variances[3] = { 0 };
    R3Triad triad = R3PrincipleAxes(centroid, points, variances);
    if (!rotate) affine.Transform(R3Affine(triad.InverseMatrix()));
    if (scale == 2) {
      if (variances[0] > 0) affine.XScale(1.0 / variances[0]);
      if (variances[1] > 0) affine.YScale(1.0 / variances[1]);
      if (variances[2] > 0) affine.ZScale(1.0 / variances[2]);
    }
    affine.Transform(R3Affine(triad.InverseMatrix()));
  }

  // Translate center of mass to origin
  affine.Translate(-(centroid.Vector()));
  
  // Return normalization transformation
  return affine;
}



static R4Matrix
SetQuaternion(RNScalar p[4])
{
  R4Matrix m(R4identity_matrix);
  RNScalar l;

  if(p[0]<0){
    p[0]=-p[0];
    p[1]=-p[1];
    p[2]=-p[2];
    p[3]=-p[3];
  }
  l=p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
  if(l<.000001){return R4identity_matrix;}

  l=sqrt(l);
  p[0]/=l;
  p[1]/=l;
  p[2]/=l;
  p[3]/=l;

  m[0][0]=p[0]*p[0]+p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
  m[0][1]=2*(p[1]*p[2]+p[0]*p[3]);
  m[0][2]=2*(p[1]*p[3]-p[0]*p[2]);

  m[1][0]=2*(p[1]*p[2]-p[0]*p[3]);
  m[1][1]=p[0]*p[0]+p[2]*p[2]-p[1]*p[1]-p[3]*p[3];
  m[1][2]=2*(p[2]*p[3]+p[0]*p[1]);

  m[2][0]=2*(p[1]*p[3]+p[0]*p[2]);
  m[2][1]=2*(p[2]*p[3]-p[0]*p[1]);
  m[2][2]=p[0]*p[0]+p[3]*p[3]-p[1]*p[1]-p[2]*p[2];

  return m;
}


RNScalar 
R3AlignError(const RNArray<R3Point *>& points1, const RNArray<R3Point *>& points2, 
  const R4Matrix& rotation, const R3Point& center1, const R3Point& center2, RNScalar s1, RNScalar s2, 
  const RNScalar* weights)
{
  // Get number of points
  int count = points1.NEntries();
  if (points2.NEntries() < count) count = points2.NEntries();
  if (count == 0) return -1;

  // Compute sum of squared distances
  RNScalar ssd = 0;
  RNScalar total_weight = 0;
  for(int i = 0; i < count; i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    R3Point *point1 = points1[i];
    R3Point *point2 = points2[i];
    R3Vector v1 = (*point1 - center1) / s1;
    R3Vector v2 = (*point2 - center2) / s2;
    R3Vector v = v1 - rotation * v2;
    ssd += weight * v.Dot(v);
    total_weight += weight;
  }

  // Return RMSD
  if (total_weight == 0) return -1;
  return sqrt(ssd / total_weight);
}



RNScalar 
R3AlignError(const RNArray<R3Point *>& points1, const RNArray<R3Point *>& points2, 
  const R4Matrix& matrix, const RNScalar* weights)
{
  // Get number of points
  int count = points1.NEntries();
  if (points2.NEntries() < count) count = points2.NEntries();
  if (count == 0) return -1;

  // Compute sum of squared distances
  RNScalar ssd = 0;
  RNScalar total_weight = 0;
  for (int i = 0; i < count; i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    R3Point *point1 = points1[i];
    R3Point point2 = *(points2[i]);
    point2 = matrix * point2;
    R3Vector v = *point1 - point2;
    ssd += weight * v.Dot(v);
    total_weight += weight;
  }

  // Return RMSD
  if (total_weight == 0) return -1;
  return sqrt(ssd / total_weight);
}



R4Matrix
R3AlignPoints(const RNArray<R3Point *>& points1, const RNArray<R3Point *>& points2, 
  const RNScalar* weights, RNBoolean align_center, RNBoolean align_rotation, int align_scale)
{
  int i,j,k;

  // Get count of points
  int count = points1.NEntries();
  if (points2.NEntries() < count) count = points2.NEntries();

  // Check number of points
  if (count < 1) align_center = 0;
  if (count < 2) align_scale = 0;
  if (count < 3) align_rotation = 0;

  // Compute centers
  R3Point center1(0.0, 0.0, 0.0);
  R3Point center2(0.0, 0.0, 0.0);
  if (align_center){
    center1 = R3Centroid(points1, weights);
    center2 = R3Centroid(points2, weights);
  }

  // Compute scales
  RNScalar s1 = 1;
  RNScalar s2 = 1;
  if (align_scale){
    s1 = R3AverageDistance(center1, points1, weights);
    s2 = R3AverageDistance(center2, points2, weights);
  }

  // Compute cross-covariance of two point sets
  R4Matrix rotation = R4identity_matrix;
  if (align_rotation) {
    R4Matrix m = R4identity_matrix;
    m[0][0] = m[1][1] = m[2][2] = 0;
    RNScalar total_weight = 0;
    for (i=0; i< count; i++){
      RNScalar weight = (weights) ? weights[i] : 1;
      total_weight += weight;
      R3Point *point1 = points1[i];
      R3Point *point2 = points2[i];
      R3Vector p1 = (*point1 - center1) / s1;
      R3Vector p2 = (*point2 - center2) / s2;
      for(j=0;j<3;j++){
        for(k=0;k<3;k++){
          m[j][k] += weight * p1[j]*p2[k];
        }
      }
    }

    // Normalize cross-covariance matrix
    if (total_weight == 0) return R4identity_matrix;
    for(j=0;j<3;j++){for(k=0;k<3;k++){ m[j][k] /= total_weight; } }

    // Make cross-covariance matrix skew-symmetric
    R4Matrix a = R4identity_matrix;
    for(j=0;j<3;j++){for(k=0;k<3;k++){a[j][k]=m[j][k]-m[k][j];}}
    
    // Compute trace of cross-covariance matrix
    RNScalar trace=m[0][0]+m[1][1]+m[2][2];
    
    // Setup symmetric matrix whose eigenvectors give quaternion terms of optimal rotation
    RNScalar M[16];
    M[0]=trace;
    M[1]=M[4]=a[1][2];
    M[2]=M[8]=a[2][0];
    M[3]=M[12]=a[0][1];
    for(j=0;j<3;j++){
      for(k=0;k<3;k++){M[4*(j+1)+(k+1)]=m[j][k]+m[k][j];}
      M[4*(j+1)+(j+1)]-=trace;
    }

    // Perform SVD to get eigenvectors (quaternion terms of optimal rotation)
    RNScalar U[16];
    RNScalar W[4];
    RNScalar Vt[16];
    RNSvdDecompose(4, 4, M, U, W, Vt);  
    
    // Look at error using all eigenvectors and keep best
    int minI=0;
    R4Matrix temp[4];
    RNScalar e[4];
    for(i=0;i<4;i++){
      RNScalar p[4];
      for(j=0;j<4;j++){p[j]=U[4*j+i];}
      if(p[0]<0){for(j=0;j<4;j++){p[j]=-p[j];}}
      temp[i] = SetQuaternion(p);
      e[i]= R3AlignError(points1, points2, temp[i], center1, center2, s1, s2, weights);
      if (e[i]<e[minI]) minI=i;
    }
    rotation = temp[minI];
  }

  // Compute result
  R4Matrix result = R4identity_matrix;
  if (align_center) result.Translate(center1.Vector());
  if (align_scale) result.Scale(s1/s2);
  if (align_rotation) result.Transform(rotation);
  if (align_center) result.Translate(-(center2.Vector()));

  // Return resulting matrix that takes points2 to points1
  return result;
}



R3Plane
R3EstimatePlaneWithPCA(const RNArray<R3Point *>& points, const RNScalar *weights)
{
  // Return best fitting plane
  RNScalar variances[3];
  R3Point centroid = R3Centroid(points);
  R3Triad axes = R3PrincipleAxes(centroid, points, weights, variances);
  if (RNIsZero(variances[1])) return R3null_plane;
  return R3Plane(centroid, axes[2]);
}



R3Plane
R3EstimatePlaneWithRansac(const RNArray<R3Point *>& points, const RNScalar *weights,
  RNScalar tolerance, int max_iterations,
  RNScalar *max_inlier_fraction, RNScalar *avg_inlier_fraction)
{
  // Check number of points
  if (points.NEntries() < 3) return R3null_plane;
  
  // Try PCA plane
  RNScalar variances[3];
  RNArray<R3Point *> inliers;
  RNScalar *w = (weights) ? new RNScalar [ points.NEntries() ] : NULL;
  R3Point centroid = R3Centroid(points);
  R3Triad axes = R3PrincipleAxes(centroid, points, weights, variances);
  if (RNIsZero(variances[1])) return R3null_plane;
  R3Plane plane(centroid, axes[2]);
  if (tolerance == 0) tolerance = sqrt(variances[2]);
  for (int i = 0; i < points.NEntries(); i++) {
    R3Point *point = points[i];
    RNScalar d = R3Distance(plane, *point);
    if (d > tolerance) continue; 
    if (w) w[inliers.NEntries()] = weights[i];
    inliers.Insert(point);
  }

  // Initialize best score
  R3Plane best_plane = plane;
  RNScalar score = (RNScalar) inliers.NEntries() / (RNScalar) points.NEntries();
  RNScalar best_score = score;
  RNScalar total_score = score;
  int nscores = 1;

  // Search for a best normal
  for (int i = 0; i < max_iterations; i++) {
    // Guess normal by selecting three random points
    R3Point *p0 = points[(int) (RNRandomScalar() * points.NEntries()) ];
    R3Point *p1 = points[(int) (RNRandomScalar() * points.NEntries()) ];
    R3Point *p2 = points[(int) (RNRandomScalar() * points.NEntries()) ];
    if (R3Contains(*p0, *p1) || R3Contains(*p1, *p2) || R3Contains(*p2, *p0)) continue;
    R3Plane plane(*p0, *p1, *p2);
    
    // Find inliers
    inliers.Empty();
    for (int i = 0; i < points.NEntries(); i++) {
      R3Point *point = points[i];
      RNScalar d = R3Distance(plane, *point);
      if (d > tolerance) continue; 
      if (w) w[inliers.NEntries()] = weights[i];
      inliers.Insert(point);
    }

    // Compute  score
    RNScalar score = (RNScalar) inliers.NEntries() / (RNScalar) points.NEntries();
    total_score += score;
    nscores++;

    // Check score
    if (score <= best_score) continue;
    
    // Recompute plane with just inliers
    plane = R3EstimatePlaneWithPCA(inliers, w);
    if (plane == R3null_plane) continue;
    
    // Update best stuff
    best_plane = plane;
    best_score = score;
  }

  // Update returned results
  if (max_inlier_fraction) *max_inlier_fraction = best_score;
  if (avg_inlier_fraction) *avg_inlier_fraction = total_score / nscores;
  
  // Return best plane
  return best_plane;
}



////////////////////////////////////////////////////////////////////////



R3Box 
R3BoundingBox(int npoints, R3Point *points)
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3BoundingBox(array);
}



R3Point 
R3Centroid(int npoints, R3Point *points, const RNScalar *weights)
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3Centroid(array, weights);
}



R3Triad 
R3PrincipleAxes(const R3Point& centroid, int npoints, R3Point *points, const RNScalar *weights, RNScalar *variances)
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3PrincipleAxes(centroid, array, weights, variances);
}



RNLength
R3AverageDistance(const R3Point& center, int npoints, R3Point *points, const RNScalar *weights)
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3AverageDistance(center, array, weights);
}



R3Affine 
R3NormalizationTransformation(int npoints, R3Point *points,  RNBoolean translate, RNBoolean rotate, int scale) 
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3NormalizationTransformation(array, translate, rotate, scale);
}



RNScalar
R3AlignError(int npoints, R3Point *points1, R3Point *points2, const R4Matrix& matrix, const RNScalar *weights)
{
  RNArray<R3Point *> array1, array2;
  for (int i = 0; i < npoints; i++) { array1.Insert(&points1[i]); array2.Insert(&points2[i]); }
  return R3AlignError(array1, array2, matrix, weights);
}


R4Matrix
R3AlignPoints(int npoints, R3Point *points1, R3Point *points2, const RNScalar *weights, RNBoolean align_center, RNBoolean align_rotation, int align_scale)
{
  RNArray<R3Point *> array1, array2;
  for (int i = 0; i < npoints; i++) { array1.Insert(&points1[i]); array2.Insert(&points2[i]); }
  return R3AlignPoints(array1, array2, weights, align_center, align_rotation, align_scale);
}



R3Plane
R3EstimatePlaneWithPCA(int npoints, R3Point *points, const RNScalar *weights)
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3EstimatePlaneWithPCA(array, weights);

}



R3Plane
R3EstimatePlaneWithRansac(int npoints, R3Point *points, const RNScalar *weights,
  RNScalar tolerance, int max_iterations, RNScalar *max_inlier_fraction, RNScalar *avg_inlier_fraction)
{
  RNArray<R3Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R3EstimatePlaneWithRansac(array, weights, tolerance, max_iterations, max_inlier_fraction, avg_inlier_fraction);

}

