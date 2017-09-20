// Source file for alignment utility functions



// Include files

#include "R2Shapes/R2Shapes.h"



R2Point 
R2Centroid(const RNArray<R2Point *>& points, const RNScalar *weights)
{
  // Compute center of mass
  R2Point centroid(0.0, 0.0);
  RNScalar total_weight = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    R2Point *point = points[i];
    centroid += *point * weight;
    total_weight += weight;
  }

  // Compute average
  if (total_weight > 0) centroid /= total_weight;

  // Return center of mass
  return centroid;
}



R2Diad 
R2PrincipleAxes(const R2Point& centroid, const RNArray<R2Point *>& points, const RNScalar *weights, RNScalar *variances)
{
  // Compute covariance matrix
  RNScalar m[4] = { 0 };
  RNScalar total_weight = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    RNScalar x = points[i]->X() - centroid[0];
    RNScalar y = points[i]->Y() - centroid[1];
    m[0] += weight * x*x;
    m[1] += weight * x*y;
    m[2] += weight * y*x;
    m[3] += weight * y*y;
    total_weight += weight;
  }

  // Normalize covariance matrix
  if (total_weight == 0) return R2xy_diad;
  for (int i = 0; i < 4; i++) {
    m[i] /= total_weight;
  }

  // Calculate SVD of second order moments
  RNScalar U[4];
  RNScalar W[2];
  RNScalar Vt[4];
  RNSvdDecompose(2, 2, m, U, W, Vt);  // m == U . DiagonalMatrix(W) . Vt

  // Principle axes are in Vt
  R2Vector axes[2];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      axes[i][j] = Vt[2*i+j];
    }
  }
  
  // Normalize all axis vectors (probably not necessary)
  RNLength length0 = axes[0].Length();
  RNLength length1 = axes[1].Length();
  if (RNIsPositive(length0)) axes[0] /= length0;
  if (RNIsPositive(length1)) axes[1] /= length1;

  // Flip axes so that "heavier" on positive side
  int positive_count[2] = { 0, 0 };
  int negative_count[2] = { 0, 0 };
  for (int i = 0; i < points.NEntries(); i++) {
    R2Point *point = points[i];
    for (int j = 0; j < 2; j++) {
      RNScalar dot = axes[j].Dot(point->Vector());
      if (dot > 0.0) positive_count[j]++;
      else negative_count[j]++;
    }
  }
  for (int j =0; j < 2; j++) {
    if (positive_count[j] < negative_count[j]) {
      axes[j].Flip();
    }
  }

  // Just checking
  assert(RNIsEqual(axes[0].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsEqual(axes[1].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsZero(axes[0].Dot(axes[1]), RN_BIG_EPSILON));

  // Return variances (eigenvalues)
  if (variances) {
    variances[0] = W[0];
    variances[1] = W[1];
  }

  // Return diad of axes
  return R2Diad(axes[0], axes[1]);
}



RNLength
R2AverageDistance(const R2Point& center, const RNArray<R2Point *>& points, const RNScalar *weights)
{
  // Compute sum of distances between a position on the surface and a center point
  RNScalar distance = 0.0;
  RNScalar total_weight = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RNScalar weight = (weights) ? weights[i] : 1;
    R2Point *point = points[i];
    distance += weight * R2Distance(*point, center);
    total_weight += weight;
  }

  // Compute average distance
  if (total_weight > 0) distance /= total_weight;

  // Return average distance
  return distance;
}



R2Affine 
R2NormalizationTransformation(const RNArray<R2Point *>& points,  RNBoolean translate, RNBoolean rotate, int scale) 
{
  // Initialize transformation
  R2Affine affine(R2identity_affine);

  // Compute center of mass
  R2Point centroid = R2Centroid(points);

  // Translate center of mass back to original (if not translating)
  if (!translate) {
    affine.Translate(centroid.Vector());
  }

  // Scale by inverse of radius
  if ((scale != 0) && (scale != 2)) {
    RNScalar radius = R2AverageDistance(centroid, points);
    if (RNIsPositive(radius)) affine.Scale(1.0 / radius);
  }

  // Rotate to align principal axes with XYZ
  if (rotate || (scale == 2)) {
    RNScalar variances[3] = { 0 };
    R2Diad diad = R2PrincipleAxes(centroid, points, variances);
    if (!rotate) affine.Transform(R2Affine(diad.InverseMatrix()));
    if (scale == 2) {
      if (variances[0] > 0) affine.XScale(1.0 / variances[0]);
      if (variances[1] > 0) affine.YScale(1.0 / variances[1]);
    }
    affine.Transform(R2Affine(diad.InverseMatrix()));
  }

  // Translate center of mass to origin
  affine.Translate(-(centroid.Vector()));
  
  // Return normalization transformation
  return affine;
}



RNScalar 
R2AlignError(const RNArray<R2Point *>& points1, const RNArray<R2Point *>& points2, 
  const R3Matrix& matrix, const RNScalar* weights)
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
    R2Point *point1 = points1[i];
    R2Point point2 = *(points2[i]);
    point2 = matrix * point2;
    R2Vector v = *point1 - point2;
    ssd += weight * v.Dot(v);
    total_weight += weight;
  }

  // Return RMSD
  if (total_weight == 0) return -1;
  else return sqrt(ssd / total_weight);
}



R3Matrix
R2AlignPoints(const RNArray<R2Point *>& points1, const RNArray<R2Point *>& points2, 
  const RNScalar *weights, RNBoolean align_center, RNBoolean align_rotation, int align_scale)
{
  // Get count of points
  int count = points1.NEntries();
  if (points2.NEntries() < count) count = points2.NEntries();

  // Check number of points
  if (count < 1) align_center = 0;
  if (count < 2) align_scale = 0;
  if (count < 3) align_rotation = 0;

  // Compute centers (note: this should be weighted)
  R2Point center1(0.0, 0.0);
  R2Point center2(0.0, 0.0);
  if (align_center){
    center1 = R2Centroid(points1, weights);
    center2 = R2Centroid(points2, weights);
  }

  // Compute scales
  RNScalar s1 = 1;
  RNScalar s2 = 1;
  if (align_scale){
    s1 = R2AverageDistance(center1, points1, weights);
    s2 = R2AverageDistance(center2, points2, weights);
  }

  // Compute cross-covariance of two point sets
  R3Matrix rotation = R3identity_matrix;
  if (align_rotation) {
    // Compute covariance matrix
    RNScalar m[4] = { 0 };
    RNScalar total_weight = 0;
    for (int i = 0; i < count; i++){
      R2Point *point1 = points1[i];
      R2Point *point2 = points2[i];
      R2Vector p1 = (*point1 - center1) / s1;
      R2Vector p2 = (*point2 - center2) / s2;
      RNScalar weight = (weights) ? weights[i] : 1;
      total_weight += weight;
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++) {
          m[k*2 + j] += weight * p1[j]*p2[k];
        }
      }
    }

    // Normalize covariance matrix
    if (total_weight == 0) return R3identity_matrix;
    for (int j = 0; j < 4; j++) m[j] /= total_weight;

    // Calculate SVD of covariance matrix
    RNScalar Um[4];
    RNScalar Wm[2];
    RNScalar Vmt[4];
    RNSvdDecompose(2, 2, m, Um, Wm, Vmt);

    // https://sakai.rutgers.edu/access/content/group/7bee3f05-9013-4fc2-8743-3c5078742791/material/svd_ls_rotation.pdf
    R3Matrix Ut(Um[0], Um[2], 0, Um[1], Um[3], 0, 0, 0, 1); 
    R3Matrix V(Vmt[0], Vmt[2], 0, Vmt[1], Vmt[3], 0, 0, 0, 1); 
    R3Matrix VUt = V * Ut;
    R3Matrix D = R3identity_matrix;
    D[1][1] = R3MatrixDet2(VUt[0][0], VUt[0][1], VUt[1][0], VUt[1][1]);
    rotation = V * D * Ut;
  }

  // Compute result
  R3Matrix result = R3identity_matrix;
  if (align_center) result.Translate(center1.Vector());
  if (align_scale) result.Scale(s1/s2);
  if (align_rotation) result.Transform(rotation);
  if (align_center) result.Translate(-(center2.Vector()));

  // Return resulting matrix that takes points2 to points1
  return result;
}



////////////////////////////////////////////////////////////////////////



R2Point 
R2Centroid(int npoints, R2Point *points, const RNScalar *weights)
{
  RNArray<R2Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R2Centroid(array, weights);
}



R2Diad 
R2PrincipleAxes(const R2Point& centroid, int npoints, R2Point *points, const RNScalar *weights, RNScalar *variances)
{
  RNArray<R2Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R2PrincipleAxes(centroid, array, weights, variances);
}



RNLength
R2AverageDistance(const R2Point& center, int npoints, R2Point *points, const RNScalar *weights)
{
  RNArray<R2Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R2AverageDistance(center, array, weights);
}



R2Affine 
R2NormalizationTransformation(int npoints, R2Point *points,  RNBoolean translate, RNBoolean rotate, int scale) 
{
  RNArray<R2Point *> array;
  for (int i = 0; i < npoints; i++) array.Insert(&points[i]);
  return R2NormalizationTransformation(array, translate, rotate, scale);
}



RNScalar
R2AlignError(int npoints, R2Point *points1, R2Point *points2, const R3Matrix& matrix, const RNScalar *weights)
{
  RNArray<R2Point *> array1, array2;
  for (int i = 0; i < npoints; i++) { array1.Insert(&points1[i]); array2.Insert(&points2[i]); }
  return R2AlignError(array1, array2, matrix, weights);
}


R3Matrix
R2AlignPoints(int npoints, R2Point *points1, R2Point *points2, const RNScalar *weights, RNBoolean align_center, RNBoolean align_rotation, int align_scale)
{
  RNArray<R2Point *> array1, array2;
  for (int i = 0; i < npoints; i++) { array1.Insert(&points1[i]); array2.Insert(&points2[i]); }
  return R2AlignPoints(array1, array2, weights, align_center, align_rotation, align_scale);
}



