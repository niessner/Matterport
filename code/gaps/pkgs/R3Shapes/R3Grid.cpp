// Source file for GAPS scalar grid class



////////////////////////////////////////////////////////////////////////
// NOTE:
// Grid values are defined as samples at the grid positions ranging from
// (0, 0, 0) to (xres-1, yres-1, zres-1).  Grid values outside this range
// are undefined.
////////////////////////////////////////////////////////////////////////



// Include files

#include "R3Shapes/R3Shapes.h"



// Useful constants

const float R3_GRID_KEEP_VALUE = R2_GRID_KEEP_VALUE;



R3Grid::
R3Grid(int xresolution, int yresolution, int zresolution)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_resolution[2] = zresolution;
  grid_row_size = xresolution;
  grid_sheet_size = grid_row_size * yresolution;
  grid_size = grid_sheet_size * zresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  grid_to_world_transform = R3identity_affine;
  world_to_grid_transform = R3identity_affine;
  world_to_grid_scale_factor = 1.0;
  grid_to_world_scale_factor = 1.0;
}



R3Grid::
R3Grid(int xresolution, int yresolution, int zresolution, const R3Box& bbox)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_resolution[2] = zresolution;
  grid_row_size = xresolution;
  grid_sheet_size = grid_row_size * yresolution;
  grid_size = grid_sheet_size * zresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  SetWorldToGridTransformation(bbox);
}



R3Grid::
R3Grid(const R3Box& bbox, RNLength spacing, int min_resolution, int max_resolution)
{
  // Check for empty bounding box
  if (bbox.IsEmpty() || (RNIsZero(spacing))) { *this = R3Grid(); return; }
  
  // Enforce max resolution
  if (max_resolution > 0) {
    if (bbox.XLength() / spacing > max_resolution) spacing = bbox.XLength() / max_resolution;
    if (bbox.YLength() / spacing > max_resolution) spacing = bbox.YLength() / max_resolution;
    if (bbox.ZLength() / spacing > max_resolution) spacing = bbox.ZLength() / max_resolution;
  }

  // Compute resolution
  grid_resolution[0] = (int) (bbox.XLength() / spacing + 0.5);
  grid_resolution[1] = (int) (bbox.YLength() / spacing + 0.5);
  grid_resolution[2] = (int) (bbox.ZLength() / spacing + 0.5);

  // Enforce min resolution
  if (min_resolution > 0) {
    if (grid_resolution[0] < min_resolution) grid_resolution[0] = min_resolution;
    if (grid_resolution[1] < min_resolution) grid_resolution[1] = min_resolution;
    if (grid_resolution[2] < min_resolution) grid_resolution[2] = min_resolution;
  }

  // Set grid resolution
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  SetWorldToGridTransformation(bbox);
}



R3Grid::
R3Grid(const R3Grid& voxels)
  : grid_values(NULL)
{
  // Copy everything
  *this = voxels;
}



R3Grid::
~R3Grid(void)
{
  // Deallocate memory for grid values
  if (grid_values) delete [] grid_values;
}



RNScalar R3Grid::
Variance(void) const
{
  // Return the variance of the values in the grid
  RNScalar sum = 0;
  RNScalar mean = Mean();
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    RNScalar delta = (*(grid_valuep++) - mean);
    sum += delta * delta;
  }
  return sum / grid_size;
}



RNScalar R3Grid::
Percentile(RNScalar percentile) const
{
  // Return value at given percentile 
  if (grid_size == 0) return 0.0;
  RNScalar *tmp_values = new RNScalar [ grid_size ];
  for (int i = 0; i < grid_size; i++) tmp_values[i] = grid_values[i];
  qsort(tmp_values, grid_size, sizeof(RNScalar), RNCompareScalars);
  int percentile_index = percentile * grid_size;
  RNScalar value = tmp_values[percentile_index];
  delete [] tmp_values;
  return value;
}



RNInterval R3Grid::
Range(void) const
{
  // Find smallest and largest values
  RNScalar minimum = FLT_MAX;
  RNScalar maximum = -FLT_MAX;
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    if (*grid_valuep < minimum) minimum = *grid_valuep;
    if (*grid_valuep > maximum) maximum = *grid_valuep;
    grid_valuep++;
  }
  return RNInterval(minimum, maximum);
}



RNScalar R3Grid::
L1Norm(void) const
{
  // Return L1 norm of grid
  RNScalar sum = 0.0;
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) 
    sum += *(grid_valuep++);
  return sum;
}



RNScalar R3Grid::
L2NormSquared(void) const
{
  // Return L2 norm of grid
  RNScalar sum = 0.0;
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = *(grid_valuep++);
    sum += value * value;
  }
  return sum;
}




int R3Grid::
Cardinality(void) const
{
  // Return number of non-zero grid values
  int count = 0;
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = *(grid_valuep++);
    if (value == 0.0) continue;
    count++;
  }
  return count;
}



RNScalar R3Grid::
GridValue(RNScalar x, RNScalar y, RNScalar z) const
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return 0.0;
  if ((y < 0) || (y > grid_resolution[1]-1)) return 0.0;
  if ((z < 0) || (z > grid_resolution[2]-1)) return 0.0;

  // Trilinear interpolation
  int ix1 = (int) x;
  int iy1 = (int) y;
  int iz1 = (int) z;
  int ix2 = ix1 + 1;
  int iy2 = iy1 + 1;
  int iz2 = iz1 + 1;
  if (ix2 >= grid_resolution[0]) ix2 = ix1;
  if (iy2 >= grid_resolution[1]) iy2 = iy1;
  if (iz2 >= grid_resolution[2]) iz2 = iz1;
  RNScalar dx = x - ix1;
  RNScalar dy = y - iy1;
  RNScalar dz = z - iz1;
  RNScalar value = 0.0;
  value += GridValue(ix1, iy1, iz1) * (1.0-dx) * (1.0-dy) * (1.0-dz);
  value += GridValue(ix1, iy1, iz2) * (1.0-dx) * (1.0-dy) * dz;
  value += GridValue(ix1, iy2, iz1) * (1.0-dx) * dy * (1.0-dz);
  value += GridValue(ix1, iy2, iz2) * (1.0-dx) * dy * dz;
  value += GridValue(ix2, iy1, iz1) * dx * (1.0-dy) * (1.0-dz);
  value += GridValue(ix2, iy1, iz2) * dx * (1.0-dy) * dz;
  value += GridValue(ix2, iy2, iz1) * dx * dy * (1.0-dz);
  value += GridValue(ix2, iy2, iz2) * dx * dy * dz;
  return value;
}



RNScalar R3Grid::
GridValue(RNScalar x, RNScalar y, RNScalar z, RNLength sigma) const
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return 0.0;
  if ((y < 0) || (y > grid_resolution[1]-1)) return 0.0;
  if ((z < 0) || (z > grid_resolution[2]-1)) return 0.0;

  // Determine Gaussian filter extent
  int f = (int) (3 * sigma + 0.5);
  int xmin = (int) (x - f + 0.5);
  int xmax = (int) (x + f + 0.5);
  int ymin = (int) (y - f + 0.5);
  int ymax = (int) (y + f + 0.5);
  int zmin = (int) (z - f + 0.5);
  int zmax = (int) (z + f + 0.5);
  if (xmin < 0) xmin = 0;
  if (xmax > XResolution()-1) xmax = XResolution()-1;
  if (ymin < 0) ymin = 0;
  if (ymax > YResolution()-1) ymax = YResolution()-1;
  if (zmin < 0) zmin = 0;
  if (zmax > ZResolution()-1) zmax = ZResolution()-1;

  // Compute Gaussian weighting variables
  const RNScalar sqrt_two_pi = sqrt(RN_TWO_PI);
  double a = sqrt_two_pi * sigma;
  double fac = 1.0 / (a * a * a);
  double denom = 2.0 * sigma * sigma;

  // Filter sample with Gaussian
  RNScalar value = 0.0;
  RNScalar weight = 0.0;
  for (int k = zmin; k <= zmax; k++) {
    for (int j = ymin; j <= ymax; j++) {
      for (int i = xmin; i <= xmax; i++) {
        RNScalar dx = x - i;
        RNScalar dy = y - j;
        RNScalar dz = z - k;
        RNScalar d = dx*dx + dy*dy + dz*dz;
        RNScalar w = fac * exp(-d / denom);
        value += w * GridValue(i, j, k);
        weight += w;
      }
    }
  }

  // Normalize value based on total weight of Gaussian samples
  if (RNIsNotZero(weight)) value /= weight;

  // Return value
  return value;
}



RNScalar R3Grid::
WorldValue(RNCoord x, RNCoord y, RNCoord z, RNLength sigma) const
{
  // Return value at world point using Gaussian filtering
  R3Point p = GridPosition(x, y, z);
  return GridValue(p[0], p[1], p[2], sigma * WorldToGridScaleFactor());
}



R3Point R3Grid::
GridCentroid(void) const
{
  // Compute weighted sum
  RNScalar total_value = 0;
  R3Point centroid(0,0,0);
  RNScalar *grid_valuesp = grid_values;
  for (int k = 0; k < grid_resolution[2]; k++) {
    for (int j = 0; j < grid_resolution[1]; j++) {
      for (int i = 0; i < grid_resolution[0]; i++) {
        R3Vector position(i, j, k);
        RNScalar value = *(grid_valuesp++);
        centroid += value * position;
        total_value += value;
      }
    }
  }

  // Divide by total value
  if (total_value > 0) centroid /= total_value;

  // Return centroid
  return centroid;
}



R3Triad R3Grid::
GridPrincipleAxes(const R3Point *grid_center, RNScalar *variances) const
{
  // Get centroid
  R3Point center = (grid_center) ? *grid_center : GridCentroid();

  // Compute covariance matrix
  RNScalar m[9] = { 0 };
  RNScalar total_value = 0;
  RNScalar *grid_valuesp = grid_values;
  for (int k = 0; k < grid_resolution[2]; k++) {
    for (int j = 0; j < grid_resolution[1]; j++) {
      for (int i = 0; i < grid_resolution[0]; i++) {
        R3Point position(i, j, k);
        RNScalar value = *(grid_valuesp++);
        RNScalar x = position[0] - center[0];
        RNScalar y = position[1] - center[1];
        RNScalar z = position[2] - center[2];
        m[0] += value * x*x;
        m[4] += value * y*y;
        m[8] += value * z*z;
        m[1] += value * x*y;
        m[3] += value * x*y;
        m[2] += value * x*z;
        m[6] += value * x*z;
        m[5] += value * y*z;
        m[7] += value * y*z;
        total_value += value;
      }
    }
  }

  // Normalize covariance matrix
  if (total_value == 0) return R3xyz_triad;
  for (int i = 0; i < 9; i++) m[i] /= total_value;

  // Compute eigenvalues and eigenvectors
  RNScalar U[9];
  RNScalar W[3];
  RNScalar Vt[9];
  RNSvdDecompose(3, 3, m, U, W, Vt);  // m == U . DiagonalMatrix(W) . Vt

  // Copy principle axes into more convenient form
  // W has eigenvalues (greatest to smallest) and Vt has eigenvectors (normalized)
  R3Vector axes[3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      axes[i][j] = Vt[3*i+j];
    }
  }

  // Flip first two axes so that "heavier" on positive side 
  grid_valuesp = grid_values;
  RNScalar positive_count[2] = { 0, 0 };
  RNScalar negative_count[2] = { 0, 0 };
  for (int k = 0; k < grid_resolution[2]; k++) {
    for (int j = 0; j < grid_resolution[1]; j++) {
      for (int i = 0; i < grid_resolution[0]; i++) {
        R3Point position(i, j, k);
        RNScalar value = *(grid_valuesp++);
        R3Vector vector = position - center;
        for (int j = 0; j < 2; j++) {
          RNScalar dot = axes[j].Dot(vector);
          if (dot > 0.0) positive_count[j] += value;
          else negative_count[j] += value;
        }
      }
    }
  }
  for (int j =0; j < 2; j++) {
    if (positive_count[j] < negative_count[j]) {
      axes[j].Flip();
    }
  }

  // Set third axis to form orthonormal triad with other two
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

  // Return triad
  return R3Triad(axes[0], axes[1], axes[2]);
}



R2Grid *R3Grid::
Slice(int dim, int grid_coordinate) const
{
  // Extract 2D grid along slice at given coordinate in given dimension

  // Compute parameters
  int dim1 = (dim+1)%3;
  int dim2 = (dim+2)%3;
  R3Box world_box = WorldBox();
  R2Box slice_box(world_box[0][dim1], world_box[0][dim2], world_box[1][dim1], world_box[1][dim2]);

  // Allocate 2D grid for slice
  R2Grid *slice = new R2Grid(grid_resolution[dim1], grid_resolution[dim2], slice_box);
  if (!slice) {
    fprintf(stderr, "Unable to allocate slice\n");
    return NULL;
  }

  // Fill values
  for (int i = 0; i < grid_resolution[dim1]; i++) {
    for (int j = 0; j < grid_resolution[dim2]; j++) {
      RNScalar value = 0.0;
      if (dim == RN_X) value = GridValue(grid_coordinate, i, j);
      else if (dim == RN_Y) value = GridValue(j, grid_coordinate, i);
      else value = GridValue(i, j, grid_coordinate);
      slice->SetGridValue(i, j, value);
    }
  }

  // Return slice
  return slice;
}



R3Grid& R3Grid::
operator=(const R3Grid& voxels) 
{
  // Delete old grid values
  if (grid_values) {
    if (grid_size != voxels.grid_size) {
      delete [] grid_values;
      grid_values = NULL;
    }
  }

  // Allocate new grid values
  if (!grid_values && (voxels.grid_size > 0)) {
    grid_values = new RNScalar [ voxels.grid_size ];
    assert(grid_values);
  }

  // Copy grid resolution
  grid_resolution[0] = voxels.grid_resolution[0];
  grid_resolution[1] = voxels.grid_resolution[1];
  grid_resolution[2] = voxels.grid_resolution[2];
  grid_row_size = voxels.grid_row_size;
  grid_sheet_size = voxels.grid_sheet_size;
  grid_size = voxels.grid_size;

  // Copy grid values
  for (int i = 0; i < grid_size; i++) {
    grid_values[i] = voxels.grid_values[i];
  }

  // Copy transforms
  grid_to_world_transform = voxels.grid_to_world_transform;
  world_to_grid_transform = voxels.world_to_grid_transform;
  world_to_grid_scale_factor = voxels.world_to_grid_scale_factor;
  grid_to_world_scale_factor = voxels.grid_to_world_scale_factor;

  // Return this
  return *this;
}



void R3Grid::
Abs(void) 
{
  // Take the absolute value of every grid value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    *grid_valuep = fabs(*grid_valuep);
    grid_valuep++;
  }
}



void R3Grid::
Sqrt(void) 
{
  // Take sqrt of every grid value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    *grid_valuep = sqrt(*grid_valuep);
    grid_valuep++;
  }
}



void R3Grid::
Square(void) 
{
  // Square every grid value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    *grid_valuep = (*grid_valuep) * (*grid_valuep);
    grid_valuep++;
  }
}



void R3Grid::
Negate(void) 
{
  // Negate every grid value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    *grid_valuep = -(*grid_valuep);
    grid_valuep++;
  }
}



void R3Grid::
Invert(void) 
{
  // Invert every grid value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    if (RNIsNotZero(*grid_valuep), 1.0E-20) *grid_valuep = 1.0/(*grid_valuep);
    grid_valuep++;
  }
}



void R3Grid::
Transpose(void)
{
  // Transpose values
  R3Grid copy(*this);
  int xres = XResolution();
  int yres = YResolution();
  int zres = ZResolution();
  for (int iz = 0; iz < zres; iz++) {
    for (int iy = 0; iy < yres; iy++) {
      for (int ix = 0; ix < xres; ix++) {
        RNScalar value = copy.GridValue(xres-1 - ix, yres-1 - iy, zres-1 - iz);
        SetGridValue(ix, iy, iz, value);
      }
    }
  }
}



void R3Grid::
Normalize(void) 
{
  // Scale so that length of "vector" is one
  Divide(L2Norm());
}



void R3Grid::
Dilate(RNLength grid_distance) 
{
  // Set voxels (to one) within grid_distance from some non-zero voxel
  SquaredDistanceTransform();
  Threshold(grid_distance * grid_distance, 1, 0);
}



void R3Grid::
Erode(RNLength grid_distance) 
{
  // Keep only voxels at least distance from some zero voxel
  Threshold(1.0E-20, 1, 0);
  SquaredDistanceTransform();
  Threshold(grid_distance * grid_distance, 0, 1);
}



void R3Grid::
Blur(RNLength grid_sigma) 
{
  // Check sigma
  if (RNIsZero(grid_sigma)) return;

  // Build filter
  RNScalar sigma = grid_sigma;
  int filter_radius = (int) (3 * sigma + 0.5);
  RNScalar *filter = new RNScalar [ filter_radius + 1 ];
  assert(filter);

  // Make buffer for temporary copy of row
  int res = XResolution();
  if (res < YResolution()) res = YResolution();
  if (res < ZResolution()) res = ZResolution();
  RNScalar *buffer = new RNScalar [ res ];
  assert(buffer);

  // Fill filter with Gaussian 
  const RNScalar sqrt_two_pi = sqrt(RN_TWO_PI);
  double a = sqrt_two_pi * sigma;
  double fac = 1.0 / (a * a * a);
  double denom = 2.0 * sigma * sigma;
  for (int i = 0; i <= filter_radius; i++) {
    filter[i] = fac * exp(-i * i / denom);
  }

  // Convolve grid with filter in X direction
  for (int k = 0; k < ZResolution(); k++) {
    for (int j = 0; j < YResolution(); j++) { 
      for (int i = 0; i < XResolution(); i++) 
        buffer[i] = GridValue(i, j, k); 
      for (int i = 0; i < XResolution(); i++) { 
        RNScalar sum = filter[0] * buffer[i];
        RNScalar weight = filter[0];
        int nsamples = i;
        if (nsamples > filter_radius) nsamples = filter_radius;
        for (int m = 1; m <= nsamples; m++) {
          sum += filter[m] * buffer[i - m];
          weight += filter[m];
        }
        nsamples = XResolution() - 1 - i;
        if (nsamples > filter_radius) nsamples = filter_radius;
        for (int m = 1; m <= nsamples; m++) {
          sum += filter[m] * buffer[i + m];
          weight += filter[m];
        }
        SetGridValue(i, j, k, sum / weight);
      }
    }
  }

  // Convolve grid with filter in Y direction
  for (int k = 0; k < ZResolution(); k++) {
    for (int j = 0; j < XResolution(); j++) { 
      for (int i = 0; i < YResolution(); i++) 
        buffer[i] = GridValue(j, i, k); 
      for (int i = 0; i < YResolution(); i++) { 
        RNScalar sum = filter[0] * buffer[i];
        RNScalar weight = filter[0];
        int nsamples = i;
        if (nsamples > filter_radius) nsamples = filter_radius;
        for (int m = 1; m <= nsamples; m++) {
          sum += filter[m] * buffer[i - m];
          weight += filter[m];
        }
        nsamples = YResolution() - 1 - i;
        if (nsamples > filter_radius) nsamples = filter_radius;
        for (int m = 1; m <= nsamples; m++) {
          sum += filter[m] * buffer[i + m];
          weight += filter[m];
        }
        SetGridValue(j, i, k, sum / weight);
      }
    }
  }

  // Convolve grid with filter in Z direction
  for (int k = 0; k < YResolution(); k++) {
    for (int j = 0; j < XResolution(); j++) { 
      for (int i = 0; i < ZResolution(); i++) 
        buffer[i] = GridValue(j, k, i); 
      for (int i = 0; i < ZResolution(); i++) { 
        RNScalar sum = filter[0] * buffer[i];
        RNScalar weight = filter[0];
        int nsamples = i;
        if (nsamples > filter_radius) nsamples = filter_radius;
        for (int m = 1; m <= nsamples; m++) {
          sum += filter[m] * buffer[i - m];
          weight += filter[m];
        }
        nsamples = ZResolution() - 1 - i;
        if (nsamples > filter_radius) nsamples = filter_radius;
        for (int m = 1; m <= nsamples; m++) {
          sum += filter[m] * buffer[i + m];
          weight += filter[m];
        }
        SetGridValue(j, k, i, sum / weight);
      }
    }
  }

  // Deallocate memory
  delete [] filter;
  delete [] buffer;
}



void R3Grid::
BilateralFilter(RNLength grid_sigma, RNLength value_sigma)
{
  // Make copy of grid
  R3Grid copy(*this);

  // Determine reasonable value sigma
  if (value_sigma == -1) {
    RNInterval range = Range();
    value_sigma = 0.01 * (range.Max() - range.Min());
  }

  // Get convenient variables
  double grid_denom = -2.0 * grid_sigma * grid_sigma;
  double value_denom = -2.0 * value_sigma * value_sigma;
  RNScalar grid_radius = 3 * grid_sigma;
  int r = (int) (grid_radius + 1);
  int r_squared = r * r;

  // Set every sample to be median of surrounding region in input grid
  for (int cz = 0; cz < ZResolution(); cz++) {
    for (int cy = 0; cy < YResolution(); cy++) {
      for (int cx = 0; cx < XResolution(); cx++) {
        // Get current value
        RNScalar value = copy.GridValue(cx, cy, cz);

        // Compute new value
        RNScalar sum = 0;
        RNScalar weight = 0;
        int zmin = cz - r;
        int zmax = cz + r;
        if (zmin < 0) zmin = 0;
        if (zmax >= ZResolution()) zmax = ZResolution() - 1;
        for (int z = zmin; z <= zmax; z++) {
          int ymin = cy - r;
          int ymax = cy + r;
          if (ymin < 0) ymin = 0;
          if (ymax >= YResolution()) ymax = YResolution() - 1;
          int dz = z - cz;
          for (int y = ymin; y <= ymax; y++) {
            int xmin = cx - r;
            int xmax = cx + r;
            if (xmin < 0) xmin = 0;
            if (xmax >= XResolution()) xmax = XResolution() - 1;
            int dy = y - cy;
            for (int x = xmin; x <= xmax; x++) {
              int dx = x - cx;
              int grid_distance_squared = dx*dx + dy*dy + dz*dz;
              if (grid_distance_squared > r_squared) continue;
              RNScalar sample = copy.GridValue(x, y, z);
              RNScalar value_distance_squared = value - sample;
              value_distance_squared *= value_distance_squared;
              RNScalar w = exp(grid_distance_squared/grid_denom) * exp(value_distance_squared/value_denom);
              sum += w * sample;
              weight += w;
            }
          }
        }

        // Set new value
        if (weight > 0) SetGridValue(cx, cy, cz, sum / weight);
      }
    }
  }
}



void R3Grid::
Convolve(const RNScalar filter[3][3][3])
{
  // Make temporary copy of grid
  R3Grid copy(*this);

  // Mark boundaries zero
  for (int j = 0; j < YResolution(); j++) { 
    for (int k = 0; k < ZResolution(); k++) { 
      SetGridValue(0,               j, k, 0);
      SetGridValue(XResolution()-1, j, k, 0);
    }
  }
  for (int i = 0; i < XResolution(); i++) { 
    for (int k = 0; k < ZResolution(); k++) { 
      SetGridValue(i, 0,               k, 0);
      SetGridValue(i, YResolution()-1, k, 0);
    }
  }
  for (int i = 0; i < XResolution(); i++) { 
    for (int j = 0; j < YResolution(); j++) { 
      SetGridValue(i, j, 0,               0);
      SetGridValue(i, j, ZResolution()-1, 0);
    }
  }

  // Convolve grid with filter
  for (int k = 1; k < ZResolution()-1; k++) {
    for (int j = 1; j < YResolution()-1; j++) { 
      for (int i = 1; i < XResolution()-1; i++) { 
        RNScalar sum = 0;
        for (int dk = -1; dk <= 1; dk++) {
          for (int dj = -1; dj <= 1; dj++) {
            for (int di = -1; di <= 1; di++) {
              RNScalar value = copy.GridValue(i + di, j + dj, k + dk);
              sum += filter[dk+1][dj+1][di+1] * value;
            }
          }
        }
        SetGridValue(i, j, k, sum);
      }
    }
  }
}



void R3Grid::
Laplacian(void)
{
  // Just a simple method for now
  const RNScalar w = 1.0 / 26.0;
  const RNScalar filter[3][3][3] = { 
    { { -w, -w, -w }, { -w,  -w, -w }, { -w, -w, -w } }, 
    { { -w, -w, -w }, { -w, 1.0, -w }, { -w, -w, -w } }, 
    { { -w, -w, -w }, { -w,  -w, -w }, { -w, -w, -w } } 
  };

  // Convolve with laplacian filter
  Convolve(filter);
}



void R3Grid::
Gradient(RNDimension dim)
{
  // Set up filters
  const RNScalar w1 = 1.0 / 12.0;
  const RNScalar w2 = 2.0 / 12.0;
  const RNScalar xfilter[3][3][3] = { 
    { { -w1, 0, w1 }, { -w1, 0, w1 }, { -w1, 0, w1 } }, 
    { { -w2, 0, w2 }, { -w2, 0, w2 }, { -w2, 0, w2 } }, 
    { { -w1, 0, w1 }, { -w1, 0, w1 }, { -w1, 0, w1 } } 
  };
  const RNScalar yfilter[3][3][3] = { 
    { { -w1, -w1, -w1 }, { 0, 0, 0 }, { w1, w1, w1 } }, 
    { { -w2, -w2, -w2 }, { 0, 0, 0 }, { w2, w2, w2 } }, 
    { { -w1, -w1, -w1 }, { 0, 0, 0 }, { w1, w1, w1 } } 
  };
  const RNScalar zfilter[3][3][3] = { 
    { { -w1, -w2, -w1 }, { -w1, -w2, -w1 }, { -w1, -w2, -w1 } }, 
    { {  0,    0,   0 }, {  0,    0,   0 }, {  0,    0,   0 } }, 
    { {  w1,  w2,  w1 }, {  w1,  w2,  w1 }, {  w1,  w2,  w1 } }, 
  };

  // Convolve with filter
  if (dim == RN_X) Convolve(xfilter);
  else if (dim == RN_Y) Convolve(yfilter);
  else if (dim == RN_Z) Convolve(zfilter);
  else RNAbort("Invalid dimension");
}



void R3Grid::
GradientMagnitude(void)
{
  // Compute magnitude of gradient
  R3Grid gx(*this); gx.Gradient(RN_X);
  R3Grid gy(*this); gy.Gradient(RN_Y);
  R3Grid gz(*this); gz.Gradient(RN_Z);
  for (int i = 0; i < grid_size; i++) {
    RNScalar x = gx.GridValue(i);
    RNScalar y = gy.GridValue(i);
    RNScalar z = gz.GridValue(i);
    grid_values[i] = sqrt(x*x + y*y + z*z);
  }
}



void R3Grid::
DetectEdges(void)
{
  Laplacian();
  Abs();
}



void R3Grid::
PercentileFilter(RNLength grid_radius, RNScalar percentile)
{
  // Make copy of grid
  R3Grid copy(*this);

  // Get convenient variables
  RNScalar grid_radius_squared = grid_radius * grid_radius;
  int r = (int) grid_radius;
  assert(r >= 0);
  int max_samples = (2*r+1) * (2*r+1) * (2*r+1);
  RNScalar *samples = new RNScalar [ max_samples ];
  assert(samples);

  // Set every sample to be Kth percentile of surrounding region in input grid
  for (int cz = 0; cz < ZResolution(); cz++) {
    for (int cy = 0; cy < YResolution(); cy++) {
      for (int cx = 0; cx < XResolution(); cx++) {
        // Build list of grid values in neighborhood
        int nsamples = 0;
        int zmin = cz - r;
        int zmax = cz + r;
        if (zmin < 0) zmin = 0;
        if (zmax >= ZResolution()) zmax = ZResolution() - 1;
        for (int z = zmin; z <= zmax; z++) {
          int dz = z - cz;
          int ymin = cy - r;
          int ymax = cy + r;
          if (ymin < 0) ymin = 0;
          if (ymax >= YResolution()) ymax = YResolution() - 1;
          for (int y = ymin; y <= ymax; y++) {
            int dy = y - cy;
            int xmin = cx - r;
            int xmax = cx + r;
            if (xmin < 0) xmin = 0;
            if (xmax >= XResolution()) xmax = XResolution() - 1;
            for (int x = xmin; x <= xmax; x++) {
              int dx = x - cx;
              int d_squared = dx*dx + dy*dy + dz*dz;
              if (d_squared > grid_radius_squared) continue;
              RNScalar sample = copy.GridValue(x, y, z);
              samples[nsamples++] = sample;
            }
          }
        }

        // Check number of grid values in neighborhood
        if (nsamples == 0) {
          SetGridValue(cx, cy, cz, 0);
        }
        else {
          // Sort samples found in neighborhood
          qsort(samples, nsamples, sizeof(RNScalar), RNCompareScalars);

          // Set grid value to percentile of neighborhood
          int index = (int) (percentile * nsamples);
          if (index < 0) index = 0;
          else if (index >= nsamples) index = nsamples-1;
          SetGridValue(cx, cy, cz, samples[index]);
        }
      }
    }
  }

  // Delete temporary memory
  delete [] samples;
}



void R3Grid::
MaskNonMinima(RNLength grid_radius)
{
  // Create grid with local minima
  R3Grid copy(*this);
  copy.MinFilter(grid_radius);

  // Mask values that are not minima
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] > copy.grid_values[i]) grid_values[i] = 0;
  }
}



void R3Grid::
MaskNonMaxima(RNLength grid_radius)
{
  // Create grid with local maxima
  R3Grid copy(*this);
  copy.MaxFilter(grid_radius);

  // Mask values that are not maxima
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] < copy.grid_values[i]) grid_values[i] = 0;
  }
}



void R3Grid::
FillHoles(int max_hole_size)
{
  // Interpolate in z
  for (int ix = 0; ix < XResolution(); ix++) {
    for (int iy = 0; iy < YResolution(); iy++) {
      int iz0 = -1;
      for (int iz1 = 0; iz1 < ZResolution(); iz1++) {
        if (GridValue(ix, iy, iz1) != 0.0) { 
          if ((iz0 >= 0) && (iz0 < iz1-1) && (iz1-iz0 < max_hole_size)) {
            RNScalar value0 = GridValue(ix, iy, iz0);
            if (value0 == 0.0) continue;
            RNScalar value1 = GridValue(ix, iy, iz1);
            if (value1 == 0.0) continue;
            for (int iz = iz0+1; iz < iz1; iz++) {
              RNScalar t = (double) (iz - iz0) / (double) (iz1 - iz0);
              RNScalar value = (1-t)*value0 + t*value1;
              SetGridValue(ix, iy, iz, value);
            }
          }
          iz0 = iz1;
        }
      }
    }
  }
  
  // Interpolate in y
  for (int iz = 0; iz < ZResolution(); iz++) {
    for (int ix = 0; ix < XResolution(); ix++) {
      int iy0 = -1;
      for (int iy1 = 0; iy1 < YResolution(); iy1++) {
        if (GridValue(ix, iy1, iz) != 0.0) { 
          if ((iy0 >= 0) && (iy0 < iy1-1) && (iy1-iy0 < max_hole_size)) {
            RNScalar value0 = GridValue(ix, iy0, iz);
            if (value0 == 0.0) continue;
            RNScalar value1 = GridValue(ix, iy1, iz);
            if (value1 == 0.0) continue;
            for (int iy = iy0+1; iy < iy1; iy++) {
              RNScalar t = (double) (iy - iy0) / (double) (iy1 - iy0);
              RNScalar value = (1-t)*value0 + t*value1;
              SetGridValue(ix, iy, iz, value);
            }
          }
          iy0 = iy1;
        }
      }
    }
  }
  
  // Interpolate in x
  for (int iy = 0; iy < YResolution(); iy++) {
    for (int iz = 0; iz < ZResolution(); iz++) {
      int ix0 = -1;
      for (int ix1 = 0; ix1 < XResolution(); ix1++) {
        if (GridValue(ix1, iy, iz) != 0.0) { 
          if ((ix0 >= 0) && (ix0 < ix1-1) && (ix1-ix0 < max_hole_size)) {
            RNScalar value0 = GridValue(ix0, iy, iz);
            if (value0 == 0.0) continue;
            RNScalar value1 = GridValue(ix1, iy, iz);
            if (value1 == 0.0) continue;
            for (int ix = ix0+1; ix < ix1; ix++) {
              RNScalar t = (double) (ix - ix0) / (double) (ix1 - ix0);
              RNScalar value = (1-t)*value0 + t*value1;
              SetGridValue(ix, iy, iz, value);
            }
          }
          ix0 = ix1;
        }
      }
    }
  }
}



void R3Grid::
Clear(RNScalar value) 
{
  // Set all grid values to value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) 
    *(grid_valuep++) = value;
}



void R3Grid::
Substitute(RNScalar old_value, RNScalar new_value) 
{
  // Replace all instances of old_value with new_value
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == old_value) {
      grid_values[i] = new_value;
    }
  }
}



void R3Grid::
Add(RNScalar value) 
{
  // Add value to all grid values 
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) 
    *(grid_valuep++) += value;
}



void R3Grid::
Copy(const R3Grid& voxels) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == voxels.grid_resolution[0]);
  assert(grid_resolution[1] == voxels.grid_resolution[1]);
  assert(grid_resolution[2] == voxels.grid_resolution[2]);

  // Copy passed grid values to corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) 
    grid_values[i] = voxels.grid_values[i];
}



void R3Grid::
Add(const R3Grid& voxels) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == voxels.grid_resolution[0]);
  assert(grid_resolution[1] == voxels.grid_resolution[1]);
  assert(grid_resolution[2] == voxels.grid_resolution[2]);

  // Add passed grid values to corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) 
    grid_values[i] += voxels.grid_values[i];
}



void R3Grid::
Add(const R3Grid& filter, const R3Point& grid_position, const R3Point& filter_position, RNScalar amplitude)
{
  // Determine extent of filter in grid coordinates
  int x1 = (int) (grid_position.X() - filter_position.X() + 1);
  int y1 = (int) (grid_position.Y() - filter_position.Y() + 1);
  int z1 = (int) (grid_position.Z() - filter_position.Z() + 1);
  int x2 = (int) (grid_position.X() + filter_position.X());
  int y2 = (int) (grid_position.Y() + filter_position.Y());
  int z2 = (int) (grid_position.Z() + filter_position.Z());
  if (x1 < 0) x1 = 0;
  if (y1 < 0) y1 = 0;
  if (z1 < 0) z1 = 0;
  if (x2 < 0) return;
  if (y2 < 0) return;
  if (z2 < 0) return;
  if (x1 >= XResolution()) return;
  if (y1 >= YResolution()) return;
  if (z1 >= ZResolution()) return;
  if (x2 >= XResolution()) x2 = XResolution()-1;
  if (y2 >= YResolution()) y2 = YResolution()-1;
  if (z2 >= ZResolution()) z2 = ZResolution()-1;

  // Add amplitude*filter to grid (aligning grid_position and filter_position)
  RNScalar sz = z1 - grid_position.Z() + filter_position.Z();
  assert((sz >= 0) && (sz < filter.ZResolution()));
  for (int gz = z1; gz <= z2; gz++, sz += 1) {
    RNScalar sy = y1 - grid_position.Y() + filter_position.Y();
    assert((sy >= 0) && (sy < filter.YResolution()));
    for (int gy = y1; gy <= y2; gy++, sy += 1) {
      RNScalar sx = x1 - grid_position.X() + filter_position.X();
      assert((sx >= 0) && (sx < filter.XResolution()));
      RNScalar *grid_valuesp = &grid_values[gz * grid_sheet_size + gy * grid_row_size + x1];
      for (int gx = x1; gx <= x2; gx++, sx += 1) {
        (*grid_valuesp++) += amplitude * filter.GridValue(sx, sy, sz);
      }
    }
  }
}



void R3Grid::
Subtract(RNScalar value) 
{
  // Add the opposite
  Add(-value);
}



void R3Grid::
Subtract(const R3Grid& voxels) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == voxels.grid_resolution[0]);
  assert(grid_resolution[1] == voxels.grid_resolution[1]);
  assert(grid_resolution[2] == voxels.grid_resolution[2]);

  // Subtract passed grid values from corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) 
    grid_values[i] -= voxels.grid_values[i];
}



void R3Grid::
Multiply(RNScalar value) 
{
  // Multiply grid values by value
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) 
    *(grid_valuep++) *= value;
}



void R3Grid::
Multiply(const R3Grid& voxels) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == voxels.grid_resolution[0]);
  assert(grid_resolution[1] == voxels.grid_resolution[1]);
  assert(grid_resolution[2] == voxels.grid_resolution[2]);

  // Multiply passed grid values by corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) 
    grid_values[i] *= voxels.grid_values[i];
}



void R3Grid::
Divide(RNScalar value) 
{
  // Just checking
  if (RNIsZero(value, 1.0E-20)) return;

  // Multiply by recipricol
  Multiply(1.0 / value);
}



void R3Grid::
Divide(const R3Grid& voxels) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == voxels.grid_resolution[0]);
  assert(grid_resolution[1] == voxels.grid_resolution[1]);
  assert(grid_resolution[2] == voxels.grid_resolution[2]);

  // Divide passed grid values by corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = voxels.grid_values[i];
    if (RNIsNotZero(value, 1.0E-20)) grid_values[i] /= value;
  }
}



void R3Grid::
Pow(RNScalar exponent) 
{
  // Raise each grid value to exponent
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = *grid_valuep;
    if (value < 0) value = -value;
    *(grid_valuep++) = pow(value, exponent);
  }
}



void R3Grid::
Mask(const R3Grid& mask) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == mask.grid_resolution[0]);
  assert(grid_resolution[1] == mask.grid_resolution[1]);
  assert(grid_resolution[2] == mask.grid_resolution[2]);

  // Set values to zero where mask is zero
  for (int i = 0; i < grid_size; i++) 
    if (mask.grid_values[i] == 0) grid_values[i] = 0;
}



void R3Grid::
Threshold(RNScalar threshold, RNScalar low, RNScalar high) 
{
  // Set grid value to low (high) if less/equal (greater) than threshold
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    if (*grid_valuep <= threshold) {
      if (low != R3_GRID_KEEP_VALUE) *grid_valuep = low;
    }
    else {
      if (high != R3_GRID_KEEP_VALUE) *grid_valuep = high;
    }
    grid_valuep++;
  }
}



void R3Grid::
Voronoi(R3Grid *squared_distance_grid)
{
  int dist;
  int* old_dist;
  int* new_dist;
  RNScalar value;
  RNScalar *old_value;
  RNScalar *new_value;
  R3Grid *dgrid;
  int res,square,tmp_dist,first;
  int x,y,z,s,t,i;

  // Allocate distance grid
  if (squared_distance_grid) dgrid = squared_distance_grid;
  else dgrid = new R3Grid(XResolution(), YResolution(), ZResolution());
  assert(dgrid);
  dgrid->SetWorldToGridTransformation(WorldToGridTransformation());

  // Allocate temporary buffers
  res = XResolution();
  if (res < YResolution()) res = YResolution();
  if (res < ZResolution()) res = ZResolution();
  old_dist = new int[res];
  new_dist = new int[res];
  old_value = new RNScalar[res];
  new_value = new RNScalar[res];
  assert(old_dist && new_dist && old_value && new_value);

  // Initalize distance grid values (0 if was set, max_value if not)
  RNScalar max_value = 3 * (res+1) * (res+1) * (res+1);
  for (i = 0; i < grid_size; i++) {
    if (grid_values[i] == 0.0) dgrid->grid_values[i] = max_value;
    else dgrid->grid_values[i] = 0.0;
  }

  // Scan along z axis
  for (x = 0; x < XResolution(); x++) {
    for (y = 0; y < YResolution(); y++) {
      first = 1;
      value = 0;
      dist = 0;
      for (z = 0; z < ZResolution(); z++) {
        if (dgrid->GridValue(x,y,z) == 0.0) {
          first=0;
          dist=0;
          value = GridValue(x,y,z);
          assert(value != 0);
        }
        else if (first == 0) {
          dist++;
          square = dist*dist;
          dgrid->SetGridValue(x, y, z, square);
          SetGridValue(x, y, z, value);
          assert(value != 0);
        }
      }
			
      // backward scan
      dist = 0;
      first = 1;
      for (z = ZResolution()-1; z >= 0; z--) {
        if (dgrid->GridValue(x,y,z) == 0.0){
          dist = 0;
          first = 0;
          value = GridValue(x,y,z);
          assert(value != 0);
        }
        else if (first == 0) {
          dist++;
          square = dist*dist;
          if (square < dgrid->GridValue(x, y, z)) {
            dgrid->SetGridValue(x, y, z, square);
            SetGridValue(x, y, z, value);
            assert(value != 0);
          }
        }
      }
    }
  } 

  // Scan along x axis
  for (z = 0; z < ZResolution(); z++) {
    for (y = 0; y < YResolution(); y++) {
      // Copy grid values
      for(x = 0; x < XResolution(); x++) {
        old_dist[x] = (int) (dgrid->GridValue(x, y, z) + 0.5);
        old_value[x] = GridValue(x, y, z);
      }
		
      // forward scan
      s = 0;
      for (x = 0; x < XResolution(); x++) {
        dist = old_dist[x];
        value = old_value[x];
        if (dist) {
          for(t = s; t <= x; t++) {
            tmp_dist = old_dist[t] + (x - t) * (x - t);
            if (tmp_dist <= dist) {
              dist = tmp_dist;
              value = old_value[t];
              s = t;
            }
          }
        }
        else {
          s = x;
        }
        new_dist[x] = dist;
        new_value[x] = value;
      }
			
      // backwards scan
      s = XResolution() - 1;
      for (x = XResolution()-1; x >= 0 ; x--) {
        dist = new_dist[x];
        value = new_value[x];
        if (dist) {
          for (t = s; t >= x; t--) {
            tmp_dist = old_dist[t] + (x - t) * (x - t);
            if (tmp_dist <= dist) {
              dist = tmp_dist;
              value = old_value[t];
              s = t;
            }
          }
          dgrid->SetGridValue(x, y, z, dist);
          SetGridValue(x, y, z, value);
        }
        else {
          s=x;
        }
      }
    }
  }
		
  // along y axis
  for (z = 0; z < ZResolution(); z++) {
    for (x = 0; x < XResolution(); x++) {
      // Copy grid values
      for (y = 0; y < YResolution(); y++) {
        old_dist[y] = (int) (dgrid->GridValue(x, y, z));
        old_value[y] = GridValue(x, y, z);
      }
			
      // forward scan
      s = 0;
      for (y = 0; y < YResolution(); y++) {
        dist = old_dist[y];
        value = old_value[y];
        if (dist) {
          for (t = s; t <= y ; t++) {
            tmp_dist = old_dist[t] + (y - t) * (y - t);
            if (tmp_dist <= dist){
              dist = tmp_dist;
              value = old_value[t];
              s = t;
            }
          }
        }
        else { 
          s = y;
        }
        new_dist[y] = dist;
        new_value[y] = value;
      }

      // backward scan
      s = YResolution() - 1;
      for(y = YResolution()-1; y >=0 ; y--) {
        dist = new_dist[y];
        value = new_value[y];
        if (dist) {
          for (t = s; t > y ; t--) {
            tmp_dist = old_dist[t] + (y - t) * (y - t);
            if (tmp_dist <= dist){
              dist = tmp_dist;
              value = old_value[t];
              s = t;
            }
          }
          dgrid->SetGridValue(x, y, z, dist);
          SetGridValue(x, y, z, value);
        }
        else { 
          s = y; 
        }
      }
    }
  }
	
  // Delete temporary buffers
  if (!squared_distance_grid) delete dgrid;
  delete[] old_dist;
  delete[] new_dist;
  delete[] old_value;
  delete[] new_value;
}



void R3Grid::
SignedDistanceTransform(void)
{
  // Compute distance from boundary into interior (negative) and into exterior (positive)
  R3Grid copy(*this);
  SquaredDistanceTransform();
  Sqrt();
  copy.Threshold(0, 1, 0);
  copy.Substitute(R2_GRID_UNKNOWN_VALUE, 1);
  copy.SquaredDistanceTransform();
  copy.Sqrt();
  Subtract(copy);
}



void R3Grid::
SquaredDistanceTransform(void)
{
  int x,y,z,s,t;
  long long dist,square,new_dist;
  long long* oldBuffer;
  long long* newBuffer;
  int first;
  int i;

  // Allocate temporary buffers
  int res = XResolution();
  if (res < YResolution()) res = YResolution();
  if (res < ZResolution()) res = ZResolution();
  oldBuffer = new long long [res];
  assert(oldBuffer);
  newBuffer = new long long[res];
  assert(newBuffer);

  // Initalize values (0 if was set, max_value if not)
  RNScalar max_value = 3.0 * (res+1) * (res+1);
  RNScalar *grid_valuesp = grid_values;
  for (i = 0; i < grid_size; i++) {
    if (*grid_valuesp == 0.0) *grid_valuesp = max_value;
    else *grid_valuesp = 0.0;
    grid_valuesp++;
  }

  // Scan along z axis
  for (x = 0; x < XResolution(); x++) {
    for (y = 0; y < YResolution(); y++) {
      first = 1;
      dist = 0;
      for (z = 0; z < ZResolution(); z++) {
        if (GridValue(x,y,z) == 0.0) {
          dist=0;
          first=0;
          SetGridValue(x, y, z, 0);
        }
        else if (first == 0) {
          dist++;
          square = dist*dist;
          SetGridValue(x, y, z, square);
        }
      }
			
      // backward scan
      dist = 0;
      first = 1;
      for (z = ZResolution()-1; z >= 0; z--) {
        if (GridValue(x,y,z) == 0.0){
          dist = 0;
          first = 0;
          SetGridValue(x, y, z, 0);
        }
        else if (first == 0) {
          dist++;
          square = dist*dist;
          if (square < GridValue(x, y, z)) {
            SetGridValue(x, y, z, square);
          }
        }
      }
    }
  }

  // Scan along x axis
  for (z = 0; z < ZResolution(); z++) {
    for (y = 0; y < YResolution(); y++) {
      // Copy grid values
      for(x = 0; x < XResolution(); x++) 
        oldBuffer[x] = (int) (GridValue(x, y, z) + 0.5);
		
      // forward scan
      s = 0;
      for (x = 0; x < XResolution(); x++) {
        dist = oldBuffer[x];
        if (dist) {
          for(t = s; t <= x; t++) {
            new_dist = oldBuffer[t] + (x - t) * (x - t);
            if (new_dist <= dist) {
              dist = new_dist;
              s = t;
            }
          }
        }
        else {
          s = x;
        }
        newBuffer[x] = dist;
      }
			
      // backwards scan
      s = XResolution() - 1;
      for (x = XResolution()-1; x >= 0 ; x--) {
        dist = newBuffer[x];
        if (dist) {
          for (t = s; t >= x; t--) {
            new_dist = oldBuffer[t] + (x - t) * (x - t);
            if (new_dist <= dist) {
              dist = new_dist;
              s = t;
            }
          }
          SetGridValue(x, y, z, dist);
        }
        else {
          s=x;
        }
      }
    }
  }
		
  // along y axis
  for (z = 0; z < ZResolution(); z++) {
    for (x = 0; x < XResolution(); x++) {
      // Copy grid values
      for (y = 0; y < YResolution(); y++)
        oldBuffer[y] = (int) (GridValue(x, y, z) + 0.5);
			
      // forward scan
      s = 0;
      for (y = 0; y < YResolution(); y++) {
        dist = oldBuffer[y];
        if (dist) {
          for (t = s; t <= y ; t++) {
            new_dist = oldBuffer[t] + (y - t) * (y - t);
            if (new_dist <= dist){
              dist = new_dist;
              s = t;
            }
          }
        }
        else { 
          s = y;
        }
        newBuffer[y] = dist;
      }

      // backward scan
      s = YResolution() - 1;
      for(y = YResolution()-1; y >=0 ; y--) {
        dist = newBuffer[y];
        if (dist) {
          for (t = s; t > y ; t--) {
            new_dist = oldBuffer[t] + (y - t) * (y - t);
            if (new_dist <= dist){
              dist = new_dist;
              s = t;
            }
          }
          SetGridValue(x, y, z, dist);
        }
        else { 
          s = y; 
        }
      }
    }
  }
	
  // Delete temporary buffers
  delete[] oldBuffer;
  delete[] newBuffer;
}



void R3Grid::
Gauss(RNLength sigma, RNBoolean square)
{
  // Check sigma
  if (sigma == 0.0) return;    

  // Replace each grid value with Gaussian of what was there before
  const RNScalar sqrt_two_pi = sqrt(RN_TWO_PI);
  RNScalar a = sqrt_two_pi * sigma;
  RNScalar fac = 1.0 / (a * a * a);
  RNScalar denom = -2.0 * sigma * sigma;
  if (RNIsZero(denom, 1.0E-6)) return;
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (square) value *= value;
    grid_values[i] = fac * exp( value / denom );
  }
}



void R3Grid::
Resample(int xresolution, int yresolution, int zresolution)
{
  // Resample grid values at new resolution
  RNScalar *new_grid_values = NULL;
  int new_grid_size = xresolution * yresolution * zresolution;
  if (new_grid_size > 0) {
    new_grid_values = new RNScalar [ new_grid_size ];
    assert(new_grid_values);
    RNScalar *new_grid_valuesp = new_grid_values;
    RNScalar xscale = (RNScalar) (grid_resolution[0]-1) / (RNScalar) (xresolution - 1);
    RNScalar yscale = (RNScalar) (grid_resolution[1]-1) / (RNScalar) (yresolution - 1);
    RNScalar zscale = (RNScalar) (grid_resolution[2]-1) / (RNScalar) (zresolution - 1);
    for (int k = 0; k < zresolution; k++) {
      RNScalar z = (k == zresolution-1) ? grid_resolution[2]-1 : k * zscale;
      for (int j = 0; j < yresolution; j++) {
        RNScalar y = (j == yresolution-1) ? grid_resolution[1]-1 : j * yscale;
        for (int i = 0; i < xresolution; i++) {
          RNScalar x = (i == xresolution-1) ? grid_resolution[0]-1 : i * xscale;
          *(new_grid_valuesp++) = GridValue(x, y, z);
        }
      }
    }
  }

  // Reset grid variables
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_resolution[2] = zresolution;
  grid_row_size = xresolution;
  grid_sheet_size = grid_row_size * yresolution;
  grid_size = grid_sheet_size * zresolution;
  if (grid_values) delete [] grid_values;
  grid_values = new_grid_values;
}



void R3Grid::
PadWithZero(int xresolution, int yresolution, int zresolution,
            int xoffset, int yoffset, int zoffset)
{
  // Add zeros to achieve desired resolution
  if ((XResolution() >= xresolution) && (YResolution() >= yresolution) && (ZResolution() >= zresolution)) return;

  // Copy this grid
  R3Grid copy(*this);

  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_resolution[2] = zresolution;
  grid_row_size = xresolution;
  grid_sheet_size = grid_row_size * yresolution;
  grid_size = grid_sheet_size * zresolution;

  // Allocate grid values
  if (grid_values) delete [] grid_values;
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Copy original values
  for (int iz = 0; iz < copy.ZResolution(); iz++) {
    if (iz + zoffset >= zresolution) continue; 
    for (int iy = 0; iy < copy.YResolution(); iy++) {
      if (iy + yoffset >= yresolution) continue; 
      for (int ix = 0; ix < copy.XResolution(); ix++) {
        if (ix + xoffset >= xresolution) continue; 
        RNScalar value = copy.GridValue(ix, iy, iz);
       SetGridValue(ix + xoffset, iy + yoffset, iz + zoffset, value);
      }
    }
  }
}



void R3Grid::
ClusterWithMeanShift(void)
{
  // Determine maximum number of iterations
  int max_iterations = XResolution();
  if (max_iterations < YResolution()) max_iterations = YResolution();
  if (max_iterations < ZResolution()) max_iterations = ZResolution();

  // Create gaussian kernel
  static RNScalar kernel[5][5][5];
  const double sigma = 1;
  const double denom = -2 * sigma * sigma;
  for (int jz = -2; jz <=2; jz++) {
    for (int jy = -2; jy <= 2; jy++) {
      for (int jx = -2; jx <= 2; jx++) {
        double dd = jx*jx + jy*jy + jz*jz;
        kernel[jx+2][jy+2][jz+2] = exp(dd/denom);
      }
    }
  }

  // Create flow
  int index = 0;
  int *flow = new int [ NEntries() ];
  for (int iz = 0; iz < ZResolution(); iz++) {
    for (int iy = 0; iy < YResolution(); iy++) {
      for (int ix = 0; ix < XResolution(); ix++) {
        flow[index] = -1;
        if (GridValue(ix, iy, iz) > 0) {
          R3Point position(ix, iy, iz);
          for (int k = 0; k < max_iterations; k++) {
            // Compute movement vector
            RNScalar weight = 0;
            R3Vector movement(0, 0, 0);
            for (int jz = -2; jz <= 2; jz++) {
              double z = position[2] + jz;
              for (int jy = -2; jy <= 2; jy++) {
                double y = position[1] + jy;
                for (int jx = -2; jx <= 2; jx++) {
                  double x = position[0] + jx;
                  RNScalar neighbor_value = GridValue(x, y, z);
                  if (neighbor_value > 0) {
                    RNScalar w = kernel[jx+2][jy+2][jz+2] * neighbor_value;
                    movement += R3Vector(w*jx, w*jy, w*jz);
                    weight += w;
                  }
                }
              }
            }
            
            // Move to centroid of neighbors
            if (weight == 0) break;
            if (RNIsZero(movement.Dot(movement))) break;
            position += movement / weight;
          }
          
          // Assign flow
          int fx = (int) (position[0] + 0.5);
          int fy = (int) (position[1] + 0.5);
          int fz = (int) (position[2] + 0.5);
          if (fx < 0) fx = 0;
          if (fy < 0) fy = 0;
          if (fz < 0) fz = 0;
          if (fx >= XResolution()) fx = XResolution() - 1;
          if (fy >= YResolution()) fx = YResolution() - 1;
          if (fz >= ZResolution()) fx = ZResolution() - 1;
          IndicesToIndex(fx, fy, fz, flow[index]);
        }
        index++;
      }
    }
  }

  // Compute cluster IDs
  int nclusters = 0;
  for (int i = 0; i < NEntries(); i++) {
    if (flow[i] == i) {
      SetGridValue(i, nclusters + 1);
      nclusters++;
    }
  }

  // Assign grid cells to cluster IDs
  for (int i = 0; i < NEntries(); i++) {
    if (flow[i] >= 0) {
      int index = flow[i];
      for (int j = 0; j < NEntries(); j++) {
        if (index == flow[index]) break;
        index = flow[index];
      }
      
      // Set grid value
      RNScalar cluster = GridValue(index);
      SetGridValue(i, cluster);
    }
  }

  // Delete flow
  delete [] flow;
}



////////////////////////////////////////////////////////////////////////
// Marching cubes code
// From http://astronomy.swin.edu.au/~pbourke/modelling/polygonise/
////////////////////////////////////////////////////////////////////////

static int 
MarchingCubes(R3Point corner_points[8], RNScalar corner_levels[8], RNScalar isolevel, R3Point *points)
{
  // Given a grid cell (defined by its corner points and corner_levels)
  // and an isolevel, calculate the triangular facets required to 
  // represent the isosurface through the cell.
  // Return the number of points added to the array "points", which
  // will be loaded up with the vertices for at most 5 triangular facets.
  // 0 will be returned if the grid cell is either totally above
  // or totally below the isolevel.

  // Initialize marching cubes edge table
  static int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   
  };

  // Initialize marching cubes triangle table
  static int triTable[256][16] =
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
     {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
     {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
     {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
     {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
     {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
     {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
     {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
     {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
     {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
     {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
     {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
     {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
     {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
     {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
     {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
     {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
     {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
     {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
     {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
     {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
     {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
     {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
     {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
     {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
     {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
     {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
     {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
     {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
     {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
     {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
     {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
     {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
     {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
     {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
     {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
     {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
     {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
     {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
     {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
     {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
     {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
     {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
     {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
     {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
     {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
     {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
     {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
     {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
     {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
     {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
     {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
     {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
     {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
     {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
     {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
     {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
     {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
     {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
     {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
     {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
     {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
     {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
     {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
     {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
     {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
     {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
     {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
     {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
     {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
     {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
     {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
     {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
     {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
     {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
     {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
     {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
     {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
     {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
     {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
     {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
     {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
     {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
     {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
     {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
     {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
     {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
     {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
     {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
     {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
     {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
     {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
     {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
     {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
     {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
     {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
     {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
     {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
     {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
     {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
     {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
     {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
     {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
     {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
     {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
     {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
     {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
     {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
     {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
     {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
     {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
     {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
     {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
     {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
     {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
     {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
     {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
     {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
     {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
     {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
     {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
     {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
     {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
     {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
     {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
     {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
     {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
     {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
     {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
     {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
     {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
     {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
     {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
     {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
     {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
     {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
     {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
     {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
     {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
     {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
     {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
     {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
     {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
     {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
     {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
     {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
     {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
     {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
     {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
     {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
     {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
     {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
     {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
     {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
     {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
     {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  /*
    Determine the index into the edge table which
    tells us which vertices are inside of the surface
  */
  int cubeindex = 0;
  if (corner_levels[0] < isolevel) cubeindex |= 1;
  if (corner_levels[1] < isolevel) cubeindex |= 2;
  if (corner_levels[2] < isolevel) cubeindex |= 4;
  if (corner_levels[3] < isolevel) cubeindex |= 8;
  if (corner_levels[4] < isolevel) cubeindex |= 16;
  if (corner_levels[5] < isolevel) cubeindex |= 32;
  if (corner_levels[6] < isolevel) cubeindex |= 64;
  if (corner_levels[7] < isolevel) cubeindex |= 128;

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0)
    return(0);

  /* Find the vertices where the surface intersects the cube */
  R3Point vertlist[12];
  if (edgeTable[cubeindex] & 1)
    vertlist[0] =
      R3Interpolate(corner_points[0],corner_points[1],corner_levels[0],corner_levels[1],isolevel);
  if (edgeTable[cubeindex] & 2)
    vertlist[1] =
      R3Interpolate(corner_points[1],corner_points[2],corner_levels[1],corner_levels[2],isolevel);
  if (edgeTable[cubeindex] & 4)
    vertlist[2] =
      R3Interpolate(corner_points[2],corner_points[3],corner_levels[2],corner_levels[3],isolevel);
  if (edgeTable[cubeindex] & 8)
    vertlist[3] =
      R3Interpolate(corner_points[3],corner_points[0],corner_levels[3],corner_levels[0],isolevel);
  if (edgeTable[cubeindex] & 16)
    vertlist[4] =
      R3Interpolate(corner_points[4],corner_points[5],corner_levels[4],corner_levels[5],isolevel);
  if (edgeTable[cubeindex] & 32)
    vertlist[5] =
      R3Interpolate(corner_points[5],corner_points[6],corner_levels[5],corner_levels[6],isolevel);
  if (edgeTable[cubeindex] & 64)
    vertlist[6] =
      R3Interpolate(corner_points[6],corner_points[7],corner_levels[6],corner_levels[7],isolevel);
  if (edgeTable[cubeindex] & 128)
    vertlist[7] =
      R3Interpolate(corner_points[7],corner_points[4],corner_levels[7],corner_levels[4],isolevel);
  if (edgeTable[cubeindex] & 256)
    vertlist[8] =
      R3Interpolate(corner_points[0],corner_points[4],corner_levels[0],corner_levels[4],isolevel);
  if (edgeTable[cubeindex] & 512)
    vertlist[9] =
      R3Interpolate(corner_points[1],corner_points[5],corner_levels[1],corner_levels[5],isolevel);
  if (edgeTable[cubeindex] & 1024)
    vertlist[10] =
      R3Interpolate(corner_points[2],corner_points[6],corner_levels[2],corner_levels[6],isolevel);
  if (edgeTable[cubeindex] & 2048)
    vertlist[11] =
      R3Interpolate(corner_points[3],corner_points[7],corner_levels[3],corner_levels[7],isolevel);

  /* Create the triangle */
  int npoints = 0;
  for (int i = 0; triTable[cubeindex][i] != -1; i+=3) {
    points[npoints++] = vertlist[triTable[cubeindex][i  ]];
    points[npoints++] = vertlist[triTable[cubeindex][i+1]];
    points[npoints++] = vertlist[triTable[cubeindex][i+2]];
  }

  // Just checking
  assert((npoints % 3) == 0);
  assert(npoints <= 3*5);

  // Return number of points
  return npoints;
}



void R3Grid::
RasterizeGridValue(int ix, int iy, int iz, RNScalar value, int operation)
{
  // Check if within bounds
  if ((ix < 0) || (ix > grid_resolution[0]-1)) return;
  if ((iy < 0) || (iy > grid_resolution[1]-1)) return;
  if ((iz < 0) || (iz > grid_resolution[2]-1)) return;

  // Update grid based on operation
  if (operation == R3_GRID_ADD_OPERATION) AddGridValue(ix, iy, iz, value);
  else if (operation == R3_GRID_SUBTRACT_OPERATION) AddGridValue(ix, iy, iz, -value);
  else if (operation == R3_GRID_REPLACE_OPERATION) SetGridValue(ix, iy, iz, value);
  else RNAbort("Unrecognized grid rasterization operation\n");
}



void R3Grid::
RasterizeGridPoint(RNScalar x, RNScalar y, RNScalar z, RNScalar value, int operation)
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return;
  if ((y < 0) || (y > grid_resolution[1]-1)) return;
  if ((z < 0) || (z > grid_resolution[2]-1)) return;

  // Trilinear interpolation
  int ix1 = (int) x;
  int iy1 = (int) y;
  int iz1 = (int) z;
  int ix2 = ix1 + 1;
  int iy2 = iy1 + 1;
  int iz2 = iz1 + 1;
  if (ix2 >= grid_resolution[0]) ix2 = ix1;
  if (iy2 >= grid_resolution[1]) iy2 = iy1;
  if (iz2 >= grid_resolution[2]) iz2 = iz1;
  RNScalar dx = x - ix1;
  RNScalar dy = y - iy1;
  RNScalar dz = z - iz1;
  RasterizeGridValue(ix1, iy1, iz1, value * (1.0-dx) * (1.0-dy) * (1.0-dz), operation);
  RasterizeGridValue(ix1, iy1, iz2, value * (1.0-dx) * (1.0-dy) * dz, operation);
  RasterizeGridValue(ix1, iy2, iz1, value * (1.0-dx) * dy * (1.0-dz), operation);
  RasterizeGridValue(ix1, iy2, iz2, value * (1.0-dx) * dy * dz, operation);
  RasterizeGridValue(ix2, iy1, iz1, value * dx * (1.0-dy) * (1.0-dz), operation);
  RasterizeGridValue(ix2, iy1, iz2, value * dx * (1.0-dy) * dz, operation);
  RasterizeGridValue(ix2, iy2, iz1, value * dx * dy * (1.0-dz), operation);
  RasterizeGridValue(ix2, iy2, iz2, value * dx * dy * dz, operation);
}



void R3Grid::
RasterizeGridPoint(RNScalar x, RNScalar y, RNScalar z, RNScalar value, RNLength sigma, int operation)
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return;
  if ((y < 0) || (y > grid_resolution[1]-1)) return;
  if ((z < 0) || (z > grid_resolution[2]-1)) return;

  // Determine filter extent
  int f = (int) (3 * sigma + 0.5);
  int xmin = (int) (x - f + 0.5);
  int xmax = (int) (x + f + 0.5);
  int ymin = (int) (y - f + 0.5);
  int ymax = (int) (y + f + 0.5);
  int zmin = (int) (z - f + 0.5);
  int zmax = (int) (z + f + 0.5);
  if (xmin < 0) xmin = 0;
  if (xmax > XResolution()-1) xmax = XResolution()-1;
  if (ymin < 0) ymin = 0;
  if (ymax > YResolution()-1) ymax = YResolution()-1;
  if (zmin < 0) zmin = 0;
  if (zmax > ZResolution()-1) zmax = ZResolution()-1;

  // Compute Gaussian weighting variables
  double a = sqrt(2.0 * RN_PI) * sigma;
  double fac = 1.0 / (a * a * a);
  double denom = 2.0 * sigma * sigma;

  // Compute total weight of Gaussian filter
  RNScalar weight = 0.0;
  for (int k = zmin; k <= zmax; k++) {
    for (int j = ymin; j <= ymax; j++) {
      for (int i = xmin; i <= xmax; i++) {
        RNScalar dx = x - i;
        RNScalar dy = y - j;
        RNScalar dz = z - k;
        RNScalar d = sqrt(dx*dx + dy*dy + dz*dz);
        RNScalar w = fac * exp(-d * d / denom);
        value += w * GridValue(i, j, k);
        weight += w;
      }
    }
  }

  // Splat values into grid (scaled by total weight of Gaussian filter)
  if (RNIsPositive(weight)) {
    for (int k = zmin; k <= zmax; k++) {
      for (int j = ymin; j <= ymax; j++) {
        for (int i = xmin; i <= xmax; i++) {
          RNScalar dx = x - i;
          RNScalar dy = y - j;
          RNScalar dz = z - k;
          RNScalar d = sqrt(dx*dx + dy*dy + dz*dz);
          RNScalar w = fac * exp(-d * d / denom);
          RasterizeGridValue(i, j, k, w * value / weight, operation);
        }
      }
    }
  }
}



void R3Grid::
RasterizeGridSpan(const int p1[3], const int p2[3], RNScalar value, int operation)
{
  // Get some convenient variables
  int d[3],p[3],dd[3],s[3];
  for (int i = 0; i < 3; i++) {
    d[i]= p2[i] - p1[i];
    if(d[i]<0){
      dd[i] = -d[i];
      s[i] = -1;
    }
    else{
      dd[i] = d[i];
      s[i] = 1;
    }
    p[i] = p1[i];
  }

  // Choose dimensions
  int i1=0;
  if(dd[1]>dd[i1]){i1=1;}
  if(dd[2]>dd[i1]){i1=2;}
  int i2=(i1+1)%3;
  int i3=(i1+2)%3;

  // Check span extent
  if(dd[i1]==0){
    // Span is a point - rasterize it
    if (((p[0] >= 0) && (p[0] < grid_resolution[0])) &&
        ((p[1] >= 0) && (p[1] < grid_resolution[1])) &&
        ((p[2] >= 0) && (p[2] < grid_resolution[2]))) {
      RasterizeGridValue(p[0], p[1], p[2], value, operation);
    }
  }
  else {
    // Step along span
    int off[3] = { 0, 0, 0 };
    for (int i = 0; i <= dd[i1]; i++) {
      if (((p[0] >= 0) && (p[0] < grid_resolution[0])) &&
          ((p[1] >= 0) && (p[1] < grid_resolution[1])) &&
          ((p[2] >= 0) && (p[2] < grid_resolution[2]))) {
        RasterizeGridValue(p[0], p[1], p[2], value, operation);
      }
      off[i2]+=dd[i2];
      off[i3]+=dd[i3];
      p[i1]+=s[i1];
      p[i2]+=s[i2]*off[i2]/dd[i1];
      p[i3]+=s[i3]*off[i3]/dd[i1];
      off[i2]%=dd[i1];
      off[i3]%=dd[i1];
    }
  }
}



void R3Grid::
RasterizeGridTriangle(const int p1[3], const int p2[3], const int p3[3], RNScalar value, int operation)
{
  int i,j;

  // Figure out the min, max, and delta in each dimension
  int mn[3], mx[3], delta[3];
  for (i = 0; i < 3; i++) {
    mx[i]=mn[i]=p1[i];
    if (p2[i] < mn[i]) mn[i]=p2[i];
    if (p3[i] < mn[i]) mn[i]=p3[i];
    if (p2[i] > mx[i]) mx[i]=p2[i];
    if (p3[i] > mx[i]) mx[i]=p3[i];
    delta[i] = mx[i] - mn[i];
  }

  // Determine direction of maximal delta
  int d = 0;
  if ((delta[1] > delta[0]) && (delta[1] > delta[2])) d = 1;
  else if (delta[2] > delta[0]) d = 2;

  // Sort by d-value
  const int *q1,*q2,*q3;
  if(p1[d]>=p2[d] && p1[d]>=p3[d]){
    q1=p1;
    if(p2[d]>=p3[d]){
      q2=p2;
      q3=p3;
    }
    else{
      q2=p3;
      q3=p2;
    }
  }
  else if(p2[d]>=p1[d] && p2[d]>=p3[d]){
    q1=p2;
    if(p1[d]>=p3[d]){
      q2=p1;
      q3=p3;
    }
    else{
      q2=p3;
      q3=p1;
    }
  }
  else{
    q1=p3;
    if(p1[d]>=p2[d]){
      q2=p1;
      q3=p2;
    }
    else{
      q2=p2;
      q3=p1;
    }
  }

  // Init state
  int dx,dx1,dx2,ddx;
  dx=q1[d]-q2[d];
  dx1=q1[d]-q2[d];
  dx2=q1[d]-q3[d];
  ddx=dx1*dx2;

  int r1[3],r2[3];
  int last1[3],last2[3];
  int off1[3],off2[3];
  for(i=0;i<3;i++){
    last1[i]=q1[i];
    last2[i]=q1[i];

    off1[i]=0;
    off2[i]=0;

    r1[i]=(-q1[i]+q2[i])*dx2;
    r2[i]=(-q1[i]+q3[i])*dx1;
  }

  // Draw Top triangle
  if(dx==0){
    for(i=0;i<3;i++){
      last1[i]=q1[i];
      last2[i]=q2[i];
    }
  }
  else{
    for(i=0;i<dx;i++){
      RasterizeGridSpan(last1,last2,value);
      for(j=0;j<3;j++){
        off1[j]+=r1[j];
        off2[j]+=r2[j];

        last1[j]+=off1[j]/ddx;
        if(off1[j]<0){off1[j]=-((-off1[j])%ddx);}
        else{off1[j]%=ddx;}

        last2[j]+=off2[j]/ddx;
        if(off2[j]<0){off2[j]=-((-off2[j])%ddx);}
        else{off2[j]%=ddx;}
      }
    }
  }

  // Init
  dx=q2[d]-q3[d];
  dx1=last1[d]-q3[d];
  dx2=last2[d]-q3[d];
  ddx=dx1*dx2;
  if(dx==0){
    RasterizeGridSpan(q2,q3,value);
    return;
  }

  for(i=0;i<3;i++){
    off1[i]=0;
    off2[i]=0;
    r1[i]=(-last1[i]+q3[i])*dx2;
    r2[i]=(-last2[i]+q3[i])*dx1;
  }

  // Draw Bottom parrallelogram
  for(i=0;i<=dx;i++){
    RasterizeGridSpan(last1,last2,value);
    for(j=0;j<3;j++){
      off1[j]+=r1[j];
      off2[j]+=r2[j];

      last1[j]+=off1[j]/ddx;
      if(off1[j]<0){off1[j]=-((-off1[j])%ddx);}
      else{off1[j]%=ddx;}

      last2[j]+=off2[j]/ddx;
      if(off2[j]<0){off2[j]=-((-off2[j])%ddx);}
      else{off2[j]%=ddx;}
    }
  }
}



void R3Grid::
RasterizeGridSphere(const R3Point& center, RNLength radius, RNScalar value, RNBoolean solid, int operation)
{
  // Figure out the min and max in each dimension
  int mn[3], mx[3];
  for (int i = 0; i < 3; i++) {
    mx[i]= (int) (center[i]+radius);
    if (mx[i] < 0) return;
    if (mx[i] > Resolution(i)-1) mx[i] = Resolution(i)-1;
    mn[i]= (int) (center[i]-radius);
    if (mn[i] > Resolution(i)-1) return;
    if (mn[i] < 0) mn[i] = 0;
  }

  // Rasterize sphere
  RNScalar radius_squared = radius * radius;
  for (int k = mn[2]; k <= mx[2]; k++) {
    RNCoord z = (int) (k - center[2]);
    RNCoord z_squared = z*z;
    RNLength xy_radius_squared = radius_squared - z_squared;
    RNLength y = sqrt(xy_radius_squared);
    int y1 = (int) (center[1] - y + 0.5);
    int y2 = (int) (center[1] + y + 0.5);
    if (y1 < mn[1]) y1 = mn[1];
    if (y2 > mx[1]) y2 = mx[1];
    for (int j = y1; j <= y2; j++) {
      RNCoord y = (int) (j - center[1]);
      RNCoord y_squared = y*y;
      RNLength x_squared = xy_radius_squared - y_squared;
      RNLength x = sqrt(x_squared);
      int x1 = (int) (center[0] - x + 0.5);
      int x2 = (int) (center[0] + x + 0.5);
      if (x1 < mn[0]) x1 = mn[0];
      if (x2 > mx[0]) x2 = mx[0];
      if (solid || (j == y1) || (j == y2)) {
        for (int i = x1; i <= x2; i++) {
          RasterizeGridValue(i, j, k, value, operation);
        }
      }
      else {
        RasterizeGridValue(x1, j, k, value, operation);
        RasterizeGridValue(x2, j, k, value, operation);
      }
    }
  }
}



void R3Grid::
RasterizeGridPlane(const R3Plane& plane, RNScalar value, int operation)
{
  // Load info to pass to MarchingCubes
  R3Point corner_points[8];
  RNScalar corner_levels[8];
  R3Point triangle_points[3 * 8];
  R3Box grid_box = GridBox();
  for (int octant = 0; octant < 8; octant++) {
    corner_points[octant] = grid_box.Corner(octant);
    corner_levels[octant] = R3SignedDistance(plane, corner_points[octant]);
  }

  // Determine triangles 
  int npoints = MarchingCubes(corner_points, corner_levels, 0, triangle_points);

  // Rasterize triangles into temporary grid
  R3Grid tmp(grid_resolution[0], grid_resolution[1], grid_resolution[2]);
  for (int i = 0; i < npoints; i += 3) {
    tmp.RasterizeGridTriangle(triangle_points[i], triangle_points[i+1], triangle_points[i+2], 1);
  }
  
  // Add value to this grid where temporary grid is nonzero (avoids double counting at triangle boundaries)
  for (int iz = 0; iz < grid_resolution[2]; iz++) {
    for (int iy = 0; iy < grid_resolution[1]; iy++) {
      for (int ix = 0; ix < grid_resolution[0]; ix++) {
        if (tmp.GridValue(ix, iy, iz) > 0) {
          RasterizeGridValue(ix, iy, iz, value, operation);
        }
      }
    }
  }
}



void R3Grid::
RasterizeGridBox(const R3Box& box, RNScalar value, int operation)
{
  // Get corner coordinates
  int ix1 = (int) (box.XMin() + 0.5);
  int ix2 = (int) (box.XMax() + 0.5);
  int iy1 = (int) (box.YMin() + 0.5);
  int iy2 = (int) (box.YMax() + 0.5);
  int iz1 = (int) (box.ZMin() + 0.5);
  int iz2 = (int) (box.ZMax() + 0.5);

  // Check for intersection
  if (ix1 >= XResolution()) return;
  if (iy1 >= YResolution()) return;
  if (iz1 >= ZResolution()) return;
  if (ix2 < 0) return;
  if (iy2 < 0) return;
  if (iz2 < 0) return;

  // Clip coordinates
  if (ix1 < 0) ix1 = 0;
  if (iy1 < 0) iy1 = 0;
  if (iz1 < 0) iz1 = 0;
  if (ix2 >= XResolution()) ix2 = XResolution()-1;
  if (iy2 >= YResolution()) iy2 = YResolution()-1;
  if (iz2 >= ZResolution()) iz2 = ZResolution()-1;

  // Add value to grid entries
  for (int iz = iz1; iz <= iz2; iz++) {
    for (int iy = iy1; iy <= iy2; iy++) {
      for (int ix = ix1; ix <= ix2; ix++) {
        RasterizeGridValue(ix, iy, iz, value, operation);
      }
    }
  }
}




RNScalar R3Grid::
Dot(const R3Grid& voxels) const
{
  // Resolutions and transforms must be the same (for now)
  assert(grid_resolution[0] == voxels.grid_resolution[0]);
  assert(grid_resolution[1] == voxels.grid_resolution[1]);
  assert(grid_resolution[2] == voxels.grid_resolution[2]);

  // Compute dot product between this and grid
  RNScalar dot = 0.0;
  for (int i = 0; i < grid_size; i++) 
    dot += grid_values[i] * voxels.grid_values[i];
  return dot;
}



RNScalar R3Grid::
L1Distance(const R3Grid& voxels) const
{
  // Compute distance between this and grid
  RNScalar distance = 0.0;
  for (int i = 0; i < grid_size; i++) 
    distance += fabs(grid_values[i] - voxels.grid_values[i]);
  return distance;
}



RNScalar R3Grid::
L2DistanceSquared(const R3Grid& voxels) const
{
  // Compute distance between this and grid
  RNScalar distance_squared = 0.0;
  for (int i = 0; i < grid_size; i++) {
    RNScalar delta = (grid_values[i] - voxels.grid_values[i]);
    distance_squared += delta * delta;
  }

  // Return result
  return distance_squared;
}



void R3Grid::
SetWorldToGridTransformation(const R3Affine& affine)
{
  // Set transformations
  world_to_grid_transform = affine;
  grid_to_world_transform = affine.Inverse();
  world_to_grid_scale_factor = affine.ScaleFactor();
  grid_to_world_scale_factor = (world_to_grid_scale_factor != 0) ? 1 / world_to_grid_scale_factor : 1.0;
}



void R3Grid::
SetWorldToGridTransformation(const R3Box& world_box)
{
  // Just checking
  if (grid_size == 0) return;
  if (world_box.NDimensions() < 3) return;

  // Compute grid origin
  R3Vector grid_diagonal(XResolution()-1, YResolution()-1, ZResolution()-1);
  R3Vector grid_origin = 0.5 * grid_diagonal;

  // Compute world origin
  R3Vector world_diagonal(world_box.XLength(), world_box.YLength(), world_box.ZLength());
  R3Vector world_origin = world_box.Centroid().Vector();

  // Compute scale
  RNScalar scale = FLT_MAX;
  RNScalar xscale = (world_diagonal[0] > 0) ? grid_diagonal[0] / world_diagonal[0] : FLT_MAX;
  if (xscale < scale) scale = xscale;
  RNScalar yscale = (world_diagonal[1] > 0) ? grid_diagonal[1] / world_diagonal[1] : FLT_MAX;
  if (yscale < scale) scale = yscale;
  RNScalar zscale = (world_diagonal[2] > 0) ? grid_diagonal[2] / world_diagonal[2] : FLT_MAX;
  if (zscale < scale) scale = zscale;
  if (scale == FLT_MAX) scale = 1;

  // Compute world-to-grid transformation
  R3Affine affine(R3identity_affine);
  affine.Translate(grid_origin);
  if (scale != 1) affine.Scale(scale);
  affine.Translate(-world_origin);

  // Set transformations
  SetWorldToGridTransformation(affine);
}



void R3Grid::
SetWorldToGridTransformation(const R3Point& world_origin, const R3Vector& world_axis1, const R3Vector& world_axis2, RNLength world_radius)
{
  // Just checking
  if (grid_size == 0) return;

  // Compute grid origin
  R3Vector grid_diagonal(XResolution()-1, YResolution()-1, ZResolution()-1);
  R3Vector grid_origin = 0.5 * grid_diagonal;
  RNScalar grid_radius = grid_origin[0];
  if (grid_origin[1] < grid_radius) grid_radius = grid_origin[1];
  if (grid_origin[2] < grid_radius) grid_radius = grid_origin[2];

  // Compute scale
  if (RNIsZero(world_radius)) return;
  if (RNIsZero(grid_radius)) return;
  RNScalar scale = grid_radius / world_radius;

  // Compute rotation
  R3Triad world_triad(world_axis1, world_axis2);
  R3Affine rotation(world_triad.InverseMatrix());

  // Compute world-to-grid transformation
  R3Affine affine(R3identity_affine);
  affine.Translate(grid_origin);
  affine.Transform(rotation);
  affine.Scale(scale);
  affine.Translate(-(world_origin.Vector()));

  // Set transformations
  SetWorldToGridTransformation(affine);
}



RNScalar R3Grid::
WorldSpacing(RNDimension dim) const
{
  // Return distance between grid samples in dimension
  assert((0 <= dim) && (dim <= 2));
  RNScalar m0 = grid_to_world_transform.Matrix()[0][dim];
  RNScalar m1 = grid_to_world_transform.Matrix()[1][dim];
  RNScalar m2 = grid_to_world_transform.Matrix()[2][dim];
  RNScalar spacing_squared = m0*m0 + m1*m1 + m2*m2;
  return sqrt(spacing_squared);
}



R3Point R3Grid::
WorldPosition(RNCoord x, RNCoord y, RNCoord z) const
{
  // Transform point from grid coordinates to world coordinates
  R3Point world_point(x, y, z);
  world_point.Transform(grid_to_world_transform);
  return world_point;

}



R3Point R3Grid::
GridPosition(RNCoord x, RNCoord y, RNCoord z) const
{
  // Transform point from world coordinates to grid coordinates
  R3Point grid_point(x, y, z);
  grid_point.Transform(world_to_grid_transform);
  return grid_point;
}



int R3Grid::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  int status = 0;
  if (!strncmp(extension, ".grd", 4)) 
    status = ReadGridFile(filename);
  else if (!strncmp(extension, ".phi", 4)) 
    status = ReadDelphiFile(filename);
  else if (!strncmp(extension, ".vox", 4)) 
    status = ReadVoxelFile(filename);
  else if (!strncmp(extension, ".ins", 4)) 
    status = ReadInsightFile(filename);
  else if (!strncmp(extension, ".dx", 3)) 
    status = ReadDXFile(filename);
  else if (!strncmp(extension, ".map", 4)) 
    status = ReadCCP4File(filename);
  else if (!strncmp(extension, ".raw", 4)) 
    status = ReadRawFile(filename);
  else {
    RNFail("Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return status
  return status;
}



int R3Grid::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  int status = 0;
  if (!strncmp(extension, ".grd", 4)) 
    status = WriteGridFile(filename);
  else if (!strncmp(extension, ".dx", 3)) 
    status = WriteDXFile(filename);
  else if (!strncmp(extension, ".raw", 4)) 
    status = WriteRawFile(filename);
  else if (!strncmp(extension, ".pdb", 4)) 
    status = WritePDBFile(filename);
  else {
    RNFail("Unable to write file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return status
  return status;
}



////////////////////////////////////////////////////////////////////////

int R3Grid::
ReadGridFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open voxel file %s", filename);
    return 0;
  }

  // Read 
  int status = ReadGrid(fp);

  // Close file
  fclose(fp);

  // Return status
  return status;
}



int R3Grid::
WriteGridFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    RNFail("Unable to open voxel file %s", filename);
    return 0;
  }

  // Write
  int status = WriteGrid(fp);

  // Close file
  fclose(fp);

  // Return status
  return status;
}



int R3Grid::
ReadGrid(FILE *fp)
{
  // Check file
  if (!fp) fp = stdin;

  // Read grid resolution from file
  int res[3];
  if (fread(res, sizeof(int), 3, fp) != 3) {
    RNFail("Unable to read resolution from file");
    return 0;
  }

  // Re-allocate grid values
  int new_size = res[0] * res[1] * res[2];
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }

  // Update grid resolution variables
  grid_resolution[0] = res[0];
  grid_resolution[1] = res[1];
  grid_resolution[2] = res[2];
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read world_to_grid transformation to file
  RNScalar *m = (RNScalar *) &(world_to_grid_transform.Matrix()[0][0]);
  for (int i = 0; i < 16; i++) {
    float value;
    if (fread(&value, sizeof(float), 1, fp) != 1) {
      RNFail("Unable to read transformation matrix value to file");
      return 0;
    }
    else {
      *(m++) = (RNScalar) value;
    }
  }

  // Read grid values
  RNScalar *grid_valuesp = grid_values;
  for (int i = 0; i < grid_size; i++) {
    float value;
    if (fread(&value, sizeof(float), 1, fp) != 1) {
      RNFail("Unable to read grid value %d of %d from file", i, grid_size);
      return 0;
    }
    else {
      *(grid_valuesp++) = (RNScalar) value;
    }
  }

  // Update transformation variables
  world_to_grid_scale_factor = world_to_grid_transform.ScaleFactor();
  grid_to_world_scale_factor = (world_to_grid_scale_factor != 0) ? 1 / world_to_grid_scale_factor : 1.0;
  grid_to_world_transform = world_to_grid_transform.Inverse();

  // Return number of grid values read
  return grid_size;
}



int R3Grid::
WriteGrid(FILE *fp) const
{
  // Check file
  if (!fp) fp = stdout;

  // Write grid resolution from file
  if (fwrite(&grid_resolution, sizeof(int), 3, fp) != 3) {
    RNFail("Unable to write resolution to file");
    return 0;
  }

  // Write world_to_grid transformation to file
  const RNScalar *m = &(world_to_grid_transform.Matrix()[0][0]);
  for (int i = 0; i < 16; i++) {
    float value = (float) *(m++);
    if (fwrite(&value, sizeof(float), 1, fp) != 1) {
      RNFail("Unable to write transformation matrix value to file");
      return 0;
    }
  }

  // Write grid values
  const RNScalar *grid_valuesp = grid_values;
  for (int i = 0; i < grid_size; i++) {
    float value = (float) *(grid_valuesp++);
    if (fwrite(&value, sizeof(float), 1, fp) != 1) {
      RNFail("Unable to write grid value to file");
      return 0;
    }
  }

  // Return number of grid values written
  return grid_size;
}



////////////////////////////////////////////////////////////////////////

static int 
ReadBinaryRawSizeFile(const char *filename, int& xres, int& yres, int& zres, char *format)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open size file %s", filename);
    return 0;
  }

  // Read dimensions
  unsigned int values[3];
  if (fread(values, sizeof(unsigned int), 3, fp) != 3) {
    RNFail("Unable to read size file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Assign return values
  xres = values[0];
  yres = values[1];
  zres = values[2];
  strncpy(format, "Double64", 4096);

  // Return success
  return 1;
}



static int 
ReadAsciiRawSizeFile(const char *filename, int& xres, int& yres, int& zres, char *format)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open size file %s", filename);
    return 0;
  }

  // Read line from file
  char buffer[4096];
  if (!fgets(buffer, 4096, fp)) {
    RNFail("Unable to read size file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Check line (may be binary file)
  if (!strstr(buffer, "(")) return 0;
  if (!strstr(buffer, ")")) return 0;

  // Parse format
  char *token = strtok(buffer, "(");
  if (!token) { RNFail("Error parsing type name in size file %s", filename); return 0; }
  strncpy(format, token, 4096);

  // Parse xres
  token = strtok(NULL, ",)");
  if (!token) { RNFail("Error parsing xres in size file %s", filename); return 0; }
  xres = atoi(token);

  // Parse yres
  token = strtok(NULL, ",)");
  if (!token) { RNFail("Error parsing yres in size file %s", filename); return 0; }
  yres = atoi(token);

  // Parse zres
  token = strtok(NULL, ",)");
  if (!token) { RNFail("Error parsing zres in size file %s", filename); return 0; }
  zres = atoi(token);

  // Return success
  return 1;
}



static int 
ReadRawSizeFile(const char *filename, int& xres, int& yres, int& zres, char *format)
{
  // Try both binary and ascii files
  if (ReadAsciiRawSizeFile(filename, xres, yres, zres, format)) return 1;
  if (ReadBinaryRawSizeFile(filename, xres, yres, zres, format)) return 1;
  return 0;
}



static int 
ReadRawTransformationFile(const char *filename, R3Affine& transformation)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) return 0;

  // Read transformation matrix
  RNScalar m[16];
  if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
    &m[0], &m[1], &m[2], &m[3],     &m[4], &m[5], &m[6], &m[7],
    &m[8], &m[9], &m[10], &m[11],   &m[12], &m[13], &m[14], &m[15]) != 16) {
    RNFail("Unable to read transformation file %s", filename);
    return 0;
  }
 
  // Close file
  fclose(fp);
   
  // Fill in transformation
  transformation.Reset(R4Matrix(m), 0);

  // Return success
  return 1;
}



template <class T>
static int
ReadRawValues(FILE *fp, RNScalar *values, int nvalues)
{
  // Read values of type T into array of RNScalar
  for (int i = 0; i < nvalues; i++) {
    T value;
    if (fread(&value, sizeof(T), 1, fp) != 1) return 0;
    values[i] = (RNScalar) value;
  }
  return 1;
}



int R3Grid::
ReadRawFile(const char *filename)
{
  // Get size file name
  char size_name[4096];
  strncpy(size_name, filename, 4096);
  char *ext = strrchr(size_name, '.');
  if (ext) *ext = '\0';
  strncat(size_name, ".size", 4096);

  // Get transformation file name
  char transformation_name[4096];
  strncpy(transformation_name, filename, 4096);
  ext = strrchr(transformation_name, '.');
  if (ext) *ext = '\0';
  strncat(transformation_name, ".xf", 4096);

  // Read size file
  char format[4096];
  int xres, yres, zres;
  if (!ReadRawSizeFile(size_name, xres, yres, zres, format)) return 0;
  int new_size = xres * yres * zres;
  if (new_size == 0) {
    RNFail("Zero dimensions in size file %s", size_name);
    return 0;
  }

  // Open raw file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open raw file %s", filename);
    return 0;
  }

  // Re-allocate grid values
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }

  // Update grid resolution variables
  grid_resolution[0] = xres;
  grid_resolution[1] = yres;
  grid_resolution[2] = zres;
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read raw grid values
  if (!strcmp(format, "Double64") || !strcmp(format, "double64")) {
    if (!ReadRawValues<RNScalar64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Float64") || !strcmp(format, "float64")) {
    if (!ReadRawValues<RNScalar64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Float32") || !strcmp(format, "float32")) {
    if (!ReadRawValues<RNScalar32>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uint64") || !strcmp(format, "Uint64")) {
    if (!ReadRawValues<RNUInt64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Int64") || !strcmp(format, "int64")) {
    if (!ReadRawValues<RNInt64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uint32") || !strcmp(format, "uint32")) {
    if (!ReadRawValues<RNUInt32>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Int32") || !strcmp(format, "int32")) {
    if (!ReadRawValues<RNInt32>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uint16") || !strcmp(format, "uint16")) {
    if (!ReadRawValues<RNUInt16>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Int16") || !strcmp(format, "int16")) {
    if (!ReadRawValues<RNInt16>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uchar8") || !strcmp(format, "uchar8")) {
    if (!ReadRawValues<RNUChar8>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Char8") || !strcmp(format, "char8")) {
    if (!ReadRawValues<RNChar8>(fp, grid_values, grid_size)) return 0;
  } else {
    fprintf(stderr, "Unrecognized format %s in %s\n", format, size_name);
    return 0;
  }

  // Close file
  fclose(fp);

  // Read transformation file
  R3Affine transformation;
  if (ReadRawTransformationFile(transformation_name, transformation)) {
    SetWorldToGridTransformation(transformation);
  }

  // Return success
  return 1;
}



static int 
WriteBinaryRawSizeFile(const char *filename, int xres, int yres, int zres, const char *format)
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    RNFail("Unable to open size file %s", filename);
    return 0;
  }

  // Write dimensions
  int ngrids = 1;
  fwrite(&xres, sizeof(int), 1, fp);
  fwrite(&yres, sizeof(int), 1, fp);
  fwrite(&zres, sizeof(int), 1, fp);
  fwrite(&ngrids, sizeof(int), 1, fp);

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int 
WriteAsciiRawSizeFile(const char *filename, int xres, int yres, int zres, const char *format)
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open size file %s", filename);
    return 0;
  }

  // Write dimensions
  fprintf(fp, "%s(%d,%d,%d,1)\n", format, xres, yres, zres);

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int 
WriteRawSizeFile(const char *filename, int xres, int yres, int zres, const char *format)
{
  // Check format
  if (!strcmp(format, "ZNN")) {
    if (WriteBinaryRawSizeFile(filename, xres, yres, zres, "Double64")) return 1;
  }
  else {
    if (WriteAsciiRawSizeFile(filename, xres, yres, zres, format)) return 1;
  }
  return 0;
}



static int 
WriteRawTransformationFile(const char *filename, const R3Affine& transformation)
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open transformation file %s", filename);
    return 0;
  }

  // Write transformation matrix
  const R4Matrix& m = transformation.Matrix();
  fprintf(fp, "%18.12f %18.12f %18.12f %18.12f\n", m[0][0], m[0][1], m[0][2], m[0][3]);
  fprintf(fp, "%18.12f %18.12f %18.12f %18.12f\n", m[1][0], m[1][1], m[1][2], m[1][3]);
  fprintf(fp, "%18.12f %18.12f %18.12f %18.12f\n", m[2][0], m[2][1], m[2][2], m[2][3]);
  fprintf(fp, "%18.12f %18.12f %18.12f %18.12f\n", m[3][0], m[3][1], m[3][2], m[3][3]);

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



template <class T>
static int
WriteRawValues(FILE *fp, RNScalar *values, int nvalues)
{
  // Write values of type T from array of RNScalar
  for (int i = 0; i < nvalues; i++) {
    T value = (T) values[i];
    if (fwrite(&value, sizeof(T), 1, fp) != 1) return 0;
  }
  return 1;
}



int R3Grid::
WriteRawFile(const char *filename, const char *format) const
{
  // Get size file name
  char size_name[4096];
  strncpy(size_name, filename, 4096);
  char *ext = strrchr(size_name, '.');
  if (ext) *ext = '\0';
  strncat(size_name, ".size", 4096);

  // Get transformation file name
  char transformation_name[4096];
  strncpy(transformation_name, filename, 4096);
  ext = strrchr(transformation_name, '.');
  if (ext) *ext = '\0';
  strncat(transformation_name, ".xf", 4096);

  // Write size file
  if (!WriteRawSizeFile(size_name, grid_resolution[0], grid_resolution[1], grid_resolution[2], format)) return 0;

  // Write transformation file
  if (!WriteRawTransformationFile(transformation_name, world_to_grid_transform)) return 0;

  // Open raw file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    RNFail("Unable to open raw file %s", filename);
    return 0;
  }

  // Write raw grid values
  if (!strcmp(format, "Double64") || !strcmp(format, "double64")) {
    if (!WriteRawValues<RNScalar64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Float64") || !strcmp(format, "float64")) {
    if (!WriteRawValues<RNScalar64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Float32") || !strcmp(format, "float32")) {
    if (!WriteRawValues<RNScalar32>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uint64") || !strcmp(format, "Uint64")) {
    if (!WriteRawValues<RNUInt64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Int64") || !strcmp(format, "int64")) {
    if (!WriteRawValues<RNInt64>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uint32") || !strcmp(format, "uint32")) {
    if (!WriteRawValues<RNUInt32>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Int32") || !strcmp(format, "int32")) {
    if (!WriteRawValues<RNInt32>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uint16") || !strcmp(format, "uint16")) {
    if (!WriteRawValues<RNUInt16>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Int16") || !strcmp(format, "int16")) {
    if (!WriteRawValues<RNInt16>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Uchar8") || !strcmp(format, "uchar8")) {
    if (!WriteRawValues<RNUChar8>(fp, grid_values, grid_size)) return 0;
  } else if (!strcmp(format, "Char8") || !strcmp(format, "char8")) {
    if (!WriteRawValues<RNChar8>(fp, grid_values, grid_size)) return 0;
  } else {
    fprintf(stderr, "Unrecognized format %s in %s\n", format, size_name);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////

int R3Grid::
ReadDelphiFile(const char *filename)
{
  // Delphi writes phi files as follows.
  // Each fortran write appears as <leading_record_size> data <trailing_record_size>
  // write(14)'now starting phimap '
  // write(14)nxtlbl,toplbl
  // write(14)phimap
  // write(14)' end of phimap  '
  // write(14)scale,oldmid,igrid
 
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open delphi file %s\n", filename);
    return 0;
  }

  // Read leading record size for uplbl
  int uplbl_record_size;
  if (fread(&uplbl_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read uplbl_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (uplbl_record_size != 20) {
    RNFail("Wrong uplbl record size (%d) in delphi file: %s\n", uplbl_record_size, filename);
    return 0;
  }

  // Read uplbl from file 
  char uplbl[32];
  if (fread(uplbl, sizeof(char), 20, fp) != 20) {
    RNFail("Unable to read uplbl from delphi file: %s\n", filename);
    return 0;
  }

  // Read trailing record size for uplbl
  if (fread(&uplbl_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read uplbl_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (uplbl_record_size != 20) {
    RNFail("Wrong uplbl record size (%d) in delphi file: %s\n", uplbl_record_size, filename);
    return 0;
  }

  // Read leading record size for nextlbl and toplbl
  int nextlbl_record_size;
  if (fread(&nextlbl_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read nextlbl_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (nextlbl_record_size != 70) {
    RNFail("Wrong nextlbl record size (%d) in delphi file: %s\n", nextlbl_record_size, filename);
    return 0;
  }

  // Read nextlbl from file 
  char nextlbl[32];
  if (fread(nextlbl, sizeof(char), 10, fp) != 10) {
    RNFail("Unable to read nextlbl from delphi file: %s\n", filename);
    return 0;
  }

  // Read toplbl from file 
  char toplbl[64];
  if (fread(toplbl, sizeof(char), 60, fp) != 60) {
    RNFail("Unable to read toplbl from delphi file: %s\n", filename);
    return 0;
  }

  // Read trailing record size for nextlbl and toplbl
  if (fread(&nextlbl_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read nextlbl_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (nextlbl_record_size != 70) {
    RNFail("Wrong nextlbl record size (%d) in delphi file: %s\n", nextlbl_record_size, filename);
    return 0;
  }

  // Read leading record size for grid values
  int res = 65;
  unsigned int values_record_size;
  if (fread(&values_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read values_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (values_record_size != res*res*res*sizeof(float)) {
    RNFail("Wrong values record size (%d) in delphi file: %s\n", values_record_size, filename);
    return 0;
  }

  // Re-allocate grid values
  int new_size = res * res * res;
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }
  // Update grid resolution variables
  grid_resolution[0] = res;
  grid_resolution[1] = res;
  grid_resolution[2] = res;
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read grid values
  RNScalar *grid_valuesp = grid_values;
  for (int i = 0; i < grid_size; i++) {
    float value;
    if (fread(&value, sizeof(float), 1, fp) != 1) {
      RNFail("Unable to read grid value %d of %d from file", i, grid_size);
      return 0;
    }
    else {
      *(grid_valuesp++) = (RNScalar) value;
    }
  }

  // Read trailing record size for grid values
  if (fread(&values_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read values_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (values_record_size != res*res*res*sizeof(float)) {
    RNFail("Wrong values record size (%d) in delphi file: %s\n", values_record_size, filename);
    return 0;
  }

  // Read leading record size for trailer text
  int trailer_record_size;
  if (fread(&trailer_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read trailer_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (trailer_record_size != 16) {
    RNFail("Wrong trailer record size (%d) in delphi file: %s\n", trailer_record_size, filename);
    return 0;
  }

  // Read trailer from file 
  char trailer[128];
  if (fread(trailer, sizeof(char), 16, fp) != 16) {
    RNFail("Unable to read trailer from delphi file: %s\n", filename);
    return 0;
  }

  // Read trailing record size for trailer text
  if (fread(&trailer_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read trailer_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (trailer_record_size != 16) {
    RNFail("Wrong trailer record size (%d) in delphi file: %s\n", trailer_record_size, filename);
    return 0;
  }

  // Read leading record size for scale, center, and igrid
  int scale_record_size;
  if (fread(&scale_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read scale_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (scale_record_size != 20) {
    RNFail("Wrong scale record size (%d) in delphi file: %s\n", scale_record_size, filename);
    return 0;
  }

  // Read scale
  float scale;
  if (fread(&scale, sizeof(float), 1, fp) != 1) {
    RNFail("Unable to read scale from delphi file");
    return 0;
  }

  // Read center
  float center[3];
  if (fread(center, sizeof(float), 3, fp) != 3) {
    RNFail("Unable to read center from delphi file");
    return 0;
  }

  // Read grid resolution
  int igrid;
  if (fread(&igrid, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read igrid from delphi file");
    return 0;
  }

  // Check grid resolution
  if (igrid != 65) {
    RNFail("Phi grid file has invalid igrid: %d\n", igrid);
    return 0;
  }

  // Read trailing record size for scale, center, and igrid
  if (fread(&scale_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read scale_record_size from delphi file: %s", filename);
    return 0;
  }
  else if (scale_record_size != 20) {
    RNFail("Wrong scale record size (%d) in delphi file: %s\n", scale_record_size, filename);
    return 0;
  }

  // Make sure at end of file
  int count = 0;
  char dummytrailer;
  while (fread(&dummytrailer, sizeof(char), 1, fp) == 1) count++;
  if (count != 0) {
    fprintf(stderr, "%d bytes of extra data at end of delphi file %s\n", count, filename);
    return 0;
  }

  // Determine world-grid transformation
  R3Point start_point(center[0] - 32 / scale, center[1] - 32 / scale, center[2] - 32 / scale);
  R3Point end_point(center[0] + 32 / scale, center[1] + 32 / scale, center[2] + 32 / scale);
  R3Box world_box(start_point, end_point);
  SetWorldToGridTransformation(world_box);

  // Close file
  fclose(fp);

  // Return number of grid values read
  return grid_size;
}



////////////////////////////////////////////////////////////////////////

struct R3GridCCP4Header {
  int nc, nr, ns; // number of columns (fastest changing), rows, and sections (slowest changing)
  int mode; // 2=real
  int ncstart, nrstart, nsstart;  // negative?
  int nx, ny, nz; // same as nc, nr, ns?
  float xlength, ylength, zlength; 
  float alpha, beta, gamma;
  int mapc, mapr, maps;  // which axis corresponds to each dimension - 1,2,3 for X,Y,Z
  float amin, amax, amean; // min, max, and mean density values
  int ispg; // space group
  int nsymbt; // number of bytes of symmetry info
  int dummy[232];
};



int R3Grid::
ReadCCP4File(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open CCP4 file %s\n", filename);
    return 0;
  }

  // Read header
  R3GridCCP4Header header;
  if (fread(&header, sizeof(R3GridCCP4Header), 1, fp) != 1) {
    fprintf(stderr, "Unable to read header from CCP4 file: %s", filename);
    return 0;
  }

  // Check dimensions
  if ((header.nc <= 0) || (header.nr <= 0) || (header.ns <= 0)) {
    fprintf(stderr, "Invalid dimensions (%d %d %d) in CCP4 file: %s", header.nc, header.nr, header.ns, filename);
    return 0;
  }


  // Check mode
  if (header.mode != 2) { // Reals
    fprintf(stderr, "Invalid mode (%d rather than 2) in CCP4 file: %s", header.mode, filename);
    return 0;
  }

  // Check symmetry information
  if (header.nsymbt < 0) {
    fprintf(stderr, "Invalid number of symmetry bytes (%d) in CCP4 file: %s", header.nsymbt, filename);
    return 0;
  }

#if 0
  // Read symmetry information
  for (int i = 0; i < header.nsymbt; i++) {
    char c;
    if (fread(&c, sizeof(char), 1, fp) != 1) {
      fprintf(stderr, "Unable to read symmetry info from CCP4 file: %s", filename);
      return 0;
    }
  }
#endif

  // Re-allocate grid values
  int new_size = header.nc * header.nr * header.ns;
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }

  // Update grid resolution variables
  grid_resolution[0] = (header.mapc == 1) ? header.nc : ((header.mapr == 1) ? header.nr : header.ns);
  grid_resolution[1] = (header.mapc == 2) ? header.nc : ((header.mapr == 2) ? header.nr : header.ns);
  grid_resolution[2] = (header.mapc == 3) ? header.nc : ((header.mapr == 3) ? header.nr : header.ns);
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read map 
  int x = 0;
  int y = 0;
  int z = 0;
  for (int k = 0; k < header.ns; k++) {
    if (header.maps == 1) x = k;
    else if (header.maps == 2) y = k;
    else if (header.maps == 3) z = k;
    for (int j = 0; j < header.nr; j++) {
      if (header.mapr == 1) x = j;
      else if (header.mapr == 2) y = j;
      else if (header.mapr == 3) z = j;
      for (int i = 0; i < header.nc; i++) {
        if (header.mapc == 1) x = i;
        else if (header.mapc == 2) y = i;
        else if (header.mapc == 3) z = i;

        // Read float
        float value;
        if (fread(&value, sizeof(float), 1, fp) != 1) {
          fprintf(stderr, "Unable to read map from CCP4 file: %s", filename);
          return 0;
        }

        // Store in grid as double
        SetGridValue(x, y, z, value);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////

int R3Grid::
ReadVoxelFile(const char *filename)
{
  // Voxel grid file as defined by Misha Kazhdan
  // One binary int N for resolution followed by NxNxN binary floats 

  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open voxel file %s\n", filename);
    return 0;
  }

  // Read resolution
  int res;
  if (fread(&res, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read resolution from voxel file: %s", filename);
    return 0;
  }
  // Re-allocate grid values
  int new_size = res * res * res;
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }

  // Update grid resolution variables
  grid_resolution[0] = res;
  grid_resolution[1] = res;
  grid_resolution[2] = res;
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read grid values
  RNScalar *grid_valuesp = grid_values;
  for (int k = 0; k < grid_size; k++) {
    float value;
    if (fread(&value, sizeof(float), 1, fp) != 1) {
      RNFail("Unable to read value from voxel file: %s", filename);
      return 0;
    }
    *(grid_valuesp++) = (RNScalar) value;
  }

  // Close file
  fclose(fp);

  // Return number of grid values read
  return grid_size;
}



////////////////////////////////////////////////////////////////////////

int R3Grid::
ReadInsightFile(const char *filename)
{
  // Delphi writes insight files as follows.
  // Each fortran write appears as <leading_record_size> data <trailing_record_size>
  // write(14)toplbl
  // write(14)ivary,nbyte,intdat,extent,extent,extent,
  //   xang,yang,zang,xstart,xend,ystart,yend,zstart,
  //   zend,intx,inty,intz
  // do 9041 k = 1,igrid
  //   do 9042 j = 1,igrid
  //     write(14)(phimap(i,j,k),i=1,igrid)

  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open insight file %s\n", filename);
    return 0;
  }

  // Read leading record size for toplbl
  int toplbl_record_size;
  if (fread(&toplbl_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read toplbl_record_size from insight file: %s", filename);
    return 0;
  }
  else if (toplbl_record_size != 60) {
    RNFail("Wrong toplbl record size (%d) in insight file: %s\n", toplbl_record_size, filename);
    return 0;
  }

  // Read toplbl from file 
  char toplabel[128];
  if (fread(toplabel, sizeof(char), 60, fp) != 60) {
    RNFail("Unable to read toplabel from insight file: %s\n", filename);
    return 0;
  }

  // Read trailing record size for toplbl
  if (fread(&toplbl_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read toplbl_record_size from insight file: %s", filename);
    return 0;
  }
  else if (toplbl_record_size != 60) {
    RNFail("Wrong toplbl record size (%d) in insight file: %s\n", toplbl_record_size, filename);
    return 0;
  }

  // Read leading record size for ivary-intervals
  int ivary_record_size;
  if (fread(&ivary_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read ivary_record_size from insight file: %s", filename);
    return 0;
  }
  else if (ivary_record_size != 72) {
    RNFail("Wrong ivary record size (%d) in insight file: %s\n", ivary_record_size, filename);
    return 0;
  }

  // Read ivary (0 means x varies most quickly)
  int ivary;
  if (fread(&ivary, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read ivary from insight file");
    return 0;
  }

  // Read nbyte (# of bytes in data (per entry?)
  int nbyte;
  if (fread(&nbyte, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read nbyte from insight file");
    return 0;
  }

  // Read intdat (0 means floating point)
  int intdat;
  if (fread(&intdat, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read intdat from insight file");
    return 0;
  }

  // Read extent
  float extent[3];
  if (fread(extent, sizeof(float), 3, fp) != 3) {
    RNFail("Unable to read extent from insight file");
    return 0;
  }

  // Read unit cell angle (90,90,90)
  float unit_cell_angle[3];
  if (fread(unit_cell_angle, sizeof(float), 3, fp) != 3) {
    RNFail("Unable to read cell angle from insight file");
    return 0;
  }

  // Read start and end values
  float start[3], end[3];
  for (int i = 0; i < 3; i++) {
    if (fread(&start[i], sizeof(float), 1, fp) != 1) {
      RNFail("Unable to read start coordinate from insight file");
      return 0;
    }
    if (fread(&end[i], sizeof(float), 1, fp) != 1) {
      RNFail("Unable to read end coordinate from insight file");
      return 0;
    }
  }

  // Read grid resolution
  int intervals[3];
  if (fread(intervals, sizeof(int), 3, fp) != 3) {
    RNFail("Unable to read number of intervals from insight file");
    return 0;
  }

  // Read trailing record size for ivary-intervals
  if (fread(&ivary_record_size, sizeof(int), 1, fp) != 1) {
    RNFail("Unable to read ivary_record_size from insight file: %s", filename);
    return 0;
  }
  else if (ivary_record_size != 72) {
    RNFail("Wrong ivary record size (%d) in insight file: %s\n", ivary_record_size, filename);
    return 0;
  }

  // Determine resolution
  int res[3];
  res[0] = intervals[0]+1;
  res[1] = intervals[1]+1;
  res[2] = intervals[2]+1;

  // Check standard values
  if ((ivary != 0) || (nbyte != 4) || (intdat != 0) ||
      (unit_cell_angle[0] != 90) || (unit_cell_angle[1] != 90) || (unit_cell_angle[2] != 90)) {
    RNFail("Problem with header info in insight file: %d %d %d (should be 0 4 0 90 90 90)\n", 
           ivary, nbyte, intdat, unit_cell_angle[0], unit_cell_angle[1], unit_cell_angle[2]);
    return 0;
  }

#if 0
  printf("%d %d %d   (%g %g %g)   (%g %g %g)    (%g %g %g) (%g %g %g)    (%d %d %d)   (%d %d %d)\n",
    ivary, nbyte, intdat,
    extent[0], extent[1], extent[2],
    unit_cell_angle[0], unit_cell_angle[1], unit_cell_angle[2],
    start[0], start[1], start[2],
    end[0], end[1], end[2],
    intervals[0], intervals[1], intervals[2], 
    res[0], res[1], res[2]);
#endif

  // Re-allocate grid values
  int new_size = res[0] * res[1] * res[2];
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }
  // Update grid resolution variables
  grid_resolution[0] = res[0];
  grid_resolution[1] = res[1];
  grid_resolution[2] = res[2];
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read grid values
  RNScalar *grid_valuesp = grid_values;
  for (int k = 0; k < grid_resolution[2]; k++) {
    for (int j = 0; j < grid_resolution[1]; j++) {
      // Read leading record size for grid row
      unsigned int row_record_size;
      if (fread(&row_record_size, sizeof(int), 1, fp) != 1) {
        RNFail("Unable to read row_record_size from insight file: %s", filename);
        return 0;
      }
      else if (row_record_size != grid_resolution[0]*sizeof(float)) {
        RNFail("Wrong row record size (%d) in insight file: %s\n", row_record_size, filename);
        return 0;
      }

      // Read row
      for (int i = 0; i < grid_resolution[0]; i++) {
        float value;
        if (fread(&value, sizeof(float), 1, fp) != 1) {
          RNFail("Unable to read grid value %d of %d from file", i, grid_size);
          return 0;
        }
        else {
          *(grid_valuesp++) = (RNScalar) value;
        }
      }

      // Read trailing record size for grid row
      if (fread(&row_record_size, sizeof(int), 1, fp) != 1) {
        RNFail("Unable to read row_record_size from insight file: %s", filename);
        return 0;
      }
      else if (row_record_size != grid_resolution[0]*sizeof(float)) {
        RNFail("Wrong row record size (%d) in insight file: %s\n", row_record_size, filename);
        return 0;
      }
    }
  }

  // Make sure at end of file
  int count = 0;
  char dummytrailer;
  while (fread(&dummytrailer, sizeof(char), 1, fp) == 1) count++;
  if (count != 0) {
    fprintf(stderr, "%d bytes of extra data at end of delphi file %s\n", count, filename);
    return 0;
  }

  // Determine world-grid transformation
  R3Point start_point(start[0]*extent[0], start[1]*extent[1], start[2]*extent[2]);
  R3Point end_point(end[0]*extent[0], end[1]*extent[1], end[2]*extent[2]);
  R3Box world_box(start_point, end_point);
  SetWorldToGridTransformation(world_box);

  // Close file
  fclose(fp);

  // Return number of grid values read
  return grid_size;
}



////////////////////////////////////////////////////////////////////////

int R3Grid::
ReadDXFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Skip comments at top
  char buffer[1024];
  do {
    if (!fgets(buffer, 1024, fp)) {
       fprintf(stderr, "Error reading comments in %s\n", filename);
       return 0;
    }
  } while (buffer[0] == '#');

  // Read dimensions
  int res[3];
  if (sscanf(buffer, "object 1 class gridpositions counts %d %d %d", &res[0], &res[1], &res[2]) != 3) {
    fprintf(stderr, "Unable to read dimensions in %s\n", filename);
    return 0;
  }

  // Read origin
  if (!fgets(buffer, 1023, fp)) {
    fprintf(stderr, "Unable to read origin in %s\n", filename);
    return 0;
  }

  // Parse origin
  R3Point origin;
  if (sscanf(buffer, "origin %lf %lf %lf", &origin[0], &origin[1], &origin[2]) != 3) {
    fprintf(stderr, "Unable to parse origin in %s\n", filename);
    return 0;
  }

  // Read delta for each dimension
  R3Vector delta[3];
  for (RNDimension dim = RN_X; dim <= RN_Z; dim++) {
    // Read delta
    if (!fgets(buffer, 1023, fp)) {
      fprintf(stderr, "Unable to read delta in %s\n", filename);
      return 0;
    }

    // Parse delta
    if (sscanf(buffer, "delta %lf %lf %lf", &delta[dim][0], &delta[dim][1], &delta[dim][2]) != 3) {
      fprintf(stderr, "Unable to parse delta in %s\n", filename);
      return 0;
    }
  }

  // Read object 2 line
  if (!fgets(buffer, 1023, fp)) {
    fprintf(stderr, "Unable to read object 2 line in %s\n", filename);
    return 0;
  }

  // Read object 3 line
  if (!fgets(buffer, 1023, fp)) {
    fprintf(stderr, "Unable to read object 3 line in %s\n", filename);
    return 0;
  }

  // Re-allocate grid values
  int new_size = res[0] * res[1] * res[2];
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }

  // Update grid resolution variables
  grid_resolution[0] = res[0];
  grid_resolution[1] = res[1];
  grid_resolution[2] = res[2];
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read data
  RNScalar *grid_valuesp  = grid_values;
  for (int i = 0; i < grid_size; i++) {
    if (fscanf(fp, "%lf", grid_valuesp++) != 1) {
      fprintf(stderr, "Error reading grid values from %s\n", filename);
      return 0;
    }
  }

  // Determine world-grid transformation
  R3Vector extent(delta[0][0] * (res[0]-1), delta[1][1] * (res[1]-1), delta[2][2] * (res[2]-1));
  R3Box world_box(origin, origin + extent);
  SetWorldToGridTransformation(world_box);

  // Return success
  return 1;
}



int R3Grid::
WriteDXFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Print header lines
  const R3Box bbox = WorldBox();
  const RNScalar delta = GridToWorldScaleFactor();
  fprintf(fp, "object 1 class gridpositions counts %d %d %d\n", XResolution(), YResolution(), ZResolution());
  fprintf(fp, "origin %lf %lf %lf\n", bbox[0][0], bbox[0][1], bbox[0][2]);
  fprintf(fp, "delta %lf 0 0\ndelta 0 %lf 0\ndelta 0 0 %lf\n", delta, delta, delta);
  fprintf(fp, "\nobject 2 class gridconnections counts %d %d %d\n", XResolution(), YResolution(), ZResolution());
  fprintf(fp, "\nobject 3 class array type double rank 0 items %d data follows\n", (XResolution() * YResolution() * ZResolution()));

  // Print values
  for (int i = 0; i < XResolution(); i++) {
    for (int j = 0; j < YResolution(); j++) {
      for (int k = 0; k < ZResolution(); k++) {
	fprintf(fp, "%lf\n", GridValue(i, j, k));
      }
    }
  }

  // Print trailer lines
  fprintf(fp, "\nend\n");

  // Close file
  fclose(fp);

  // Return OK status
  return grid_size;
}




int R3Grid::
WritePDBFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Print values
  int n = 0;
  for (int k = 0; k < ZResolution(); k++) {
    for (int j = 0; j < YResolution(); j++) {
      for (int i = 0; i < XResolution(); i++) {
	if (GridValue(i,j, k) > 0){
	  R3Point world_position = WorldPosition(i, j, k);	  
	  fprintf(fp, "%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n", "HETATM", n, "H", ' ', "POC", "1", n, ' ', world_position.X(), world_position.Y(), world_position.Z(), GridValue(i, j, k), GridValue(i, j, k), "");
	  n++;
	}
      }
    }
  }

  // Close file
  fclose(fp);

  // Return OK status
  return grid_size;
}


////////////////////////////////////////////////////////////////////////

int R3Grid::
ReadASCIIFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read resolutions
  int res[3];
  if (fscanf(fp, "%d%d%d\n", &res[0], &res[1], &res[2]) != (unsigned int) 3) {
    fprintf(stderr, "Unable to read file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Read world to grid transformation
  double m[16];
  for (int i = 0; i < 16; i++) {
    if (fscanf(fp, "%lf", &m[i]) != (unsigned int) 1) {
      fprintf(stderr, "Invalid format in %s\n", filename);
      fclose(fp);
      return 0;
    }
  }

  // Re-allocate grid values
  int new_size = res[0] * res[1] * res[2];
  if (!grid_values || (new_size > grid_size)) { 
    if (grid_values) delete [] grid_values;
    grid_values = new RNScalar [ new_size ];
    assert(grid_values);
  }

  // Update grid resolution variables
  grid_resolution[0] = res[0];
  grid_resolution[1] = res[1];
  grid_resolution[2] = res[2];
  grid_row_size = grid_resolution[0];
  grid_sheet_size = grid_row_size * grid_resolution[1];
  grid_size = grid_sheet_size * grid_resolution[2];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    fclose(fp);
    return 0;
  }

  // Read data
  RNScalar *grid_valuesp  = grid_values;
  for (int i = 0; i < grid_size; i++) {
    if (fscanf(fp, "%lf", grid_valuesp++) != 1) {
      fprintf(stderr, "Error reading grid values from %s\n", filename);
      fclose(fp);
      return 0;
    }
  }

  // Set world-grid transformation
  SetWorldToGridTransformation(R3Affine(R4Matrix(m)));

  // Close the file
  fclose(fp);

  // Return success
  return grid_size;
}



int R3Grid::
WriteASCIIFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Print resolutions
  fprintf(fp, "%d %d %d\n", XResolution(), YResolution(), ZResolution());

  // Print world to grid transformation
  const R4Matrix& m = world_to_grid_transform.Matrix();
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      fprintf(fp, "%lf ", m[i][j]);
    }
    fprintf(fp, "\n");
  }

  // Print values
  for (int i = 0; i < XResolution(); i++) {
    for (int j = 0; j < YResolution(); j++) {
      for (int k = 0; k < ZResolution(); k++) {
	fprintf(fp, "%lf\n", GridValue(i, j, k));
      }
    }
  }

  // Close file
  fclose(fp);

  // Return OK status
  return grid_size;
}




int R3Grid::
Print(FILE *fp) const
{
  // Check file
  if (!fp) fp = stdout;

  // Print resolutions
  fprintf(fp, "%d %d %d\n", XResolution(), YResolution(), ZResolution());

  // Print values
  for (int k = 0; k < ZResolution(); k++) {
    for (int j = 0; j < YResolution(); j++) {
      for (int i = 0; i < XResolution(); i++) {
        fprintf(fp, "%g ", GridValue(i, j, k));
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }

  // Return number of grid values written
  return grid_size;
}




int R3Grid::
ConnectedComponents(RNScalar isolevel, int max_components, int *seeds, int *sizes, int *grid_components)
{
  // Allocate array of component identifiers
  int *components = grid_components;
  if (!grid_components) {
    components = new int [ grid_size ];
    assert(components);
  }

  // Initialize array of components
  for (int i = 0; i < grid_size; i++) 
    components[i] = -1;

  // Find connected components
  int ncomponents = 0;
  int next_seed = 0;
  while (TRUE){
    // Find unmarked seed
    int seed = next_seed;
    while (seed < grid_size) { 
      if ((components[seed] < 0) && (grid_values[seed] > isolevel)) break;
      seed++;
    }

    // Check if found a seed
    if (seed >= grid_size) break;

    // Flood fill marking all grid entries 6-connected to seed
    int x, y, z, neighbor;
    int size = 0;
    RNArray<RNScalar *> stack;
    stack.Insert(&grid_values[seed]);
    components[seed] = ncomponents;
    while (!stack.IsEmpty()) {
      // Pop top of stack
      RNScalar *c = stack.Tail();
      stack.RemoveTail();

      // Add grid entry to component
      int index = c - grid_values;
      assert(index >= 0);
      assert(index < grid_size);
      assert(components[index] == ncomponents);
      assert(grid_values[index] > isolevel);

      // Add neighbors to stack
      IndexToIndices(index, x, y, z);
      if (x > 0) {
        IndicesToIndex(x-1, y, z, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (x < XResolution()-1) {
        IndicesToIndex(x+1, y, z, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (y > 0) {
        IndicesToIndex(x, y-1, z, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (y < YResolution()-1) {
        IndicesToIndex(x, y+1, z, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] > isolevel)) { 
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (z > 0) {
        IndicesToIndex(x, y, z-1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (z < ZResolution()-1) {
        IndicesToIndex(x, y, z+1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }

      // Increment component_size
      size++;
    }

    // Update output variables
    if (ncomponents < max_components) {
      if (seeds) seeds[ncomponents] = seed;
      if (sizes) sizes[ncomponents] = size;
    }

    // Increment the number of components
    ncomponents++;
  }


  // Delete components
  if (!grid_components) delete [] components;

  // Return number of connected components found
  return ncomponents;
}



static R3MeshVertex *
InterpolatedVertex(const R3Grid *grid, int ix0, int iy0, int iz0, int dim,
  R3Mesh *mesh, R3MeshVertex **vertices, RNScalar isolevel)
{
  // Check indices
  assert(dim < 3);
  int grid_index;
  int ix1 = (dim == 0) ? ix0+1 : ix0;
  int iy1 = (dim == 1) ? iy0+1 : iy0;
  int iz1 = (dim == 2) ? iz0+1 : iz0;
  assert(ix0 < grid->XResolution());
  assert(iy0 < grid->YResolution());
  assert(iz0 < grid->ZResolution());
  assert(ix1 < grid->XResolution());
  assert(iy1 < grid->YResolution());
  assert(iz1 < grid->ZResolution());
  grid->IndicesToIndex(ix0, iy0, iz0, grid_index);
  if (!vertices[3*grid_index+dim]) {
    RNScalar value0 = grid->GridValue(ix0, iy0, iz0);
    RNScalar value1 = grid->GridValue(ix1, iy1, iz1);
    RNScalar delta0 = fabs(value0 - isolevel);
    RNScalar delta1 = fabs(value1 - isolevel);
    RNScalar denom = delta0 + delta1;
    RNScalar t = (RNIsNotZero(denom)) ? delta0 / denom : 0.5;
    R3Point p0 = grid->WorldPosition(ix0, iy0, iz0);
    R3Point p1 = grid->WorldPosition(ix1, iy1, iz1);
    R3Point p = (1.0-t)*p0 + t*p1;
    vertices[3*grid_index+dim] = mesh->CreateVertex(p);
  }

  // Return vertex
  return vertices[3*grid_index+dim];
}



int R3Grid::
GenerateIsoSurface(RNScalar isolevel, R3Mesh *mesh) const
{
  // Initialize marching cubes edge table
  static int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   
  };

  // Initialize marching cubes triangle table
  static int triTable[256][16] =
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
     {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
     {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
     {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
     {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
     {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
     {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
     {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
     {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
     {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
     {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
     {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
     {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
     {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
     {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
     {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
     {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
     {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
     {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
     {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
     {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
     {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
     {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
     {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
     {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
     {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
     {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
     {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
     {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
     {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
     {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
     {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
     {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
     {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
     {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
     {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
     {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
     {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
     {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
     {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
     {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
     {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
     {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
     {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
     {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
     {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
     {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
     {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
     {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
     {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
     {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
     {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
     {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
     {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
     {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
     {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
     {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
     {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
     {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
     {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
     {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
     {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
     {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
     {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
     {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
     {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
     {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
     {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
     {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
     {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
     {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
     {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
     {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
     {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
     {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
     {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
     {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
     {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
     {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
     {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
     {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
     {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
     {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
     {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
     {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
     {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
     {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
     {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
     {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
     {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
     {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
     {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
     {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
     {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
     {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
     {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
     {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
     {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
     {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
     {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
     {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
     {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
     {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
     {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
     {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
     {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
     {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
     {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
     {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
     {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
     {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
     {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
     {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
     {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
     {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
     {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
     {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
     {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
     {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
     {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
     {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
     {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
     {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
     {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
     {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
     {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
     {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
     {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
     {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
     {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
     {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
     {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
     {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
     {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
     {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
     {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
     {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
     {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
     {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
     {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
     {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
     {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
     {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
     {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
     {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
     {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
     {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
     {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
     {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
     {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
     {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
     {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
     {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
     {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
     {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
     {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
     {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
     {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
     {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
     {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
     {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
     {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
     {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
     {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
     {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
     {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  // Allocate vertex pointers
  R3MeshVertex **vertices = new R3MeshVertex * [3*NEntries()];
  for (int i = 0; i < 3*NEntries(); i++) vertices[i] = NULL;
  
  // Create faces
  RNScalar corner_levels[8];
  for (int iz0 = 0; iz0 < ZResolution()-1; iz0++) {
    for (int iy0 = 0; iy0 < YResolution()-1; iy0++) {
      for (int ix0 = 0; ix0 < XResolution()-1; ix0++) {

        // Compute cube corner values
        corner_levels[0] = GridValue(ix0,iy0,iz0);
        corner_levels[1] = GridValue(ix0+1,iy0,iz0);
        corner_levels[2] = GridValue(ix0+1,iy0,iz0+1);
        corner_levels[3] = GridValue(ix0,iy0,iz0+1);
        corner_levels[4] = GridValue(ix0,iy0+1,iz0);
        corner_levels[5] = GridValue(ix0+1,iy0+1,iz0);
        corner_levels[6] = GridValue(ix0+1,iy0+1,iz0+1);
        corner_levels[7] = GridValue(ix0,iy0+1,iz0+1);

        // Compute cube index
        int cubeindex = 0;
        if (corner_levels[0] < isolevel) cubeindex |= 1;
        if (corner_levels[1] < isolevel) cubeindex |= 2;
        if (corner_levels[2] < isolevel) cubeindex |= 4;
        if (corner_levels[3] < isolevel) cubeindex |= 8;
        if (corner_levels[4] < isolevel) cubeindex |= 16;
        if (corner_levels[5] < isolevel) cubeindex |= 32;
        if (corner_levels[6] < isolevel) cubeindex |= 64;
        if (corner_levels[7] < isolevel) cubeindex |= 128;

        // Check if cube is entirely above/below the surface
        if (edgeTable[cubeindex] == 0) continue;

        // Find the vertices where the surface intersects the cube 
        R3MeshVertex *vertlist[12];
        if (edgeTable[cubeindex] & 1) {
          vertlist[0] = InterpolatedVertex(this, ix0, iy0, iz0, 0, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 2) {
          vertlist[1] = InterpolatedVertex(this, ix0+1, iy0, iz0, 2, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 4) {
          vertlist[2] = InterpolatedVertex(this, ix0, iy0, iz0+1, 0, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 8) {
          vertlist[3] = InterpolatedVertex(this, ix0, iy0, iz0, 2, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 16) {
          vertlist[4] = InterpolatedVertex(this, ix0, iy0+1, iz0, 0, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 32) {
          vertlist[5] = InterpolatedVertex(this, ix0+1, iy0+1, iz0, 2, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 64) {
          vertlist[6] = InterpolatedVertex(this, ix0, iy0+1, iz0+1, 0, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 128) {
          vertlist[7] = InterpolatedVertex(this, ix0, iy0+1, iz0, 2, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 256) {
          vertlist[8] = InterpolatedVertex(this, ix0, iy0, iz0, 1, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 512) {
          vertlist[9] = InterpolatedVertex(this, ix0+1, iy0, iz0, 1, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 1024) {
          vertlist[10] = InterpolatedVertex(this, ix0+1, iy0, iz0+1, 1, mesh, vertices, isolevel);
        }
        if (edgeTable[cubeindex] & 2048) {
          vertlist[11] = InterpolatedVertex(this, ix0, iy0, iz0+1, 1, mesh, vertices, isolevel);
        }

        // Create the triangles
        for (int i = 0; triTable[cubeindex][i] != -1; i+=3) {
          R3MeshVertex *v0 = vertlist[triTable[cubeindex][i  ]];
          R3MeshVertex *v1 = vertlist[triTable[cubeindex][i+1]];
          R3MeshVertex *v2 = vertlist[triTable[cubeindex][i+2]];
          if (!v0 || !v1 || !v2) continue;
          mesh->CreateFace(v0, v1, v2);
        }
      }
    }
  }

  // Delete vertex pointers
  delete [] vertices;

  // Return success
  return 1;
}



int R3Grid::
GenerateIsoSurface(RNScalar isolevel, R3Point *points, int max_points) const
{
  // Initialize pointer into array of points
  R3Point *pointsp = points;

  // Fill array of points with triplets forming triangles of isosurface
  int npoints = 0;
  R3Point corner_points[8];
  RNScalar corner_levels[8];
  for (int k = 0; k < ZResolution()-1; k++) {
    corner_points[0][2] = corner_points[1][2] = corner_points[4][2] = corner_points[5][2] = k;
    corner_points[2][2] = corner_points[3][2] = corner_points[6][2] = corner_points[7][2] = k+1;
    for (int j = 0; j < YResolution()-1; j++) {
      corner_points[0][1] = corner_points[1][1] = corner_points[2][1] = corner_points[3][1] = j;
      corner_points[4][1] = corner_points[5][1] = corner_points[6][1] = corner_points[7][1] = j+1;
      for (int i = 0; i < XResolution()-1; i++) {
        corner_points[0][0] = corner_points[3][0] = corner_points[4][0] = corner_points[7][0] = i;
        corner_points[1][0] = corner_points[2][0] = corner_points[5][0] = corner_points[6][0] = i+1;

        // Check if points array is full (five tris is max for any cell)
        if (npoints >= max_points - 3*5) return npoints;

        // Fill in grid values
        corner_levels[0] = GridValue(i,j,k);
        corner_levels[1] = GridValue(i+1,j,k);
        corner_levels[2] = GridValue(i+1,j,k+1);
        corner_levels[3] = GridValue(i,j,k+1);
        corner_levels[4] = GridValue(i,j+1,k);
        corner_levels[5] = GridValue(i+1,j+1,k);
        corner_levels[6] = GridValue(i+1,j+1,k+1);
        corner_levels[7] = GridValue(i,j+1,k+1);

        // Create triangles
        int n = MarchingCubes(corner_points, corner_levels, isolevel, pointsp);
        pointsp += n;
        npoints += n;
      }
    }
  }

  // Just checking
  assert((npoints % 3) == 0);
  assert(npoints < max_points);

  // Return number of points
  return npoints;
}



void R3Grid::
DrawIsoSurface(RNScalar isolevel) const
{
  // Allocate storage for isosurface
  static const R3Grid *isosurface_grid = NULL;
  static RNScalar isosurface_level = -12345679;
  static const int isosurface_max_points = 8*1024*1024;
  static R3Point *isosurface_points = NULL;
  static int isosurface_npoints = 0;

  // Allocate points for isosurface
  if (!isosurface_points) { 
    isosurface_points = new R3Point [isosurface_max_points];
    assert(isosurface_points);
  }

  // Generate isosurface
  if ((this != isosurface_grid) || (isosurface_level != isolevel)) {
    isosurface_npoints = GenerateIsoSurface(isolevel, isosurface_points, isosurface_max_points);
    isosurface_grid = this;
    isosurface_level = isolevel;
  }

  // Draw isosurface
  R3Point *pointsp = isosurface_points;
  for (int i = 0; i < isosurface_npoints; i += 3) {
    glBegin(GL_POLYGON);
    R3LoadPoint(*(pointsp++));
    R3LoadPoint(*(pointsp++));
    R3LoadPoint(*(pointsp++));
    glEnd();
  }
}



void R3Grid::
DrawSlice(RNDimension dim, int coord) const
{
  // Check coordinates
  if ((dim < 0) || (dim > 2)) return;
  if ((coord < 0) || (coord >= Resolution(dim))) return;

  // Get useful variables
  RNDimension dim1 = (dim+1) % 3;
  RNDimension dim2 = (dim+2) % 3;
  int width = Resolution(dim1);
  int height = Resolution(dim2);
  const int max_resolution = 1024;
  if (width > max_resolution) width = max_resolution;
  if (height > max_resolution) height = max_resolution;

  // Define slice texture
  static const R3Grid *previous_grid[3] = { NULL, NULL, NULL };
  static int previous_coord[3] = { -1, -1, -1 };
  static GLuint texture_id[3] = { 0, 0, 0 };
  if ((this != previous_grid[dim]) || (coord != previous_coord[dim])) {
    if (texture_id[dim] != 0) glDeleteTextures(1, &texture_id[dim]);
    glGenTextures(1, &texture_id[dim]);
    glBindTexture(GL_TEXTURE_2D, texture_id[dim]);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    static GLfloat pixels[max_resolution * max_resolution];
    RNInterval range = Range();
    if (range.Diameter() <= 0) range = RNunit_interval;
    GLfloat *pixelsp = pixels;
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        RNCoord p[3];
        p[dim] = coord;
        p[dim1] = (i+0.5) * Resolution(dim1) / width;
        p[dim2] = (j+0.5) * Resolution(dim2) / height;
        RNScalar value = GridValue(p[0], p[1], p[2]);
        *(pixelsp++) = (GLfloat) ((value - range.Min()) / range.Diameter());
      }
    }
    glTexImage2D(GL_TEXTURE_2D, 0, 1, width, height, 0, GL_LUMINANCE, GL_FLOAT, pixels);
    previous_coord[dim] = coord;
    previous_grid[dim] = this;
  }

  // Set OpenGL modes
  assert(texture_id[dim]);
  glBindTexture(GL_TEXTURE_2D, texture_id[dim]);
  glEnable(GL_TEXTURE_2D);

  // Create quad
  R3Point p0, p1, p2, p3;
  p0[dim] = coord; p0[dim1] = 0;                  p0[dim2] = 0;
  p1[dim] = coord; p1[dim1] = Resolution(dim1)-1; p1[dim2] = 0;
  p2[dim] = coord; p2[dim1] = Resolution(dim1)-1; p2[dim2] = Resolution(dim2)-1;
  p3[dim] = coord; p3[dim1] = 0;                  p3[dim2] = Resolution(dim2)-1;

  // Draw quad 
  RNScalar one = 1.0;
  RNScalar zero = 0;
  R3BeginPolygon();
  R3LoadTextureCoords(zero, zero);
  R3LoadPoint(p0[0], p0[1], p0[2]);
  R3LoadTextureCoords(one, zero);
  R3LoadPoint(p1[0], p1[1], p1[2]);
  R3LoadTextureCoords(one, one);
  R3LoadPoint(p2[0], p2[1], p2[2]);
  R3LoadTextureCoords(zero, one);
  R3LoadPoint(p3[0], p3[1], p3[2]);
  R3EndPolygon();

  // Reset OpenGL modes
  glDisable(GL_TEXTURE_2D);
}





