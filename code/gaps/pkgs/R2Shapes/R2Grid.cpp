// Source file for GAPS scalar grid class



////////////////////////////////////////////////////////////////////////
// NOTE:
// Grid values are defined as samples at the grid positions ranging from
// (0, 0) to (xres-1, yres-1).  Grid values outside this range are undefined.
////////////////////////////////////////////////////////////////////////



// Include files

#include "R2Shapes/R2Shapes.h"



R2Grid::
R2Grid(int xresolution, int yresolution)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_row_size = xresolution;
  grid_size = grid_row_size * yresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  grid_to_world_transform = R2identity_affine;
  world_to_grid_transform = R2identity_affine;
  world_to_grid_scale_factor = 1.0;
}



R2Grid::
R2Grid(int xresolution, int yresolution, const R2Box& bbox)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_row_size = xresolution;
  grid_size = grid_row_size * yresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformation
  SetWorldToGridTransformation(bbox);
}



R2Grid::
R2Grid(int xresolution, int yresolution, const R2Affine& world_to_grid)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_row_size = xresolution;
  grid_size = grid_row_size * yresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformation
  SetWorldToGridTransformation(world_to_grid);
}



R2Grid::
R2Grid(const R2Grid& grid, int x1, int y1, int x2, int y2)
  : grid_values(NULL)
{
  // Determine grid resolution
  grid_resolution[0] = x2 - x1 + 1;
  grid_resolution[1] = y2 - y1 + 1;
  if (grid_resolution[0] < 0) grid_resolution[0] = 0;
  if (grid_resolution[1] < 0) grid_resolution[1] = 0;
  grid_row_size = grid_resolution[0];
  grid_size = grid_row_size * grid_resolution[1];

  // Allocate grid values
  if (grid_size <= 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Copy grid values
  for (int j = 0; j < grid_resolution[1]; j++) {
    int iy = y1 + j;
    for (int i = 0; i < grid_resolution[0]; i++) {
      int ix = x1 + i;
      RNScalar value = 0;
      if ((ix >= 0) && (ix < grid.XResolution())) {
        if ((iy >= 0) && (iy < grid.YResolution())) {
          value = grid.GridValue(ix, iy);
        }
      }
      SetGridValue(i, j, value);
    }
  }

  // Determine transforms
  world_to_grid_transform = grid.world_to_grid_transform;
  world_to_grid_transform.Translate(R2Vector(-x1, -y1));
  world_to_grid_scale_factor = grid.world_to_grid_scale_factor;
  grid_to_world_transform = world_to_grid_transform.Inverse();
}



R2Grid::
R2Grid(const R2Box& bbox, RNLength spacing, int min_resolution, int max_resolution)
{
  // Check for empty bounding box
  if (bbox.IsEmpty() || (RNIsZero(spacing))) { *this = R2Grid(); return; }
  
  // Enforce max resolution
  if (max_resolution > 0) {
    if (bbox.XLength() / spacing > max_resolution) spacing = bbox.XLength() / max_resolution;
    if (bbox.YLength() / spacing > max_resolution) spacing = bbox.YLength() / max_resolution;
  }

  // Compute resolution
  grid_resolution[0] = (int) (bbox.XLength() / spacing + 0.5);
  grid_resolution[1] = (int) (bbox.YLength() / spacing + 0.5);

  // Enforce min resolution
  if (min_resolution > 0) {
    if (grid_resolution[0] < min_resolution) grid_resolution[0] = min_resolution;
    if (grid_resolution[1] < min_resolution) grid_resolution[1] = min_resolution;
  }

  // Set grid resolution
  grid_row_size = grid_resolution[0];
  grid_size = grid_row_size * grid_resolution[1];

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  SetWorldToGridTransformation(bbox);
}



R2Grid::
R2Grid(const R2Grid& grid)
  : grid_values(NULL)
{
  // Copy everything
  *this = grid;
}



R2Grid::
R2Grid(const R2Image& image, int dummy)
  : grid_values(NULL)
{
  // Determine grid resolution
  grid_resolution[0] = image.Width();
  grid_resolution[1] = image.Height();
  grid_row_size = grid_resolution[0];
  grid_size = grid_row_size * grid_resolution[1];

  // Allocate grid values
  if (grid_size <= 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Copy grid values
  for (int j = 0; j < grid_resolution[1]; j++) {
    for (int i = 0; i < grid_resolution[0]; i++) {
      RNRgb rgb = image.PixelRGB(i, j);
      SetGridValue(i, j, rgb.Luminance());
    }
  }

  // Set transformations
  grid_to_world_transform = R2identity_affine;
  world_to_grid_transform = R2identity_affine;
  world_to_grid_scale_factor = 1.0;
}



R2Grid::
~R2Grid(void)
{
  // Deallocate memory for grid values
  if (grid_values) delete [] grid_values;
}



RNInterval R2Grid::
Range(void) const
{
  // Find smallest and largest values
  RNScalar minimum = FLT_MAX;
  RNScalar maximum = -FLT_MAX;
  RNScalar *grid_valuep = grid_values;
  for (int i = 0; i < grid_size; i++) {
    if (*grid_valuep != R2_GRID_UNKNOWN_VALUE) {
      if (*grid_valuep < minimum) minimum = *grid_valuep;
      if (*grid_valuep > maximum) maximum = *grid_valuep;
    }
    grid_valuep++;
  }
  return RNInterval(minimum, maximum);
}



RNScalar R2Grid::
Percentile(RNScalar percentile) const
{
  // Return value at given percentile
  if (grid_size == 0) return 0.0;
  int tmp_count = 0;
  RNScalar *tmp_values = new RNScalar [ grid_size ];
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      tmp_values[tmp_count++] = grid_values[i];
    }
  }
  if (tmp_count == 0) return R2_GRID_UNKNOWN_VALUE;
  qsort(tmp_values, tmp_count, sizeof(RNScalar), RNCompareScalars);
  int percentile_index = percentile * tmp_count;
  RNScalar value = tmp_values[percentile_index];
  delete [] tmp_values;
  return value;
}



RNScalar R2Grid::
L1Norm(void) const
{
  // Return L1 norm of grid
  RNScalar sum = 0.0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      sum += grid_values[i];
    }
  }
  return sum;
}



RNScalar R2Grid::
L2Norm(void) const
{
  // Return L2 norm of grid
  RNScalar sum = 0.0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      RNScalar value = grid_values[i];
      sum += value * value;
    }
  }
  return sqrt(sum);
}



int R2Grid::
Cardinality(void) const
{
  // Return number of non-zero grid values
  int count = 0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      if (grid_values[i] != 0) count++;
    }
  }
  return count;
}



RNScalar R2Grid::
Mean(void) const
{
  // Sum values
  int count = 0;
  RNScalar sum = 0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      sum += grid_values[i];
      count++;
    }
  }

  // Return mean
  if (count == 0) return 0.0;
  else return sum / count;
}



R2Point R2Grid::
GridCentroid(void) const
{
  // Compute weighted sum
  RNScalar total_value = 0;
  R2Point centroid(0,0);
  RNScalar *grid_valuesp = grid_values;
  for (int j = 0; j < grid_resolution[1]; j++) {
    for (int i = 0; i < grid_resolution[0]; i++) {
      R2Vector position(i, j);
      RNScalar value = *(grid_valuesp++);
      centroid += value * position;
      total_value += value;
    }
  }

  // Divide by total value
  if (total_value > 0) centroid /= total_value;

  // Return centroid
  return centroid;
}



R2Diad R2Grid::
GridPrincipleAxes(const R2Point *grid_center, RNScalar *variances) const
{
  // Get centroid
  R2Point center = (grid_center) ? *grid_center : GridCentroid();

  // Compute covariance matrix
  RNScalar m[4] = { 0, 0, 0, 0 };
  RNScalar total_value = 0;
  RNScalar *grid_valuesp = grid_values;
  for (int j = 0; j < grid_resolution[1]; j++) {
    for (int i = 0; i < grid_resolution[0]; i++) {
      R2Point position(i, j);
      RNScalar value = *(grid_valuesp++);
      RNScalar x = position[0] - center[0];
      RNScalar y = position[1] - center[1];
      m[0] += value * x*x;
      m[1] += value * x*y;
      m[2] += value * x*y;
      m[3] += value * y*y;
      total_value += value;
    }
  }

  // Normalize covariance matrix
  if (total_value == 0) return R2xy_diad;
  for (int i = 0; i < 4; i++) m[i] /= total_value;

  // Compute eigenvalues and eigenvectors
  RNScalar U[4];
  RNScalar W[2];
  RNScalar Vt[4];
  RNSvdDecompose(2, 2, m, U, W, Vt);  // m == U . DiagonalMatrix(W) . Vt

  // Copy principle axes into more convenient form
  // W has eigenvalues (greatest to smallest) and Vt has eigenvectors (normalized)
  R2Vector axes[3];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      axes[i][j] = Vt[2*i+j];
    }
  }

  // Find heavier side of first axis
  grid_valuesp = grid_values;
  RNScalar positive_count = 0;
  RNScalar negative_count = 0;
  for (int j = 0; j < grid_resolution[1]; j++) {
    for (int i = 0; i < grid_resolution[0]; i++) {
      R2Point position(i, j);
      RNScalar value = *(grid_valuesp++);
      R2Vector vector = position - center;
      RNScalar dot = axes[0].Dot(vector);
      if (dot > 0.0) positive_count += value;
      else negative_count += value;
    }
  }

  // Set second axis to form orthonormal triad with first one
  if (positive_count < negative_count) axes[0].Flip();
  axes[1] = axes[0]; axes[1].Rotate(RN_PI_OVER_TWO);

  // Just checking
  assert(RNIsEqual(axes[0].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsEqual(axes[1].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsZero(axes[0].Dot(axes[1]), RN_BIG_EPSILON));

  // Return variances (eigenvalues)
  if (variances) {
    variances[0] = W[0];
    variances[1] = W[1];
  }

  // Return triad
  return R2Diad(axes[0], axes[1]);
}



RNScalar R2Grid::
GridValue(RNScalar x, RNScalar y) const
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return 0.0;
  if ((y < 0) || (y > grid_resolution[1]-1)) return 0.0;

  // Bilinear interpolation
  int ix1 = (int) x;
  int iy1 = (int) y;
  int ix2 = ix1 + 1;
  int iy2 = iy1 + 1;
  if (ix2 >= grid_resolution[0]) ix2 = ix1;
  if (iy2 >= grid_resolution[1]) iy2 = iy1;
  RNScalar dx = x - ix1;
  RNScalar dy = y - iy1;
  RNScalar value11 = GridValue(ix1, iy1);
  RNScalar value12 = GridValue(ix1, iy2);
  RNScalar value21 = GridValue(ix2, iy1);
  RNScalar value22 = GridValue(ix2, iy2);
  RNScalar weight11 = (1.0-dx) * (1.0-dy);
  RNScalar weight12 = (1.0-dx) * dy;
  RNScalar weight21 = dx * (1.0-dy);
  RNScalar weight22 = dx * dy;
  RNScalar value = 0;
  RNScalar weight = 0;
  if (value11 != R2_GRID_UNKNOWN_VALUE) { value += weight11 * value11; weight += weight11; }
  if (value12 != R2_GRID_UNKNOWN_VALUE) { value += weight12 * value12; weight += weight12; }
  if (value21 != R2_GRID_UNKNOWN_VALUE) { value += weight21 * value21; weight += weight21; }
  if (value22 != R2_GRID_UNKNOWN_VALUE) { value += weight22 * value22; weight += weight22; }
  if (weight == 0) return R2_GRID_UNKNOWN_VALUE;
  return value / weight;
}



R2Grid& R2Grid::
operator=(const R2Grid& grid) 
{
  // Copy grid resolution
  grid_resolution[0] = grid.grid_resolution[0];
  grid_resolution[1] = grid.grid_resolution[1];
  grid_row_size = grid.grid_row_size;
  grid_size = grid.grid_size;

  // Copy grid values
  if (grid_values) delete [] grid_values;
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);
  for (int i = 0; i < grid_size; i++) {
    grid_values[i] = grid.grid_values[i];
  }

  // Copy transforms
  grid_to_world_transform = grid.grid_to_world_transform;
  world_to_grid_transform = grid.world_to_grid_transform;
  world_to_grid_scale_factor = grid.world_to_grid_scale_factor;

  // Return this
  return *this;
}



void R2Grid::
Abs(void) 
{
  // Take absolute value of every grid value
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (value != R2_GRID_UNKNOWN_VALUE) grid_values[i] = fabs(value);
  }
}



void R2Grid::
Sqrt(void) 
{
  // Take sqrt of every grid value
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (value != R2_GRID_UNKNOWN_VALUE) grid_values[i] = sqrt(value);
  }
}



void R2Grid::
Square(void) 
{
  // Square every grid value
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (value != R2_GRID_UNKNOWN_VALUE) grid_values[i] *= value;
  }
}



void R2Grid::
Negate(void) 
{
  // Square every grid value
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (value != R2_GRID_UNKNOWN_VALUE) grid_values[i] = -value;
  }
}



void R2Grid::
Invert(void) 
{
  // Invert every grid value
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if ((value != R2_GRID_UNKNOWN_VALUE) && (value != 0)) {
      grid_values[i] = 1.0 / value;
    }
  }
}



void R2Grid::
Transpose(void)
{
  // Transpose values
  R2Grid copy(*this);
  int xres = XResolution();
  int yres = YResolution();
  for (int iy = 0; iy < yres; iy++) {
    for (int ix = 0; ix < xres; ix++) {
      RNScalar value = copy.GridValue(xres-1 - ix, yres-1 - iy);
      SetGridValue(ix, iy, value);
    }
  }
}



void R2Grid::
Normalize(void) 
{
  // Scale so that length of "vector" is one
  Divide(L2Norm());
}



void R2Grid::
FillHoles(void) 
{
  // Build Gaussian filter
  RNScalar sigma = 1.5;
  double denom = -2.0 * sigma * sigma;
  RNScalar filter[5][5];
  for (int i = -2; i <= 2; i++) {
    for (int j = -2; j <= 2; j++) {
      filter[i+2][j+2] = exp((i*i+j*j)/denom);
    }
  }

  // Seed queue with border unknown values 
  RNQueue<RNScalar *> queue;
  const RNScalar on_queue_value = -46573822;
  for (int x = 0; x < grid_resolution[0]; x++) {
    for (int y = 0; y < grid_resolution[1]; y++) {
      if (GridValue(x, y) == R2_GRID_UNKNOWN_VALUE) continue;
      if (GridValue(x, y) == on_queue_value) continue;

      // Push neighbors with unknown values onto queue
      if ((x > 0) && (GridValue(x-1, y) == R2_GRID_UNKNOWN_VALUE)) {
        queue.Push(&grid_values[y*grid_row_size+(x-1)]);
        SetGridValue(x-1, y, on_queue_value);
      }
      if ((y > 0) && (GridValue(x, y-1) == R2_GRID_UNKNOWN_VALUE)) {
        queue.Push(&grid_values[(y-1)*grid_row_size+x]);
        SetGridValue(x, y-1, on_queue_value);
      }
      if ((x < grid_resolution[0]-1) && (GridValue(x+1, y) == R2_GRID_UNKNOWN_VALUE)) {
        queue.Push(&grid_values[y*grid_row_size+(x+1)]);
        SetGridValue(x+1, y, on_queue_value);
      }
      if ((y < grid_resolution[1]-1) && (GridValue(x, y+1) == R2_GRID_UNKNOWN_VALUE)) {
        queue.Push(&grid_values[(y+1)*grid_row_size+x]);
        SetGridValue(x, y+1, on_queue_value);
      }
    }
  }

  // Iteratively update border unknown values with blur of immediate neighbors
  while (!queue.IsEmpty()) {
    // Pop grid cell from queue
    RNScalar *valuep = queue.Pop();
    assert(*valuep == on_queue_value);
    int index = valuep - grid_values;
    assert((index >= 0) && (index < grid_size));
    int x, y; IndexToIndices(index, x, y);

    // Update value
    RNScalar sum = 0;
    RNScalar weight = 0;
    for (int i = -2; i <= 2; i++) {
      for (int j = -2; j <= 2; j++) {
        int nx = x + i;
        int ny = y + j;
        if ((nx < 0) || (nx >= grid_resolution[0])) continue;
        if ((ny < 0) || (ny >= grid_resolution[1])) continue;
        RNScalar value = GridValue(nx, ny);
        if (value == R2_GRID_UNKNOWN_VALUE) continue;
        if (value == on_queue_value) continue;
        RNScalar w = filter[i+2][j+2];
        sum += w * value;
        weight += w;
      }
    }

    // Divide by total weight
    if (weight > 0) *valuep = sum / weight;
    else fprintf(stderr, "Zero weight in fill holes\n"); 

    // Push neighbors with unknown values onto queue
    if ((x > 0) && (GridValue(x-1, y) == R2_GRID_UNKNOWN_VALUE)) {
      queue.Push(&grid_values[y*grid_row_size+(x-1)]);
      SetGridValue(x-1, y, on_queue_value);
    }
    if ((y > 0) && (GridValue(x, y-1) == R2_GRID_UNKNOWN_VALUE)) {
      queue.Push(&grid_values[(y-1)*grid_row_size+x]);
      SetGridValue(x, y-1, on_queue_value);
    }
    if ((x < grid_resolution[0]-1) && (GridValue(x+1, y) == R2_GRID_UNKNOWN_VALUE)) {
      queue.Push(&grid_values[y*grid_row_size+(x+1)]);
      SetGridValue(x+1, y, on_queue_value);
    }
    if ((y < grid_resolution[1]-1) && (GridValue(x, y+1) == R2_GRID_UNKNOWN_VALUE)) {
      queue.Push(&grid_values[(y+1)*grid_row_size+x]);
      SetGridValue(x, y+1, on_queue_value);
    }
  }
}



void R2Grid::
FillHoles(int max_hole_size)
{
  // Interpolate vertically
  for (int ix = 0; ix < XResolution(); ix++) {
    int iy0 = -1;
    for (int iy1 = 0; iy1 < YResolution(); iy1++) {
      if (GridValue(ix, iy1) != R2_GRID_UNKNOWN_VALUE) { 
        if ((iy0 >= 0) && (iy0 < iy1-1) && (iy1-iy0 < max_hole_size)) {
          RNScalar value0 = GridValue(ix, iy0);
          if (value0 == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar value1 = GridValue(ix, iy1);
          if (value1 == R2_GRID_UNKNOWN_VALUE) continue;
          for (int iy = iy0+1; iy < iy1; iy++) {
            RNScalar t = (double) (iy - iy0) / (double) (iy1 - iy0);
            RNScalar value = (1-t)*value0 + t*value1;
            SetGridValue(ix, iy, value);
          }
        }
        iy0 = iy1;
      }
    }
  }

  // Interpolate horizontally
  for (int iy = 0; iy < YResolution(); iy++) {
    int ix0 = -1;
    for (int ix1 = 0; ix1 < XResolution(); ix1++) {
      if (GridValue(ix1, iy) != R2_GRID_UNKNOWN_VALUE) { 
        if ((ix0 >= 0) && (ix0 < ix1-1) && (ix1-ix0 < max_hole_size)) {
          RNScalar value0 = GridValue(ix0, iy);
          if (value0 == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar value1 = GridValue(ix1, iy);
          if (value1 == R2_GRID_UNKNOWN_VALUE) continue;
          for (int ix = ix0+1; ix < ix1; ix++) {
            RNScalar t = (double) (ix - ix0) / (double) (ix1 - ix0);
            RNScalar value = (1-t)*value0 + t*value1;
            SetGridValue(ix, iy, value);
          }
        }
        ix0 = ix1;
      }
    }
  }
}



void R2Grid::
Clear(RNScalar value) 
{
  // Set all grid values to value
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    grid_values[i] = value;
  }
}



void R2Grid::
Dilate(RNLength grid_distance) 
{
  // Make copy so that can restore unknown values
  R2Grid copy(*this);

  // Set pixels (to one) within grid_distance from some non-zero pixel
  SquaredDistanceTransform();
  Threshold(grid_distance * grid_distance, 1, 0);

  // Restore unknown values
  for (int i = 0; i < grid_size; i++) {
    if (copy.grid_values[i] == R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = R2_GRID_UNKNOWN_VALUE;
    }
  }
}



void R2Grid::
Erode(RNLength grid_distance) 
{
  // Make copy so that can restore unknown values
  R2Grid copy(*this);

  // Keep only pixels at least distance from some zero pixel
  Threshold(1.0E-20, 1, 0);
  SquaredDistanceTransform();
  Threshold(grid_distance * grid_distance, 0, 1);

  // Restore unknown values
  for (int i = 0; i < grid_size; i++) {
    if (copy.grid_values[i] == R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = R2_GRID_UNKNOWN_VALUE;
    }
  }
}



void R2Grid::
Blur(RNDimension dim, RNLength grid_sigma) 
{
  // Build filter
  RNScalar sigma = grid_sigma;
  int filter_radius = (int) (3 * sigma + 0.5);
  RNScalar *filter = new RNScalar [ filter_radius + 1 ];
  assert(filter);

  // Make buffer for temporary copy of row
  int res = Resolution(dim);
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

  // Convolve grid with filter 
  for (int j = 0; j < Resolution(1-dim); j++) { 
    for (int i = 0; i < Resolution(dim); i++) 
      buffer[i] = (dim == RN_X) ? GridValue(i, j) : GridValue(j, i); 
    for (int i = 0; i < Resolution(dim); i++) { 
      if (buffer[i] == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar sum = 0;
      RNScalar weight = 0;
      RNScalar value = buffer[i];
      if (value != R2_GRID_UNKNOWN_VALUE) {
        sum += filter[0] * value;
        weight += filter[0];
      }
      int nsamples = i;
      if (nsamples > filter_radius) nsamples = filter_radius;
      for (int m = 1; m <= nsamples; m++) {
        RNScalar value = buffer[i - m];;
        if (value != R2_GRID_UNKNOWN_VALUE) {
          sum += filter[m] * value;
          weight += filter[m];
        }
      }
      nsamples = XResolution() - 1 - i;
      if (nsamples > filter_radius) nsamples = filter_radius;
      for (int m = 1; m <= nsamples; m++) {
        RNScalar value = buffer[i + m];
        if (value != R2_GRID_UNKNOWN_VALUE) {
          sum += filter[m] * value;
          weight += filter[m];
        }
      }
      if (weight > 0) {
        RNScalar value = sum / weight;
        if (dim == RN_X) SetGridValue(i, j, value);
        else SetGridValue(j, i, value);
      }
    }
  }

  // Deallocate memory
  delete [] filter;
  delete [] buffer;
}



void R2Grid::
Blur(RNLength grid_sigma) 
{
  // Build filter
  RNScalar sigma = grid_sigma;
  int filter_radius = (int) (3 * sigma + 0.5);
  RNScalar *filter = new RNScalar [ filter_radius + 1 ];
  assert(filter);

  // Make buffer for temporary copy of row
  int res = XResolution();
  if (res < YResolution()) res = YResolution();
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
  for (int j = 0; j < YResolution(); j++) { 
    for (int i = 0; i < XResolution(); i++) 
      buffer[i] = GridValue(i, j); 
    for (int i = 0; i < XResolution(); i++) { 
      if (buffer[i] == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar sum = 0;
      RNScalar weight = 0;
      RNScalar value = buffer[i];
      if (value != R2_GRID_UNKNOWN_VALUE) {
        sum += filter[0] * value;
        weight += filter[0];
      }
      int nsamples = i;
      if (nsamples > filter_radius) nsamples = filter_radius;
      for (int m = 1; m <= nsamples; m++) {
        RNScalar value = buffer[i - m];;
        if (value != R2_GRID_UNKNOWN_VALUE) {
          sum += filter[m] * value;
          weight += filter[m];
        }
      }
      nsamples = XResolution() - 1 - i;
      if (nsamples > filter_radius) nsamples = filter_radius;
      for (int m = 1; m <= nsamples; m++) {
        RNScalar value = buffer[i + m];
        if (value != R2_GRID_UNKNOWN_VALUE) {
          sum += filter[m] * value;
          weight += filter[m];
        }
      }
      if (weight > 0) SetGridValue(i, j, sum / weight);
    }
  }

  // Convolve grid with filter in Y direction
  for (int j = 0; j < XResolution(); j++) { 
    for (int i = 0; i < YResolution(); i++) 
      buffer[i] = GridValue(j, i); 
    for (int i = 0; i < YResolution(); i++) { 
      if (buffer[i] == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar sum = 0;
      RNScalar weight = 0;
      RNScalar value = buffer[i];
      if (value != R2_GRID_UNKNOWN_VALUE) {
        sum += filter[0] * value;
        weight += filter[0];
      }
      int nsamples = i;
      if (nsamples > filter_radius) nsamples = filter_radius;
      for (int m = 1; m <= nsamples; m++) {
        RNScalar value = buffer[i - m];
        if (value != R2_GRID_UNKNOWN_VALUE) {
          sum += filter[m] * value;
          weight += filter[m];
        }
      }
      nsamples = YResolution() - 1 - i;
      if (nsamples > filter_radius) nsamples = filter_radius;
      for (int m = 1; m <= nsamples; m++) {
        RNScalar value = buffer[i + m];
        if (value != R2_GRID_UNKNOWN_VALUE) {
          sum += filter[m] * value;
          weight += filter[m];
        }
      }
      if (weight > 0) SetGridValue(j, i, sum / weight);
    }
  }

  // Deallocate memory
  delete [] filter;
  delete [] buffer;
}



void R2Grid::
AddNoise(RNScalar sigma_fraction)
{
  // Add noise to every grid entry
  if (sigma_fraction == 0) return;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar r = RNRandomScalar();
    RNScalar sign = (RNRandomScalar() < 0.5) ? -1.0 : 1.0;
    RNScalar sigma = sigma_fraction * grid_values[i];
    RNScalar magnitude = sigma * (1.0 - exp(-0.5 * r * r));
    grid_values[i] += sign * magnitude;
  }
}



void R2Grid::
BilateralFilter(RNLength grid_sigma, RNLength value_sigma)
{
  // Make copy of grid
  R2Grid copy(*this);

  // Determine reasonable value sigma
  if (value_sigma == -1) {
    RNInterval range = Range();
    value_sigma = 0.01 * (range.Max() - range.Min());
  }

  // Get convenient variables
  double grid_denom = 2.0 * grid_sigma * grid_sigma;
  double value_denom = 2.0 * value_sigma * value_sigma;
  RNScalar grid_radius = 3 * grid_sigma;
  int r = (int) (grid_radius + 1);
  int r_squared = r * r;

  // Set every sample to be filter of surrounding region in input grid
  for (int cy = 0; cy < YResolution(); cy++) {
    for (int cx = 0; cx < XResolution(); cx++) {
      // Check if current value is unknown - if so, don't update
      RNScalar value = copy.GridValue(cx, cy);
      if (value != R2_GRID_UNKNOWN_VALUE) {
        RNScalar sum = 0;
        RNScalar weight = 0;
        int ymin = cy - r;
        int ymax = cy + r;
        if (ymin < 0) ymin = 0;
        if (ymax >= YResolution()) ymax = YResolution() - 1;
        for (int y = ymin; y <= ymax; y++) {
          int xmin = cx - r;
          int xmax = cx + r;
          if (xmin < 0) xmin = 0;
          if (xmax >= XResolution()) xmax = XResolution() - 1;
          int dy = y - cy;
          for (int x = xmin; x <= xmax; x++) {
            int dx = x - cx;
            int grid_distance_squared = dx*dx + dy*dy;
            if (grid_distance_squared > r_squared) continue;
            RNScalar sample = copy.GridValue(x, y);
            if (sample == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar value_distance_squared = value - sample;
            value_distance_squared *= value_distance_squared;
            RNScalar w = exp(-grid_distance_squared/grid_denom) * exp(-value_distance_squared/value_denom);
            sum += w * sample;
            weight += w;
          }
        }

        // Set grid value
        if (weight == 0) SetGridValue(cx, cy, R2_GRID_UNKNOWN_VALUE);
        else SetGridValue(cx, cy, sum / weight);
      }
    }
  }
}



void R2Grid::
AnisotropicDiffusion(RNLength grid_sigma, RNLength gradient_sigma)
{
  RNAbort("Not implemented");
}



void R2Grid::
PercentileFilter(RNLength grid_radius, RNScalar percentile)
{
  // Make copy of grid
  R2Grid copy(*this);

  // Get convenient variables
  RNScalar grid_radius_squared = grid_radius * grid_radius;
  int r = (int) grid_radius;
  assert(r >= 0);
  int max_samples = (2*r+1) * (2*r+1);
  RNScalar *samples = new RNScalar [ max_samples ];
  assert(samples);

  // Set every sample to be Kth percentile of surrounding region in input grid
  for (int cy = 0; cy < YResolution(); cy++) {
    for (int cx = 0; cx < XResolution(); cx++) {
      // Check if current value is unknown - if so, don't update
      if (GridValue(cx, cy) != R2_GRID_UNKNOWN_VALUE) {
        // Build list of grid values in neighborhood
        int nsamples = 0;
        int ymin = cy - r;
        int ymax = cy + r;
        if (ymin < 0) ymin = 0;
        if (ymax >= YResolution()) ymax = YResolution() - 1;
        for (int y = ymin; y <= ymax; y++) {
          int xmin = cx - r;
          int xmax = cx + r;
          if (xmin < 0) xmin = 0;
          if (xmax >= XResolution()) xmax = XResolution() - 1;
          int dy = y - cy;
          for (int x = xmin; x <= xmax; x++) {
            int dx = x - cx;
            int d_squared = dx*dx + dy*dy;
            if (d_squared > grid_radius_squared) continue;
            RNScalar sample = copy.GridValue(x, y);
            if (sample == R2_GRID_UNKNOWN_VALUE) continue;
            samples[nsamples++] = sample;
          }
        }

        // Check number of grid values in neighborhood
        if (nsamples == 0) {
          SetGridValue(cx, cy, R2_GRID_UNKNOWN_VALUE);
        }
        else {
          // Sort samples found in neighborhood
          qsort(samples, nsamples, sizeof(RNScalar), RNCompareScalars);

          // Set grid value to percentile of neighborhood
          int index = (int) (percentile * nsamples);
          if (index < 0) index = 0;
          else if (index >= nsamples) index = nsamples-1;
          SetGridValue(cx, cy, samples[index]);
        }
      }
    }
  }

  // Delete temporary memory
  delete [] samples;
}



static int 
RNCompareScalarPtrs(const void *value1, const void *value2)
{
  const RNScalar **scalar1pp = (const RNScalar **) value1;
  const RNScalar **scalar2pp = (const RNScalar **) value2;
  const RNScalar *scalar1p = *scalar1pp;
  const RNScalar *scalar2p = *scalar2pp;
  if (*scalar1p < *scalar2p) return -1;
  else if (*scalar1p > *scalar2p) return 1;
  else return 0;
}



static RNBoolean
IsLocalExtremum(const R2Grid& grid, int ix, int iy, RNBoolean maximum)
{
  // Check if local minimum
  RNScalar value = grid.GridValue(ix, iy);
  for (int s = -1; s <= 1; s++) {
    int i = ix + s;
    if ((i < 0) || (i >= grid.XResolution())) continue;
    for (int t = -1; t <= 1; t++) {
      if ((s == 0) && (t == 0)) continue;
      int j = iy + t;
      if ((j < 0) || (j >= grid.YResolution())) continue;
      RNScalar neighbor_value = grid.GridValue(i, j);
      if (!maximum && (neighbor_value < value)) return FALSE;
      if (maximum && (neighbor_value > value)) return FALSE;
    }
  }

  // Passed all tests
  return TRUE;
}



void R2Grid::
MaskNonMinima(RNLength grid_radius)
{
  // Check everything
  if (!GridValues()) return;
  if (NEntries() == 0) return;

  // Determine ammount to add to grid value to mark as masked
  RNInterval range = Range();
  RNScalar mask_summand = 2.0 * range.Diameter();

  // Save copy of grid
  R2Grid copy(*this);
  R2Grid mask(*this);

  // Set all values in this grid to zero
  Clear(R2_GRID_UNKNOWN_VALUE);

  // Load mask grid value pointers into array
  int nptrs = 0;
  RNScalar **ptrs = new RNScalar * [ mask.NEntries() ];
  for (int i = 0; i < mask.NEntries(); i++) {
    if (mask.grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    ptrs[nptrs] = &mask.grid_values[i];
    nptrs++;
  }

  // Check number of non-unknown values
  if (nptrs == 0) { delete [] ptrs; return; }

  // Sort mask value pointers
  qsort(ptrs, nptrs, sizeof(RNScalar *), RNCompareScalarPtrs);

  // Select values in sorted order, masking neighborhoods
  for (int i = 0; i < nptrs; i++){
    int ix, iy;
    RNScalar *ptr = ptrs[i];
    int grid_index = ptr - mask.grid_values;
    mask.IndexToIndices(grid_index, ix, iy);

    // Check if it has been masked
    if (*ptr <= range.Max()) {
      // Set this grid value if it is a local minimum
      if (IsLocalExtremum(copy, ix, iy, FALSE)) {
        SetGridValue(ix, iy, *ptr);
      }

      // Mask neighborhood
      mask.RasterizeGridCircle(R2Point(ix, iy), grid_radius, mask_summand);
    }
  }

  // Delete temporary data
  delete [] ptrs;
}



void R2Grid::
MaskNonMaxima(RNLength grid_radius)
{
  // Check everything
  if (!GridValues()) return;
  if (NEntries() == 0) return;

  // Determine ammount to add to grid value to mark as masked
  RNInterval range = Range();
  RNScalar mask_summand = -2.0 * range.Diameter();

  // Save copy of grid
  R2Grid copy(*this);
  R2Grid mask(*this);

  // Set all values in this grid to zero
  Clear(R2_GRID_UNKNOWN_VALUE);

  // Load mask grid value pointers into array
  int nptrs = 0;
  RNScalar **ptrs = new RNScalar * [ mask.NEntries() ];
  for (int i = 0; i < mask.NEntries(); i++) {
    if (mask.grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    ptrs[nptrs] = &mask.grid_values[i];
    nptrs++;
  }

  // Check number of non-unknown values
  if (nptrs == 0) { delete [] ptrs; return; }

  // Sort mask value pointers
  qsort(ptrs, nptrs, sizeof(RNScalar *), RNCompareScalarPtrs);

  // Select values in sorted order, masking neighborhoods
  for (int i = nptrs-1; i >= 0; i--){
    int ix, iy;
    RNScalar *ptr = ptrs[i];
    int grid_index = ptr - mask.grid_values;
    assert(grid_index < mask.NEntries());
    mask.IndexToIndices(grid_index, ix, iy);

    // Check if it has been masked
    if (*ptr >= range.Min()) {
      // Set this grid value if it is a local maximum
      if (IsLocalExtremum(copy, ix, iy, TRUE)) {
        SetGridValue(ix, iy, *ptr);
      }

      // Mask neighborhood
      mask.RasterizeGridCircle(R2Point(ix, iy), grid_radius, mask_summand);
    }
  }

  // Delete temporary data
  delete [] ptrs;
}



void R2Grid::
Convolve(const RNScalar filter[3][3]) 
{
  // Make temporary copy of grid
  R2Grid copy(*this);

  // Mark boundaries unknown
  for (int i = 0; i < XResolution(); i++) { 
    SetGridValue(i, 0, R2_GRID_UNKNOWN_VALUE);
    SetGridValue(i, YResolution()-1, R2_GRID_UNKNOWN_VALUE);
  }
  for (int j = 0; j < YResolution(); j++) { 
    SetGridValue(0, j, R2_GRID_UNKNOWN_VALUE);
    SetGridValue(XResolution()-1, j, R2_GRID_UNKNOWN_VALUE);
  }

  // Convolve grid with 3x3 filter
  for (int j = 1; j < YResolution()-1; j++) { 
    for (int i = 1; i < XResolution()-1; i++) { 
      RNScalar value = copy.GridValue(i, j);
      if (value != R2_GRID_UNKNOWN_VALUE) {
        RNScalar sum = 0;
        RNBoolean unknown = FALSE;
        for (int dj = -1; dj <= 1; dj++) {
          for (int di = -1; di <= 1; di++) {
            value = copy.GridValue(i + di, j + dj);
            if (value == R2_GRID_UNKNOWN_VALUE) { unknown = TRUE; break; }
            else sum += filter[dj+1][di+1] * value;
          }
          if (unknown) break; 
        }
        if (unknown) SetGridValue(i, j, R2_GRID_UNKNOWN_VALUE);
        else SetGridValue(i, j, sum);
      }
    }
  }
}



void R2Grid::
Gradient(RNDimension dim)
{
  // Set up xfilter
  const RNScalar xfilter[3][3] = { 
    { -0.25, 0, 0.25 }, 
    { -0.5,  0, 0.5}, 
    { -0.25, 0, 0.25 } 
  };
  const RNScalar yfilter[3][3] = { 
    { -0.25, -0.5, -0.25 }, 
    {  0,     0,    0 }, 
    {  0.25,  0.5,  0.25 } 
  };

  // Convolve with filter
  if (dim == RN_X) Convolve(xfilter);
  else if (dim == RN_Y) Convolve(yfilter);
  else RNAbort("Invalid dimension");
}



void R2Grid::
Hessian(RNDimension dim1, RNDimension dim2)
{
  // Compute gradient twice
  Gradient(dim1);
  Gradient(dim2);
}



void R2Grid::
GradientAngle(void)
{
  // Compute direction of gradient
  R2Grid gx(*this); gx.Gradient(RN_X);
  R2Grid gy(*this); gy.Gradient(RN_Y);
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) {
      continue;
    }
    RNScalar x = gx.GridValue(i);
    if (x == R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      continue;
    }
    RNScalar y = gy.GridValue(i);
    if (y == R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      continue;
     }
    if (x == 0) {
      grid_values[i] = RN_PI_OVER_TWO;
      continue;
    }
    grid_values[i] = atan(y/x);
  }
}



void R2Grid::
GradientMagnitude(void)
{
  // Compute magnitude of gradient (Sobel operator)
  R2Grid gx(*this); gx.Gradient(RN_X);
  R2Grid gy(*this); gy.Gradient(RN_Y);
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) {
      continue;
    }
    RNScalar x = gx.GridValue(i);
    if (x == R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      continue;
    }
    RNScalar y = gy.GridValue(i);
    if (y == R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      continue;
    }
    grid_values[i] = sqrt(x*x + y*y);
  }
}



void R2Grid::
Laplacian(void)
{
  // Set up Laplacian filter
  const RNScalar filter[3][3] = { 
    { -0.125, -0.125, -0.125 }, 
    { -0.125,  1,     -0.125 }, 
    { -0.125, -0.125, -0.125 } 
  };

  // Convolve with filter
  Convolve(filter);
}



void R2Grid::
Laplacian(RNDimension dim)
{
  // Set up 1D Laplacian filter
  const RNScalar filter[3] = { -1, 2, -1 };

  // Make temporary copy of grid
  R2Grid copy(*this);

  // Mark boundaries unknown
  for (int i = 0; i < XResolution(); i++) { 
    SetGridValue(i, 0, R2_GRID_UNKNOWN_VALUE);
    SetGridValue(i, YResolution()-1, R2_GRID_UNKNOWN_VALUE);
  }
  for (int j = 0; j < YResolution(); j++) { 
    SetGridValue(0, j, R2_GRID_UNKNOWN_VALUE);
    SetGridValue(XResolution()-1, j, R2_GRID_UNKNOWN_VALUE);
  }

  // Convolve grid with 3x1 filter
  for (int j = 1; j < YResolution()-1; j++) { 
    for (int i = 1; i < XResolution()-1; i++) { 
      RNScalar value = copy.GridValue(i, j);
      if (value != R2_GRID_UNKNOWN_VALUE) {
        RNScalar sum = 0;
        RNBoolean unknown = FALSE;
        for (int k = -1; k <= 1; k++) {
          if (dim == 0) value = copy.GridValue(i + k, j);
          else value = copy.GridValue(i, j + k);
          if (value == R2_GRID_UNKNOWN_VALUE) { unknown = TRUE; break; }
          else sum += filter[k+1] * value;
        }
        if (unknown) SetGridValue(i, j, R2_GRID_UNKNOWN_VALUE);
        else SetGridValue(i, j, sum);
      }
    }
  }
}



void R2Grid::
HarrisCornerFilter(int radius, RNScalar kappa)
{
  // Compute gradients
  R2Grid xgradient(*this);  
  R2Grid ygradient(*this);  
  xgradient.Gradient(RN_X);
  ygradient.Gradient(RN_Y);

  // Initialize this grid
  Clear(R2_GRID_UNKNOWN_VALUE);

  // Get convenient variables
  RNScalar twice_squared_sigma = 0.5 * radius * radius;

  // Create harris response function
  for (int cx = 0; cx < XResolution(); cx++) {
    for (int cy = 0; cy < YResolution(); cy++) {
      // Check this pixel
      if (xgradient.GridValue(cx, cy) == R2_GRID_UNKNOWN_VALUE) continue;
      if (ygradient.GridValue(cx, cy) == R2_GRID_UNKNOWN_VALUE) continue;

      // Compute covariance matrix of window centered at pixel
      RNScalar dx2 = 0; 
      RNScalar dy2 = 0; 
      RNScalar dxdy = 0; 
      for (int s = -radius; s <= radius; s++) {
        int s2 = s*s;
        int ix = cx + s;
        if ((ix < 0) || (ix >= XResolution())) continue;
        for (int t = -radius; t <= radius; t++) {
          int t2 = t*t;
          int iy = cy + t;
          if ((iy < 0) || (iy >= YResolution())) continue;
          RNScalar dx = xgradient.GridValue(ix, iy);
          if (dx == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar dy = ygradient.GridValue(ix, iy);
          if (dy == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar w = exp(-(s2 + t2)/twice_squared_sigma);
          dx2 += w * dx*dx;
          dy2 += w * dy*dy;
          dxdy += w * dx*dy;
        }
      }

      // Compute eigenvalues
      // RNScalar det = dx2*dy2 - dxdy*dxdy;
      // RNScalar trace = dx2 + dy2;
      // RNScalar gap_squared = trace*trace - 4.0*det;
      // if (gap_squared < 0) continue;
      // RNScalar lambda1 = 0.5*(trace + sqrt(gap_squared));
      // RNScalar lambda2 = 0.5*(trace - sqrt(gap_squared));
      
      // Compute harris response from covariance matrix
      RNScalar det = dx2*dy2 - dxdy*dxdy;
      RNScalar trace = dx2 + dy2;
      RNScalar R = det - kappa*trace*trace;

      // Set grid value
      SetGridValue(cx, cy, R);
    }
  }
}



void R2Grid::
Substitute(RNScalar old_value, RNScalar new_value) 
{
  // Replace all instances of old_value with new_value
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == old_value) {
      grid_values[i] = new_value;
    }
  }
}



void R2Grid::
Add(RNScalar value) 
{
  // Add value to all grid values 
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] += value;
    }
  }
}



void R2Grid::
Add(const R2Grid& grid) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Add passed grid values to corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      else grid_values[i] += grid.grid_values[i];
    }
  }
}



void R2Grid::
Subtract(RNScalar value) 
{
  // Add the opposite
  Add(-value);
}



void R2Grid::
Subtract(const R2Grid& grid) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Subtract passed grid values from corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      else grid_values[i] -= grid.grid_values[i];
    }
  }
}



void R2Grid::
Multiply(RNScalar value) 
{
  // Multiply grid values by value
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] *= value;
    }
  }
}



void R2Grid::
Multiply(const R2Grid& grid) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Multiply passed grid values by corresponding entries of this grid
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      else grid_values[i] *= grid.grid_values[i];
    }
  }
}



void R2Grid::
Divide(RNScalar value) 
{
  // Just checking
  if (RNIsZero(value)) return;

  // Multiply by recipricol
  Multiply(1.0 / value);
}



void R2Grid::
Divide(const R2Grid& grid) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Divide entries of this grid by by corresponding entries of passed grid 
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) {
       grid_values[i] = R2_GRID_UNKNOWN_VALUE;
      }
      else if (grid.grid_values[i] == 0) {
        if (grid_values[i] > 0) grid_values[i] = RN_INFINITY;
        else if (grid_values[i] < 0) grid_values[i] = -RN_INFINITY;
      }
      else {
        grid_values[i] /= grid.grid_values[i];
      }
    }
  }
}



void R2Grid::
Pow(RNScalar exponent) 
{
  // Apply exponent to all grid values 
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (value != R2_GRID_UNKNOWN_VALUE) {
      grid_values[i] = pow(value, exponent);
    }
  }
}



void R2Grid::
Mask(const R2Grid& grid) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Mask entries of grid
  for (int i = 0; i < grid_size; i++) {
    if (grid.grid_values[i] == 0) grid_values[i] = 0;
    else if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
  }
}



void R2Grid::
Overlay(const R2Grid& grid) 
{
  // Resolutions must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Assign non-zero entries of grid 
  for (int i = 0; i < grid_size; i++) {
    if (grid.grid_values[i] == 0) continue;
    if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    grid_values[i] = grid.grid_values[i];
  }
}



void R2Grid::
Threshold(RNScalar threshold, RNScalar low, RNScalar high) 
{
  // Set grid value to low (high) if less/equal (greater) than threshold
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] != R2_GRID_UNKNOWN_VALUE) {
      if (grid_values[i] <= threshold) {
        if (low != R2_GRID_KEEP_VALUE) grid_values[i] = low;
      }
      else {
        if (high != R2_GRID_KEEP_VALUE) grid_values[i] = high;
      }
    }
  }
}



void R2Grid::
Threshold(const R2Grid& threshold, RNScalar low, RNScalar high) 
{
  // Set grid value to low (high) if less/equal (greater) than threshold
  for (int i = 0; i < grid_size; i++) {
    if ((grid_values[i] != R2_GRID_UNKNOWN_VALUE) && (threshold.grid_values[i] != R2_GRID_UNKNOWN_VALUE)) {
      if (grid_values[i] <= threshold.grid_values[i]) {
        if (low == R2_GRID_INPUT_VALUE) grid_values [i]= threshold.grid_values[i];
        else if (low != R2_GRID_KEEP_VALUE) grid_values[i] = low;
      }
      else {
        if (high == R2_GRID_INPUT_VALUE) grid_values[i] = threshold.grid_values[i];
        else if (high != R2_GRID_KEEP_VALUE) grid_values[i] = high;
      }
    }
  }
}



void R2Grid::
SignedDistanceTransform(void)
{
  // Compute distance from boundary into interior (negative) and into exterior (positive)
  R2Grid copy(*this);
  SquaredDistanceTransform();
  Sqrt();
  copy.Threshold(0, 1, 0);
  copy.Substitute(R2_GRID_UNKNOWN_VALUE, 1);
  copy.SquaredDistanceTransform();
  copy.Sqrt();
  Subtract(copy);
}



void R2Grid::
SquaredDistanceTransform(void)
{
  int x,y,s,t;
  int dist,square,new_dist;
  int* oldBuffer;
  int* newBuffer;
  int first;
  int i;

  // Allocate temporary buffers
  int res = XResolution();
  if (res < YResolution()) res = YResolution();
  oldBuffer = new int[res];
  assert(oldBuffer);
  newBuffer = new int[res];
  assert(newBuffer);

  // Initalize values (0 if was set, max_value if not)
  RNScalar max_value = 2 * (res+1) * (res+1);
  RNScalar *grid_valuesp = grid_values;
  for (i = 0; i < grid_size; i++) {
    if (*grid_valuesp == 0.0) *grid_valuesp = max_value;
    else if (*grid_valuesp == R2_GRID_UNKNOWN_VALUE) *grid_valuesp = max_value;
    else *grid_valuesp = 0.0;
    grid_valuesp++;
  }

  // Scan along x axis
  for (y = 0; y < YResolution(); y++) {
    first = 1;
    dist = 0;
    for (x = 0; x < XResolution(); x++) {
      if (GridValue(x,y) == 0.0) {
        dist=0;
        first=0;
        SetGridValue(x, y, 0);
      }
      else if (first == 0) {
        dist++;
        square = dist*dist;
        SetGridValue(x, y, square);
      }
    }
      		
    // backward scan
    dist = 0;
    first = 1;
    for (x = XResolution()-1; x >= 0; x--) {
      if (GridValue(x,y) == 0.0) {
        dist = 0;
        first = 0;
        SetGridValue(x, y, 0);
      }
      else if (first == 0) {
        dist++;
        square = dist*dist;
        if (square < GridValue(x, y)) {
          SetGridValue(x, y, square);
        }
      }
    }
  }

  // Scan along y axis
  for (x = 0; x < XResolution(); x++) {
    // Copy grid values
    for (y = 0; y < YResolution(); y++) 
      oldBuffer[y] = (int) (GridValue(x, y) + 0.5);
      		
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
        SetGridValue(x, y, dist);
      }
      else { 
        s = y; 
      }
    }
  }
	
  // Delete temporary buffers
  delete[] oldBuffer;
  delete[] newBuffer;
}



void R2Grid::
Voronoi(R2Grid *squared_distance_grid)
{
  int dist;
  int* old_dist;
  int* new_dist;
  RNScalar value;
  RNScalar *old_value;
  RNScalar *new_value;
  R2Grid *dgrid;
  int res,square,tmp_dist,first;
  int x,y,s,t,i;

  // Allocate distance grid
  if (squared_distance_grid) dgrid = squared_distance_grid;
  else dgrid = new R2Grid(XResolution(), YResolution());
  assert(dgrid);
  dgrid->SetWorldToGridTransformation(WorldToGridTransformation());

  // Allocate temporary buffers
  res = XResolution();
  if (res < YResolution()) res = YResolution();
  old_dist = new int[res];
  new_dist = new int[res];
  old_value = new RNScalar[res];
  new_value = new RNScalar[res];
  assert(old_dist && new_dist && old_value && new_value);

  // Initalize distance grid values (0 if was set, max_value if not)
  RNScalar max_value = 3 * (res+1) * (res+1);
  for (i = 0; i < grid_size; i++) {
    if (grid_values[i] == 0.0) dgrid->grid_values[i] = max_value;
    else dgrid->grid_values[i] = 0.0;
  }

  // Scan along x axis
  for (y = 0; y < YResolution(); y++) {
    first = 1;
    value = 0;
    dist = 0;
    for (x = 0; x < XResolution(); x++) {
      if (dgrid->GridValue(x,y) == 0.0) {
        first=0;
        dist=0;
        value = GridValue(x,y);
        assert(value != 0);
      }
      else if (first == 0) {
        dist++;
        square = dist*dist;
        dgrid->SetGridValue(x, y, square);
        SetGridValue(x, y, value);
        assert(value != 0);
      }
    }
			
    // backward scan
    dist = 0;
    first = 1;
    for (x = XResolution()-1; x >= 0; x--) {
      if (dgrid->GridValue(x,y) == 0.0){
        dist = 0;
        first = 0;
        value = GridValue(x,y);
        assert(value != 0);
      }
      else if (first == 0) {
        dist++;
        square = dist*dist;
        if (square < dgrid->GridValue(x, y)) {
          dgrid->SetGridValue(x, y, square);
          SetGridValue(x, y, value);
          assert(value != 0);
        }
      }
    }
  }

  // Scan along y axis
  for (x = 0; x < XResolution(); x++) {
    // Copy grid values
    for (y = 0; y < YResolution(); y++) {
      old_dist[y] = (int) (dgrid->GridValue(x, y));
      old_value[y] = GridValue(x, y);
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
        dgrid->SetGridValue(x, y, dist);
        SetGridValue(x, y, value);
      }
      else { 
        s = y; 
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



void R2Grid::
PointSymmetryTransform(int radius)
{
  // Make copy of grid
  R2Grid copy(*this);
  copy.Normalize();

  // Brute force for now
  for (int j = 0; j < grid_resolution[1]; j++) {
    int y_radius = (j < grid_resolution[1]/2) ? j : grid_resolution[1]-1 - j;
    if ((radius > 0) && (y_radius > radius)) y_radius = radius;
    for (int i = 0; i < grid_resolution[0]; i++) {
      int x_radius = (i < grid_resolution[0]/2) ? i : grid_resolution[0]-1 - i;
      if ((radius > 0) && (x_radius > radius)) x_radius = radius;
      RNScalar dot = 0;
      RNScalar norm = 0;
      for (int t = -y_radius; t <= y_radius; t++) {
        for (int s = -x_radius; s <= x_radius; s++) {
          RNScalar value1 = GridValue(i-s, j-t);
          RNScalar value2 = GridValue(i+s, j+t);
          dot += 2 * value1 * value2;
          norm += value1 * value1; 
          norm += value2 * value2; 
        }
      }
      if (norm > 0) dot /= norm;
      SetGridValue(i, j, dot);
    }
  }
}



void R2Grid::
Gauss(RNLength sigma, RNBoolean square)
{
  // Check sigma
  if (sigma == 0.0) return;    

  // Replace each grid value with Gaussian of what was there before
  RNScalar fac = 1.0 / (2.0 * sigma * sigma);
  RNScalar denom = -2.0 * sigma * sigma;
  if (RNIsZero(denom, 1.0E-6)) return;
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (square) value *= value;
    grid_values[i] = fac * exp( value / denom );
  }
}



void R2Grid::
PadWithZero(int xresolution, int yresolution)
{
  // Add zeros to achieve desired resolution
  if ((XResolution() >= xresolution) && (YResolution() >= yresolution)) return;

  // Copy this grid
  R2Grid copy(*this);

  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_row_size = xresolution;
  grid_size = grid_row_size * yresolution;

  // Allocate grid values
  if (grid_values) delete [] grid_values;
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new RNScalar [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (int i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Copy original values
  for (int iy = 0; iy < copy.YResolution(); iy++) {
    for (int ix = 0; ix < copy.XResolution(); ix++) {
      RNScalar value = copy.GridValue(ix, iy);
      SetGridValue(ix, iy, value);
    }
  }
}



void R2Grid::
Resample(int xresolution, int yresolution)
{
  // Resample grid values at new resolution
  RNScalar *new_grid_values = NULL;
  int new_grid_size = xresolution * yresolution;
  if (new_grid_size > 0) {
    new_grid_values = new RNScalar [ new_grid_size ];
    assert(new_grid_values);
    if (grid_values && (grid_resolution[0] > 0) && (grid_resolution[1] > 0)) {
      RNScalar *new_grid_valuesp = new_grid_values;
      RNScalar xscale = (RNScalar) (grid_resolution[0]-1) / (RNScalar) (xresolution - 1);
      RNScalar yscale = (RNScalar) (grid_resolution[1]-1) / (RNScalar) (yresolution - 1);
      for (int j = 0; j < yresolution; j++) {
        RNScalar y = (j == yresolution-1) ? grid_resolution[1]-1 : j * yscale;
        for (int i = 0; i < xresolution; i++) {
          RNScalar x = (i == xresolution-1) ? grid_resolution[0]-1 : i * xscale;
          *(new_grid_valuesp++) = GridValue(x, y);
        }
      }
    }
    else {
      for (int i = 0; i < new_grid_size; i++) {
        new_grid_values[i] = 0;
      }
    }
  }

  // Reset grid variables
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_row_size = xresolution;
  grid_size = grid_row_size * yresolution;
  if (grid_values) delete [] grid_values;
  grid_values = new_grid_values;

  // Reset transformations
  SetWorldToGridTransformation(WorldBox());
}



void R2Grid::
RasterizeGridValue(int ix, int iy, RNScalar value, int operation)
{
  // Check if within bounds
  if ((ix < 0) || (ix > grid_resolution[0]-1)) return;
  if ((iy < 0) || (iy > grid_resolution[1]-1)) return;

  // Update grid based on operation
  if (operation == R2_GRID_ADD_OPERATION) AddGridValue(ix, iy, value);
  else if (operation == R2_GRID_SUBTRACT_OPERATION) AddGridValue(ix, iy, -value);
  else if (operation == R2_GRID_REPLACE_OPERATION) SetGridValue(ix, iy, value);
  else RNAbort("Unrecognized grid rasterization operation\n");
}



void R2Grid::
RasterizeGridPoint(RNScalar x, RNScalar y, RNScalar value, int operation)
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return;
  if ((y < 0) || (y > grid_resolution[1]-1)) return;

  // Bilinear interpolation
  int ix1 = (int) x;
  int iy1 = (int) y;
  int ix2 = ix1 + 1;
  int iy2 = iy1 + 1;
  if (ix2 >= grid_resolution[0]) ix2 = ix1;
  if (iy2 >= grid_resolution[1]) iy2 = iy1;
  RNScalar dx = x - ix1;
  RNScalar dy = y - iy1;
  RasterizeGridValue(ix1, iy1, value * (1.0-dx) * (1.0-dy), operation);
  RasterizeGridValue(ix1, iy2, value * (1.0-dx) * dy, operation);
  RasterizeGridValue(ix2, iy1, value * dx * (1.0-dy), operation);
  RasterizeGridValue(ix2, iy2, value * dx * dy, operation);
}



void R2Grid::
RasterizeGridSpan(const int p1[2], const int p2[2], RNScalar value1, RNScalar value2, int operation)
{
  // Resolve values
  if (value1 == R2_GRID_UNKNOWN_VALUE) return;
  if (value2 == R2_GRID_UNKNOWN_VALUE) return;

  // Get some convenient variables
  int d[2],p[2],dd[2],s[2];
  for (int i = 0; i < 2; i++) {
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
  int i2=(i1+1)%2;

  // Check span extent
  if(dd[i1]==0){
    // Span is a point - rasterize it
    RNScalar value = 0.5 * (value1 + value2);
    RasterizeGridValue(p[0], p[1], value, operation);
  }
  else {
    // Step along span
    int off[2] = { 0, 0 };
    RNScalar value = value1;
    RNScalar dvalue = (value2 - value1) / dd[i1];
    for (int i = 0; i <= dd[i1]; i++) {
      RasterizeGridValue(p[0], p[1], value, operation);
      off[i2]+=dd[i2];
      p[i1]+=s[i1];
      p[i2]+=s[i2]*off[i2]/dd[i1];
      off[i2]%=dd[i1];
      value += dvalue;
    }
  }
}



void R2Grid::
RasterizeGridBox(const int p1[2], const int p2[2], RNScalar value, int operation)
{
  // Check value
  if (value == R2_GRID_UNKNOWN_VALUE) return;

  // Clip box
  int x1 = p1[0];
  if (x1 < 0) x1 = 0;
  if (x1 > XResolution()-1) return;
  int x2 = p2[0];
  if (x2 < 0) return;
  if (x2 > XResolution()-1) x2 = XResolution()-1;
  int y1 = p1[1];
  if (y1 < 0) y1 = 0;
  if (y1 > YResolution()-1) return;
  int y2 = p2[1];
  if (y2 < 0) return;
  if (y2 > YResolution()-1) y2 = YResolution()-1;

  // Rasterize box
  for (int ix = x1; ix <= x2; ix++) {
    for (int iy = y1; iy <= y2; iy++) {
      RasterizeGridValue(ix, iy, value, operation);
    }
  }
}



void R2Grid::
RasterizeGridTriangle(const int p1[2], const int p2[2], const int p3[2], RNScalar valueA, RNScalar valueB, RNScalar valueC, int operation)
{
  // Resolve values
  if (valueA == R2_GRID_UNKNOWN_VALUE) return;
  if (valueB == R2_GRID_UNKNOWN_VALUE) return;
  if (valueC == R2_GRID_UNKNOWN_VALUE) return;

  // Get convenient point variables
  R2Point points[3];
  points[0][0] = p1[0];
  points[0][1] = p1[1];
  points[1][0] = p2[0];
  points[1][1] = p2[1];
  points[2][0] = p3[0];
  points[2][1] = p3[1];

  // Get convenient value variables
  RNScalar values[3];
  values[0] = valueA;
  values[1] = valueB;
  values[2] = valueC;

  // Sort vertex indices by Y coordinate
  int iv0, iv1, iv2;
  if (points[0].Y() < points[1].Y()) {
    if (points[0].Y() < points[2].Y()) { 
      if (points[1].Y() < points[2].Y()) { iv0 = 0; iv1 = 1; iv2 = 2; }
      else { iv0 = 0; iv1 = 2; iv2 = 1; }
    }
    else { iv0 = 2; iv1 = 0; iv2 = 1; }
  }
  else {
    if (points[1].Y() < points[2].Y()) { 
      if (points[0].Y() < points[2].Y()) { iv0 = 1; iv1 = 0; iv2 = 2; }
      else { iv0 = 1; iv1 = 2; iv2 = 0; }
    }
    else { iv0 = 2; iv1 = 1; iv2 = 0; }
  }
  
  // Sort vertex coordinates
  double dy = points[iv2].Y() - points[iv0].Y();
  double t1 = (dy > 0) ? (points[iv1].Y() - points[iv0].Y()) / dy : 0;
  double y0 = points[iv0].Y();
  double y1 = points[iv1].Y();
  double y2 = points[iv2].Y();
  double x0 = points[iv0].X();
  double x2 = points[iv2].X();
  double x1a = points[iv1].X();
  double x1b = (1-t1)*x0 + t1*x2;
  double value0 = values[iv0];
  double value2 = values[iv2];
  double value1a = values[iv1];
  double value1b = (1-t1)*value0 + t1*value2;
  if (x1a > x1b) { 
    double swap;
    swap = x1a; x1a = x1b; x1b = swap; 
    swap = value1a; value1a = value1b; value1b = swap; 
  }

  // Rasterize lower half of triangle
  int iy0 = (int) (y0 + 0.5);
  int iy1 = (int) (y1 + 0.5);
  int nysteps = (iy1 - iy0) + 1;
  double xa_step = (x1a - x0) / nysteps;
  double xb_step = (x1b - x0) / nysteps;
  double valuea_step = (value1a  - value0) / nysteps;
  double valueb_step = (value1b  - value0) / nysteps;
  double xa = x0;
  double xb = x0;
  double valuea = value0;
  double valueb = value0;
  for (int iy = iy0; iy <= iy1; iy++) {
    int ixa = (int) (xa + 0.5);
    int ixb = (int) (xb + 0.5);
    int nxsteps = (ixb - ixa) + 1;
    double value_step = (valueb - valuea) / nxsteps;
    double value = valuea;
    for (int ix = ixa; ix <= ixb; ix++) {
      RasterizeGridValue(ix, iy, value, operation);
      value += value_step;
    }
    valuea += valuea_step;
    valueb += valueb_step;
    xa += xa_step;
    xb += xb_step;
  }
  
  // Rasterize upper half of triangle
  int iy2 = (int) (y2 + 0.5);
  nysteps = (iy2 - iy1) + 1;
  xa_step = (x2 - x1a) / nysteps;
  xb_step = (x2 - x1b) / nysteps;
  valuea_step = (value2  - value1a) / nysteps;
  valueb_step = (value2  - value1b) / nysteps;
  xa = x1a;
  xb = x1b;
  valuea = value1a;
  valueb = value1b;
  for (int iy = iy1; iy <= iy2; iy++) {
    int ixa = (int) (xa + 0.5);
    int ixb = (int) (xb + 0.5);
    int nxsteps = (ixb - ixa) + 1;
    double value_step = (valueb - valuea) / nxsteps;
    double value = valuea;
    for (int ix = ixa; ix <= ixb; ix++) {
      RasterizeGridValue(ix, iy, value, operation);
      value += value_step;
    }
    valuea += valuea_step;
    valueb += valueb_step;
    xa += xa_step;
    xb += xb_step;
  }
}



void R2Grid::
RasterizeGridCircle(const R2Point& center, RNLength radius, RNScalar value, int operation)
{
  // Figure out the min and max in each dimension
  int mn[2], mx[2];
  for (int i = 0; i < 2; i++) {
    mx[i]= (int) (center[i]+radius);
    if (mx[i] < 0) return;
    if (mx[i] > Resolution(i)-1) mx[i] = Resolution(i)-1;
    mn[i]= (int) (center[i]-radius);
    if (mn[i] > Resolution(i)-1) return;
    if (mn[i] < 0) mn[i] = 0;
  }

  // Rasterize circle interior
  int y1 = (int) (center[1] - radius + 0.5);
  int y2 = (int) (center[1] + radius + 0.5);
  if (y1 < mn[1]) y1 = mn[1];
  if (y2 > mx[1]) y2 = mx[1];
  RNLength radius_squared = radius * radius;
  for (int j = y1; j <= y2; j++) {
    RNCoord y = j - center[1];
    RNCoord y_squared = y*y;
    RNLength x_squared = radius_squared - y_squared;
    RNLength x = sqrt(x_squared);
    int x1 = (int) (center[0] - x + 0.5);
    int x2 = (int) (center[0] + x + 0.5);
    if (x1 < mn[0]) x1 = mn[0];
    if (x2 > mx[0]) x2 = mx[0];
    for (int i = x1; i <= x2; i++) {
      RasterizeGridValue(i, j, value, operation);
    }
  }
}


// Used by RasterizeGridPolygon
struct ScanLineCrossing { RNCoord x; int side; };



// Used by RasterizeGridPolygon
static int CompareScanLineCrossings(const void *data1, const void *data2)
{
  ScanLineCrossing *crossing1 = (ScanLineCrossing *) data1;
  ScanLineCrossing *crossing2 = (ScanLineCrossing *) data2;
  if (RNIsLess(crossing1->x, crossing2->x)) return -1;
  else if (RNIsGreater(crossing1->x, crossing2->x)) return 1;
  else if (crossing1->side < crossing2->side) return -1;
  else if (crossing1->side > crossing2->side) return 1;
  else return 0;
}



void R2Grid::
RasterizeGridPolygon(const R2Polygon& polygon, RNScalar value, int operation) 
{
  // Clip polygon to grid bounding box
  R2Polygon clipped_polygon(polygon);
  clipped_polygon.Clip(GridBox());
  if (clipped_polygon.NPoints() == 0) return; 

  // Check if only one point
  if (clipped_polygon.NPoints() == 1) { 
    RasterizeGridPoint(clipped_polygon.Point(0), value, operation); 
    return; 
  }

  // Check if only two points
  if (clipped_polygon.NPoints() == 2) { 
    RasterizeGridSpan(clipped_polygon.Point(0), clipped_polygon.Point(1), value, operation); 
    return; 
  }

  // Allocate list of scan line crossings
  int max_crossings = 2 * YResolution() * clipped_polygon.NPoints();
  ScanLineCrossing *crossings = new ScanLineCrossing [max_crossings];
  int ncrossings = 0;

  // Scan convert polygon
  for (int iy = 0; iy < YResolution(); iy++) {
    // Build list of scan line crossings
    ncrossings = 0;    
    const R2Point *p1 = &clipped_polygon.Point(clipped_polygon.NPoints()-1);
    RNScalar d1 = p1->Y() - iy;
    for (int i = 0; i < clipped_polygon.NPoints(); i++) {
      if (ncrossings >= max_crossings) break;
      const R2Point *p2 = &clipped_polygon.Point(i);
      RNScalar d2 = p2->Y() - iy;
      if ((d1 < 0) && (d2 >= 0)) {
        RNScalar t = d2 / (d2 - d1);
        crossings[ncrossings].x = p2->X() + t * (p1->X() - p2->X());
        crossings[ncrossings].side = -1;
        ncrossings++;
      }
      else if ((d2 < 0) && (d1 >= 0)) {
        RNScalar t = d1 / (d1 - d2);
        crossings[ncrossings].x = p1->X() + t * (p2->X() - p1->X());
        crossings[ncrossings].side = 1;
        ncrossings++;
      }
      p1 = p2;
      d1 = d2;
    }

    // Check number of crossings
    if (ncrossings < 2) continue;

    // Sort scan line crossings by x coordinate
    qsort(crossings, ncrossings, sizeof(ScanLineCrossing), CompareScanLineCrossings);

    // Rasterize scan line
    int index = 0;
    while (index < ncrossings-1) {
      // Find two crossings
      ScanLineCrossing *crossing1 = &crossings[index++];
      ScanLineCrossing *crossing2 = &crossings[index++];
      while ((index < ncrossings) && (crossing2->side == crossing1->side)) crossing2 = &crossings[index++];
      if (crossing1->side == crossing2->side) break;

      // Fill pixels between crossings
      int ix1 = (int) (crossing1->x + 0.5);
      int ix2 = (int) (crossing2->x + 0.5);
      for (int ix = ix1; ix <= ix2; ix++) {
        RasterizeGridValue(ix, iy, value, operation);
      }
    }
  }

  // Delete crossings
  delete [] crossings;
}



void R2Grid::
RasterizeWorldPolygon(const R2Polygon& polygon, RNScalar value, int operation) 
{
  // Rasterize polygon into grid coordinates
  R2Polygon transformed_polygon(polygon);
  transformed_polygon.Transform(WorldToGridTransformation());
  RasterizeGridPolygon(transformed_polygon, value, operation);
}



RNScalar R2Grid::
Dot(const R2Grid& grid) const
{
  // Resolutions and transforms must be the same (for now)
  assert(grid_resolution[0] == grid.grid_resolution[0]);
  assert(grid_resolution[1] == grid.grid_resolution[1]);

  // Compute dot product between this and grid
  RNScalar dot = 0.0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    dot += grid_values[i] * grid.grid_values[i];
  }
  return dot;
}



RNScalar R2Grid::
L1Distance(const R2Grid& grid) const
{
  // Compute distance between this and grid
  RNScalar distance = 0.0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    distance += fabs(grid_values[i] - grid.grid_values[i]);
  }
  return distance;
}



RNScalar R2Grid::
L2DistanceSquared(const R2Grid& grid) const
{
  // Compute distance between this and grid
  RNScalar distance_squared = 0.0;
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    if (grid.grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar delta = (grid_values[i] - grid.grid_values[i]);
    distance_squared += delta * delta;
  }

  // Return result
  return distance_squared;
}



void R2Grid::
SetWorldToGridTransformation(const R2Affine& affine)
{
  // Set transformations
  world_to_grid_transform = affine;
  grid_to_world_transform = affine.Inverse();
  world_to_grid_scale_factor = affine.ScaleFactor();
}



void R2Grid::
SetWorldToGridTransformation(const R2Box& world_box)
{
  // Just checking
  if (grid_size == 0) return;
  if (world_box.NDimensions() < 2) return;

  // Compute grid origin
  R2Vector grid_diagonal(XResolution()-1, YResolution()-1);
  R2Vector grid_origin = 0.5 * grid_diagonal;

  // Compute world origin
  R2Vector world_diagonal(world_box.XLength(), world_box.YLength());
  R2Vector world_origin = world_box.Centroid().Vector();

  // Compute scale
  RNScalar scale = FLT_MAX;
  RNScalar xscale = (world_diagonal[0] > 0) ? grid_diagonal[0] / world_diagonal[0] : FLT_MAX;
  if (xscale < scale) scale = xscale;
  RNScalar yscale = (world_diagonal[1] > 0) ? grid_diagonal[1] / world_diagonal[1] : FLT_MAX;
  if (yscale < scale) scale = yscale;
  if (scale == FLT_MAX) scale = 1;

  // Compute world-to-grid transformation
  R2Affine affine(R2identity_affine);
  affine.Translate(grid_origin);
  if (scale != 1) affine.Scale(scale);
  affine.Translate(-world_origin);

  // Set transformations
  SetWorldToGridTransformation(affine);
}



void R2Grid::
SetWorldToGridTransformation(const R2Point& world_origin, const R2Vector& world_xaxis, RNLength world_radius)
{
  // Just checking
  if (grid_size == 0) return;

  // Compute grid origin
  R2Vector grid_diagonal(XResolution()-1, YResolution()-1);
  R2Vector grid_origin = 0.5 * grid_diagonal;
  RNScalar grid_radius = grid_origin[0];
  if (grid_origin[1] < grid_radius) grid_radius = grid_origin[1];

  // Compute scale
  if (RNIsZero(world_radius)) return;
  if (RNIsZero(grid_radius)) return;
  RNScalar scale = grid_radius / world_radius;

  // Compute rotation
  RNAngle rotation = R2InteriorAngle(world_xaxis, R2posx_vector);
  if (world_xaxis.Y() < 0.0) rotation = RN_TWO_PI - rotation;

  // Compute world-to-grid transformation
  R2Affine affine(R2identity_affine);
  affine.Translate(grid_origin);
  affine.Rotate(-rotation);
  affine.Scale(scale);
  affine.Translate(-(world_origin.Vector()));

  // Set transformations
  SetWorldToGridTransformation(affine);
}



R2Point R2Grid::
WorldPosition(RNCoord x, RNCoord y) const
{
  // Transform point from grid coordinates to world coordinates
  R2Point world_point(x, y);
  world_point.Transform(grid_to_world_transform);
  return world_point;

}



R2Point R2Grid::
GridPosition(RNCoord x, RNCoord y) const
{
  // Transform point from world coordinates to grid coordinates
  R2Point grid_point(x, y);
  grid_point.Transform(world_to_grid_transform);
  return grid_point;
}



void R2Grid::
Capture(void)
{
  // Check image size
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  assert(grid_resolution[0] >= viewport[2]);
  assert(grid_resolution[1] >= viewport[3]);
  assert(grid_values);

  // Allocate pixels
  float *pixels = new float [ grid_size ];

  // Read pixels from frame buffer 
  glReadPixels(0, 0, grid_resolution[0], grid_resolution[1], GL_LUMINANCE, GL_FLOAT, pixels); 

  // Copy pixels
  for (int i = 0; i < grid_size; i++) grid_values[i] = pixels[i];

  // Delete pixels
  delete [] pixels;
}



void R2Grid::
DrawMesh(void) const
{
  // Push transformation
  grid_to_world_transform.Push();

  // Draw triangle strip
  for (int j = 1; j < grid_resolution[1]; j++) {
    glBegin(GL_TRIANGLE_STRIP);
    for (int i = 0; i < grid_resolution[0]; i++) {
      R3LoadPoint(i, j, GridValue(i, j));
      R3LoadPoint(i, j-1, GridValue(i, j-1));
    }
    glEnd();
  }

  // Pop transformation
  grid_to_world_transform.Pop();
}



void R2Grid::
DrawImage(int x, int y) const
{
  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, grid_resolution[0], 0, grid_resolution[1]);

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Set position for image
  glRasterPos2i(x, y);

  // Allocate pixels
  float *pixels = new float [ grid_size ];

  // Copy pixels
  for (int i = 0; i < grid_size; i++) pixels[i] = (float) grid_values[i];

  // Draw pixels
  glDrawPixels(grid_resolution[0], grid_resolution[1], GL_LUMINANCE, GL_FLOAT, pixels); 

  // Delete pixels
  delete [] pixels;

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}



RNScalar R2Grid::
GridValue(RNScalar x, RNScalar y, RNLength sigma) const
{
  // Check if within bounds
  if ((x < 0) || (x > grid_resolution[0]-1)) return 0.0;
  if ((y < 0) || (y > grid_resolution[1]-1)) return 0.0;

  // Determine Gaussian filter extent
  int f = (int) (3 * sigma + 0.5);
  int xmin = (int) (x - f + 0.5);
  int xmax = (int) (x + f + 0.5);
  int ymin = (int) (y - f + 0.5);
  int ymax = (int) (y + f + 0.5);
  if (xmin < 0) xmin = 0;
  if (xmax > XResolution()-1) xmax = XResolution()-1;
  if (ymin < 0) ymin = 0;
  if (ymax > YResolution()-1) ymax = YResolution()-1;

  // Filter sample with Gaussian
  RNScalar sum = 0.0;
  RNScalar weight = 0.0;
  double denom = 2.0 * sigma * sigma;
  for (int j = ymin; j <= ymax; j++) {
    for (int i = xmin; i <= xmax; i++) {
      RNScalar value = GridValue(i, j);
      if (value != R2_GRID_UNKNOWN_VALUE) { 
        RNScalar dx = x - i;
        RNScalar dy = y - j;
        RNScalar d = dx*dx + dy*dy;
        RNScalar w = exp(-d / denom);
        sum += w * value;
        weight += w;
      }
    }
  }

  // Normalize value based on total weight of Gaussian samples
  if (weight > 0) return sum / weight;
  else return R2_GRID_UNKNOWN_VALUE;
}



void R2Grid::
ConnectedComponentLabelFilter(RNScalar isolevel)
{
  // Compute connected components
  int *components = new int [ grid_size ];
  int ncomponents = ConnectedComponents(isolevel, grid_size, NULL, NULL, components);

  // Assign component labels to grid elements
  if (ncomponents > 0) {
    // Assign component labels
    for (int i = 0; i < grid_size; i++) {
      if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      grid_values[i] = components[i];
    }
  }

  // Delete temporary memory
  delete [] components;
}



void R2Grid::
ConnectedComponentSizeFilter(RNScalar isolevel)
{
  // Compute connected components
  int *sizes = new int [ grid_size ];
  int *components = new int [ grid_size ];
  int ncomponents = ConnectedComponents(isolevel, grid_size, NULL, sizes, components);

  // Assign component sizes to grid elements
  if (ncomponents > 0) {
    for (int i = 0; i < grid_size; i++) {
      if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      int component = components[i];
      if ((component < 0) || (component >= ncomponents)) grid_values[i] = 0;
      else grid_values[i] = sizes[component];
    }
  }

  // Delete temporary memory
  delete [] components;
  delete [] sizes;
}



void R2Grid::
ConnectedComponentCentroidFilter(RNScalar isolevel)
{
  // Compute connected components
  int *sizes = new int [ grid_size ];
  int *components = new int [ grid_size ];
  int ncomponents = ConnectedComponents(isolevel, grid_size, NULL, sizes, components);
  if (ncomponents == 0) { Clear(0); return; }

  // Allocate and initalize centroids
  int *centroid_counts = new int [ ncomponents ];
  R2Point *centroid_positions = new R2Point [ ncomponents ];
  RNScalar *centroid_values = new RNScalar [ ncomponents ];
  for (int i = 0; i < ncomponents; i++) {
    centroid_counts[i] = 0;
    centroid_positions[i] = R2zero_point;
    centroid_values[i] = 0;
  }

  // Sum positions in every connected component
  for (int i = 0; i < grid_size; i++) {
    if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
    int component = components[i];
    if ((component < 0) || (component >= ncomponents)) continue;
    int ix, iy;
    IndexToIndices(i, ix, iy);
    centroid_counts[component]++;
    centroid_positions[component][0] += ix;
    centroid_positions[component][1] += iy;
  }

  // Divide by counts
  for (int i = 0; i < ncomponents; i++) {
    if (centroid_counts[i] == 0) continue;
    centroid_positions[i] /= centroid_counts[i];
  }

  // Remember original values at centroids in grid
  for (int i = 0; i < ncomponents; i++) {
    if (centroid_counts[i] == 0) continue;
    int ix = (int) (centroid_positions[i][0] + 0.5);
    int iy = (int) (centroid_positions[i][1] + 0.5);
    if ((ix < 0) || (ix >= grid_resolution[0])) continue;
    if ((iy < 0) || (iy >= grid_resolution[1])) continue;
    centroid_values[i] = GridValue(ix, iy);
  }

  // Clear out entire grid
  Clear(0);

  // Set values at centroids of connected components
  for (int i = 0; i < ncomponents; i++) {
    if (centroid_counts[i] == 0) continue;
    int ix = (int) (centroid_positions[i][0] + 0.5);
    int iy = (int) (centroid_positions[i][1] + 0.5);
    if ((ix < 0) || (ix >= grid_resolution[0])) continue;
    if ((iy < 0) || (iy >= grid_resolution[1])) continue;
    SetGridValue(ix, iy, centroid_values[i]);
  }

  // Delete temporary memory
  delete [] centroid_counts;
  delete [] centroid_positions;
  delete [] centroid_values;
  delete [] components;
  delete [] sizes;
}



void R2Grid::
ConnectedComponentFilter(RNScalar isolevel, RNArea min_grid_area, RNArea max_grid_area, 
  RNScalar under_isolevel_value, RNScalar too_small_value, RNScalar too_large_value)
{
  // Compute connected components
  int *sizes = new int [ grid_size ];
  int *components = new int [ grid_size ];
  int ncomponents = ConnectedComponents(isolevel, grid_size, NULL, sizes, components);

  // Mask values under isolevel
  if (under_isolevel_value != R2_GRID_KEEP_VALUE) {
    for (int i = 0; i < grid_size; i++) {
      if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      if (grid_values[i] < isolevel) grid_values[i] = under_isolevel_value;
    }
  }

  // Mask components with areas too small 
  if (too_small_value != R2_GRID_KEEP_VALUE) {
    for (int component = 0; component < ncomponents; component++) {
      int area = sizes[component];
      if (area < min_grid_area) {
        for (int i = 0; i < grid_size; i++) {
          if (components[i] == component) grid_values[i] = too_small_value;
        }
      }
    }
  }

  // Mask components with areas too large 
  if (too_large_value != R2_GRID_KEEP_VALUE) {
    for (int component = 0; component < ncomponents; component++) {
      int area = sizes[component];
      if (area > max_grid_area) {
        for (int i = 0; i < grid_size; i++) {
          if (components[i] == component) grid_values[i] = too_large_value;
        }
      }
    }
  }

  // Delete temporary memory
  delete [] components;
  delete [] sizes;
}



int R2Grid::
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
      if ((components[seed] < 0) && (grid_values[seed] != R2_GRID_UNKNOWN_VALUE)&& (grid_values[seed] > isolevel)) break;
      seed++;
    }

    // Check if found a seed
    if (seed >= grid_size) break;

    // Flood fill marking all grid entries 8-connected to seed
    int x, y, neighbor;
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
      IndexToIndices(index, x, y);
      if (x > 0) {
        IndicesToIndex(x-1, y, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (x < XResolution()-1) {
        IndicesToIndex(x+1, y, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE)&& (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (y > 0) {
        IndicesToIndex(x, y-1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE)&& (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if (y < YResolution()-1) {
        IndicesToIndex(x, y+1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE)&& (grid_values[neighbor] > isolevel)) { 
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if ((x > 0) && (y > 0)) {
        IndicesToIndex(x-1, y-1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if ((x > 0) && (y < YResolution()-1)) {
        IndicesToIndex(x-1, y+1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if ((x < XResolution()-1) && (y > 0)) {
        IndicesToIndex(x+1, y-1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE) && (grid_values[neighbor] > isolevel)) {
          components[neighbor] = ncomponents;
          stack.Insert(&grid_values[neighbor]);
        }
      }
      if ((x < XResolution()-1) && (y < YResolution()-1)) {
        IndicesToIndex(x+1, y+1, neighbor);
        if ((components[neighbor] < 0) && (grid_values[neighbor] != R2_GRID_UNKNOWN_VALUE) && (grid_values[neighbor] > isolevel)) {
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



#if 0

static int
NextIndexCCW(const R2Grid *grid, int cur, int prev)
{
  // Get next index in CCW direction (used for generating isocontour)
  int next, next_x, next_y;
  int cur_x, cur_y, prev_x, prev_y;
  IndexToIndices(cur, cur_x, cur_y);
  IndexToIndices(prev, prev_x, prev_y);
  int dx = cur_x - prev_x;
  int dy = cur_y - prev_y;
  if (dx == 1) { next_x = cur_x; next_y = cur_y + 1; } 
  else if (dx == -1) { next_x = cur_x; next_y = cur_y - 1; } 
  else if (dy == 1) { next_x = cur_x - 1; next_y = cur_y; } 
  else if (dy == -1) { next_x = cur_x + 1; next_y = cur_y; } 
  else { next_x = cur_x + 1; next_y = cur_y; } 
  if ((next_x < 0) || (next_x >= grid->XResolution())) return -1;
  if ((next_y < 0) || (next_y >= grid->YResolution())) return -1;
  IndicesToIndex(next_x, next_y, next);
  return next;
}



struct ContourVertex {
  ContourVertex(int x, int y, int dim, int dir, RNScalar t) : x(x), y(y), dim(dim), dir(dir), t(t), mark(0) {};
  int x, y;
  int dim;
  int dir; 
  RNScalar t;
  int mark;
};



int R2Grid::
GenerateIsoContour(RNScalar isolevel, R2Point *points, int max_points) const
{
  // Initialize array of all vertices
  ContourVertex **vertices = new ContourVertex * [ 2 * grid_size ];
  for (int i = 0; i < 2 * grid_size; i++) vertices[i] = NULL;

  // Create vertices on horizontal edges
  for (int j = 0; j < grid_resolution[1]; j++) {
    for (int i = 0; i < grid_resolution[0]-1; i++) {
      RNScalar value1 = GridValue(i, j);
      RNScalar value2 = GridValue(i+1, j);
      RNScalar delta1 = isolevel - value1;
      RNScalar delta2 = isolevel - value2;
      if (delta1 < 0) {
        if (delta2 > 0) {
          RNScalar t = -delta1 / (delta2 - delta1);
          ContourVertex *vertex = new ContourVertex(i+1,j,RN_X,-1,t);
          vertices[j*grid_resolution[0]+i] = vertex;
        }
      }
      else if (delta1 > 0) {
        if (delta2 < 0) {
          RNScalar t = delta1 / (delta1 - delta2);
          ContourVertex *vertex = new ContourVertex(i,j,RN_X,1,t);
          vertices[j*grid_resolution[0]+i] = vertex;
        }
      }
    }
  }

  // Create vertices on vertical edges
  ContourVertex **vverts = new ContourVertex * [ grid_size ];
  for (int i = 0; i < grid_resolution[0]; i++) {
    for (int j = 0; j < grid_resolution[1]-1; j++) {
      RNScalar value1 = GridValue(i, j);
      RNScalar value2 = GridValue(i, j+1);
      RNScalar delta1 = isolevel - value1;
      RNScalar delta2 = isolevel - value2;
      if (delta1 < 0) {
        if (delta2 > 0) {
          RNScalar t = -delta1 / (delta2 - delta1);
          ContourVertex *vertex = new ContourVertex(i,j+1,RN_Y,-1,t);
          vertices[grid_size + j*grid_resolution[0]+i] = vertex;
        }
      }
      else if (delta1 > 0) {
        if (delta2 < 0) {
          RNScalar t = delta1 / (delta1 - delta2);
          ContourVertex *vertex = new ContourVertex(i,j,RN_Y,1,t);
          vertices[grid_size + j*grid_resolution[0]+i] = vertex;
        }
      }
    }
  }

  // Trace contour boundaries
  for (int start = 0; start < 2 * grid_size; start++) {
    if (vertices[start] == NULL) continue;
    if (vertices[start]->mark) continue;

    // Insert start vertex 
    RNArray<ContourVertex *> contour;
    ContourVertex *start_vertex = vertices[start];
    contour.Insert(start_vertex);
    start_vertex->mark = 1;

    // Traverse boundary
    ContourVertex *vertex = start_vertex;
    do {
      // Get convenient variables
      int x = vertex->x;
      int y = vertex->y;
      int dim = vertex->dim;
      int dir = vertex->dir;

      // Get neighbor vertex
      ContourVertex *neighbor = NULL;
      if (dim == RN_X) {
        if (dir == -1) {
          assert((GridValue(x-1,y) < isolevel) && (GridValue(x,y) > isolevel));
          neighbor1 = vertices[ grid_size + (y
        }
        else {
        }
      }
      int x1 = (vertex->dim == RN_X) ? vertex->x + vertex->dir : vertex->x;
      int y1 = (vertex->dim == RN_Y) ? vertex->y + vertex->dir : vertex->y;



    } while(vertex != vertices[start]);
          
    int prev = marks[start];
    for (int i = 0; i < 4; i++) {
      int index = NextIndexCCW(this, start, prev);
      if (index < 0) break;
      if (components[index] != c) break; 
      prev = index;
    }

    // Check if found previous neighbor
    if (prev == -1) continue;
    if (prev == start) continue;

    // Walk along boundary counter-clockwise outputing points
    do {
      int neighbor = prev;
      for (int i = 0; i < 4; i++) {
        neighbor = NextIndexCCW(this, current, neighbor);
        if (components[index] != c) {
          // Output a point
        }
        else {
          // Move to next grid cell
          prev = current;
          current = neighbor;
          break;
        }
      }
    } while (current != start);
  }

  // Return number of points
  return npoints;
}

#endif



int R2Grid::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *input_extension;
  if (!(input_extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .pfm).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".raw", 4)) return ReadRAWFile(filename);
  else if (!strncmp(input_extension, ".pfm", 4)) return ReadPFMFile(filename);
  else if (!strncmp(input_extension, ".pnm", 4)) return ReadPNMFile(filename);
  else if (!strncmp(input_extension, ".png", 4)) return ReadPNGFile(filename);
  else if (!strncmp(input_extension, ".grd", 4)) return ReadGridFile(filename);
  else  return ReadImage(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Grid::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *input_extension;
  if (!(input_extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Output file has no extension (e.g., .pfm).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".raw", 4)) return WriteRAWFile(filename);
  else if (!strncmp(input_extension, ".pfm", 4)) return WritePFMFile(filename);
  else if (!strncmp(input_extension, ".png", 4)) return WritePNGFile(filename);
  else if (!strncmp(input_extension, ".grd", 4)) return WriteGridFile(filename);
  else return WriteImage(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// RAW FORMAT READ/WRITE
////////////////////////////////////////////////////////////////////////

int R2Grid::
ReadRAWFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open raw image file %s", filename);
    return 0;
  }

  // Read header
  unsigned int header[4];
  if (fread(header, sizeof(unsigned int), 4, fp) != 4) {
    fprintf(stderr, "Unable to header to raw image file %s", filename);
    return -1;
  }

  // Check magic number
  const unsigned int magic = 54321;
  if (header[0] != magic) {
    fprintf(stderr, "Invalid header in raw image file %s", filename);
    return -1;
  }

  // Check number of channels
  if (header[3] != 1) {
    fprintf(stderr, "Invalid number of channels (%d)in raw image file %s", header[3], filename);
    return -1;
  }

  // Parse header
  int width = header[1];
  int height = header[2];

  // Allocate pixels
  unsigned int npixels = width * height;
  float *pixels = new float [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate data for raw image file %s", filename);
    return -1;
  }

  // Read pixels
  int count = npixels;
  while (count > 0) {
    int n = fread(&pixels[npixels-count], sizeof(float), count, fp);
    if (n < 0) { fprintf(stderr, "Unable to read pixels from %s\n", filename); return 0; }
    count -= n;
  }

  // Fill in grid info
  grid_resolution[0] = width;
  grid_resolution[1] = height;
  grid_row_size = width;
  grid_size = npixels;
  world_to_grid_scale_factor = 1;
  world_to_grid_transform = R2identity_affine;
  grid_to_world_transform = R2identity_affine;
  if (grid_values) delete [] grid_values;
  grid_values = new RNScalar [ grid_size ];
  for (int i = 0; i < grid_size; i++) {
    if (RNIsEqual(pixels[i], R2_GRID_UNKNOWN_VALUE)) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
    else grid_values[i] = pixels[i];
  }

  // Delete pixels
  delete [] pixels;

  // Close file
  fclose(fp);

  // Return success
  return grid_size;
}


int R2Grid::
WriteRAWFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open raw image file %s", filename);
    return 0;
  }

  // Write header
  const unsigned int magic = 54321;
  unsigned int xres = (unsigned int) grid_resolution[0];
  unsigned int yres = (unsigned int) grid_resolution[1];
  unsigned int header[4] = { magic, xres, yres, 1 };
  if (fwrite(header, sizeof(unsigned int), 4, fp) != 4) {
    fprintf(stderr, "Unable to header to raw image file %s", filename);
    return -1;
  }

  // Allocate pixels
  float *pixels = new float [ grid_row_size ];

  // Write pixels (row by row to avoid large buffers)
  for (int i = 0; i < grid_resolution[1]; i++) {
    RNScalar *valuesp = &grid_values[i * grid_resolution[0]];
    for (int i = 0; i < grid_resolution[0]; i++) pixels[i] = (float) valuesp[i];
    if (fwrite(pixels, sizeof(float), grid_resolution[0], fp) != (unsigned int) grid_resolution[0]) {
      fprintf(stderr, "Unable to write grid values to file %s\n", filename);
      return 0;
    }
  }

  // Delete pixels
  delete [] pixels;

  // Close file
  fclose(fp);

  // Return success
  return grid_size;
}




////////////////////////////////////////////////////////////////////////
// PNM/PNG FORMAT READ/WRITE
////////////////////////////////////////////////////////////////////////

int R2Grid::
ReadPNMFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image: %s\n", filename);
    return 0;
  }

  // Read magic header
  int c;
  c = fgetc(fp); if (c != 'P') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }
  c = fgetc(fp); if (c != '5') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }
  c = fgetc(fp); if (c != '\n') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }

  // Read width
  int width_count = 0;
  char width_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (width_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { width_string[width_count] = '\0'; break; }
    else if (!isdigit(c)) { fprintf(stderr, "Bad width character %c in %s\n", c, filename); return 0; }
    else width_string[width_count++] = c;
  }

  // Check width
  if ((width_count == 0) || (width_count > 128)) {
    fprintf(stderr, "Error reading width in %s\n", filename); 
    return 0; 
  }

  // Read height
  int height_count = 0;
  char height_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (height_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { height_string[height_count] = '\0'; break; }
    else if (!isdigit(c)) { fprintf(stderr, "Bad height character %c in %s\n", c, filename); return 0; }
    else height_string[height_count++] = c;
  }

  // Check height
  if ((height_count == 0) || (height_count > 128)) {
    fprintf(stderr, "Error reading height in %s\n", filename); 
    return 0; 
  }

  // Read max_value
  int max_value_count = 0;
  char max_value_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (max_value_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { max_value_string[max_value_count] = '\0'; break; }
    if (!isdigit(c) && (c != '.') && (c != '-')) { fprintf(stderr, "Bad max_value character %c in %s\n", c, filename); return 0; }
    max_value_string[max_value_count++] = c;
  }

  // Check max_value
  if ((max_value_count == 0) || (max_value_count > 128)) {
    fprintf(stderr, "Error reading max_value in %s\n", filename); 
    return 0; 
  }

  // Parse values
  int width = atoi(width_string);
  int height = atoi(height_string);
  unsigned int max_value = atoi(max_value_string);

  // Fill in grid info
  grid_resolution[0] = width;
  grid_resolution[1] = height;
  grid_row_size = width;
  grid_size = width*height;
  world_to_grid_scale_factor = 1;
  world_to_grid_transform = R2identity_affine;
  grid_to_world_transform = R2identity_affine;
  if (grid_values) delete [] grid_values;
  grid_values = new RNScalar [ grid_size ];
  if (!grid_values) {
    fprintf(stderr, "Unable to allocate %d pixels for %s\n", grid_size, filename);
    return 0;
  }

  // Check endian of this computer
  // PNM files are big endian
  const unsigned short endian_test = 1;
  RNBoolean little_endian = (*((unsigned char *) &endian_test) == '\0') ? FALSE : TRUE;
  printf("Endian = %d\n", little_endian);
  
  // Read pixels
  for (int i = 0; i < grid_size; i++) {
    if (max_value <= 0xFF) {
      unsigned char value;
      fread(&value, sizeof(unsigned char), 1, fp);
      grid_values[i] = value;
    }
    else if (max_value <= 0xFFFF) {
      unsigned short value;
      fread(&value, sizeof(unsigned short), 1, fp);
      if (little_endian) { unsigned short x = value; value = (x<<8) | (x>>8); }
      grid_values[i] = value;
    }
    else if (max_value <= 0xFFFFFFFF) {
      unsigned int value;
      fread(&value, sizeof(unsigned int), 1, fp);
      if (little_endian) { unsigned int x = value; value = (x<<24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x>>24); }
      grid_values[i] = value;
    }
    else {
      fprintf(stderr, "Unrecognized maximum value in %s\n", filename);
      return 0;
    }
  }

  // Close image
  fclose(fp);

  // Return success
  return grid_size;
}



////////////////////////////////////////////////////////////////////////
// PFM FORMAT READ/WRITE
////////////////////////////////////////////////////////////////////////

int R2Grid::
ReadPFMFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image: %s\n", filename);
    return 0;
  }

  // Read magic header
  int c;
  c = fgetc(fp); if (c != 'P') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }
  c = fgetc(fp); if (c != 'f') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }
  c = fgetc(fp); if (c != '\n') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }

  // Read width
  int width_count = 0;
  char width_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (width_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { width_string[width_count] = '\0'; break; }
    else if (!isdigit(c)) { fprintf(stderr, "Bad width character %c in %s\n", c, filename); return 0; }
    else width_string[width_count++] = c;
  }

  // Check width
  if ((width_count == 0) || (width_count > 128)) {
    fprintf(stderr, "Error reading width in %s\n", filename); 
    return 0; 
  }

  // Read height
  int height_count = 0;
  char height_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (height_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { height_string[height_count] = '\0'; break; }
    else if (!isdigit(c)) { fprintf(stderr, "Bad height character %c in %s\n", c, filename); return 0; }
    else height_string[height_count++] = c;
  }

  // Check height
  if ((height_count == 0) || (height_count > 128)) {
    fprintf(stderr, "Error reading height in %s\n", filename); 
    return 0; 
  }

  // Read endian
  int endian_count = 0;
  char endian_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (endian_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { endian_string[endian_count] = '\0'; break; }
    if (!isdigit(c) && (c != '.') && (c != '-')) { fprintf(stderr, "Bad endian character %c in %s\n", c, filename); return 0; }
    endian_string[endian_count++] = c;
  }

  // Check endian
  if ((endian_count == 0) || (endian_count > 128)) {
    fprintf(stderr, "Error reading endian in %s\n", filename); 
    return 0; 
  }

  // Parse values
  int width = atoi(width_string);
  int height = atoi(height_string);
  float endian = (float) atof(endian_string);
  if (endian == -999.0F) fprintf(stderr, "Just trying to avoid compiler warning for unused variable\n");

  // Allocate pixels
  int npixels = width * height;
  float *pixels = new float [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate pixels for %s\n", filename);
    return 0;
  }

  // Read pixels
  int count = npixels;
  while (count > 0) {
    int n = fread(&pixels[npixels-count], sizeof(float), count, fp);
    if (n <= 0) { fprintf(stderr, "Unable to read pixels from %s\n", filename); return 0; }
    count -= n;
  }

  // Fill in grid info
  grid_resolution[0] = width;
  grid_resolution[1] = height;
  grid_row_size = width;
  grid_size = npixels;
  world_to_grid_scale_factor = 1;
  world_to_grid_transform = R2identity_affine;
  grid_to_world_transform = R2identity_affine;
  if (grid_values) delete [] grid_values;
  grid_values = new RNScalar [ grid_size ];
  for (int i = 0; i < grid_size; i++) {
    if (RNIsEqual(pixels[i], R2_GRID_UNKNOWN_VALUE)) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
    else grid_values[i] = pixels[i];
  }

  // Delete pixels
  delete [] pixels;

  // Close image
  fclose(fp);

  // Return success
  return grid_size;
}



int R2Grid::
WritePFMFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open pfm image file %s", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "Pf\n");
  fprintf(fp, "%d %d\n", grid_resolution[0], grid_resolution[1]);
  fprintf(fp, "-1.0\n");

  // Allocate pixels
  float *pixels = new float [ grid_resolution[0] ];

  // Write pixels (row by row to avoid large buffers)
  for (int i = 0; i < grid_resolution[1]; i++) {
    RNScalar *valuesp = &grid_values[i * grid_resolution[0]];
    for (int i = 0; i < grid_resolution[0]; i++) pixels[i] = valuesp[i];
    if (fwrite(pixels, sizeof(float), grid_resolution[0], fp) != (unsigned int) grid_resolution[0]) {
      fprintf(stderr, "Unable to write grid values to file %s\n", filename);
      return 0;
    }
  }

  // Delete pixels
  delete [] pixels;

  // Close file
  fclose(fp);

  // Return success
  return grid_size;
}




////////////////////////////////////////////////////////////////////////
// GRD FORMAT READ/WRITE
////////////////////////////////////////////////////////////////////////

int R2Grid::
ReadGridFile(const char *filename)
{
  // Open file
  FILE *fp = stdin;
  if (filename) {
    fp = fopen(filename, "rb");
    if (!fp) {
      RNFail("Unable to open file %s", filename);
      return 0;
    }
  }

  // Read file
  int status = ReadGrid(fp);
  if (!status) return 0;

  // Close file
  fclose(fp);

  // Return number of grid values read
  return status;
}



int R2Grid::
WriteGridFile(const char *filename) const
{
  // Open file
  FILE *fp = stdout;
  if (filename) {
    fp = fopen(filename, "wb");
    if (!fp) {
      RNFail("Unable to open file %s", filename);
      return 0;
    }
  }

  // Write file
  int status = WriteGrid(fp);
  if (!status) return 0;

  // Close file
  fclose(fp);

  // Return number of grid values written
  return status;
}



int R2Grid::
ReadGrid(FILE *fp)
{
  // Read grid resolution from file
  if (fread(&grid_resolution, sizeof(int), 2, fp) != 2) {
    RNFail("Unable to read resolution from grid file");
    return 0;
  }

  // Update grid resolution variables
  grid_row_size = grid_resolution[0];
  grid_size = grid_row_size * grid_resolution[1];
  if (grid_size <= 0) {
    RNFail("Invalid grid size (%d) in file", grid_size);
    return 0;
  }

  // Read world_to_grid transformation from file
  RNScalar m[9];
  if (fread(m, sizeof(RNScalar), 9, fp) != 9) {
    RNFail("Unable to read transformation matrix from file");
    return 0;
  }

  // Update transformation variables
  world_to_grid_transform.Reset(R3Matrix(m));
  world_to_grid_scale_factor = world_to_grid_transform.ScaleFactor();
  grid_to_world_transform = world_to_grid_transform.Inverse();

  // Allocate grid values
  grid_values = new RNScalar [ grid_size ];
  assert(grid_values);

  // Read values
  int count = grid_size;
  while (count > 0) {
    int n = fread(&grid_values[grid_size-count], sizeof(RNScalar), count, fp);
    if (n < 0) { fprintf(stderr, "Unable to read pixels from grid\n"); return 0; }
    count -= n;
  }

  // Just to be sure
  for (int i = 0; i < grid_size; i++) {
    if (RNIsEqual(grid_values[i], R2_GRID_UNKNOWN_VALUE)) grid_values[i] = R2_GRID_UNKNOWN_VALUE;
  }

  // Return number of grid values read
  return grid_size;
}



int R2Grid::
WriteGrid(FILE *fp) const
{
  // Write grid resolution from file
  if (fwrite(&grid_resolution, sizeof(int), 2, fp) != 2) {
    RNFail("Unable to write resolution to file");
    return 0;
  }

  // Write world_to_grid transformation to file
  const RNScalar *m = &(world_to_grid_transform.Matrix()[0][0]);
  if (fwrite(m, sizeof(RNScalar), 9, fp) != 9) {
    RNFail("Unable to write transformation matrix to file");
    return 0;
  }

  // Write grid values (row by row to avoid large buffers)
  for (int i = 0; i < grid_resolution[1]; i++) {
    RNScalar *row_values = &grid_values[i * grid_resolution[0]];
    if (fwrite(row_values, sizeof(RNScalar), grid_resolution[0], fp) != (unsigned int) grid_resolution[0]) {
      RNFail("Unable to write grid values to grid\n");
      return 0;
    }
  }

  // Return number of grid values written
  return grid_size;
}



////////////////////////////////////////////////////////////////////////
// PNG READ/WRITE
////////////////////////////////////////////////////////////////////////

#define RN_USE_PNG
#ifdef RN_NO_PNG
#undef RN_USE_PNG
#endif



#ifdef RN_USE_PNG
# include "png/png.h"
#endif



int R2Grid::
ReadPNGFile(const char *filename)
{
#ifdef RN_USE_PNG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open PNG file %s\n", filename);
    return 0;
   }

  // Create and initialize the png_struct 
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fclose(fp);
    return 0;
  }

  // Allocate/initialize the memory for image information. 
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    return 0;
  }

  // Set up the input control if you are using standard C streams 
  png_init_io(png_ptr, fp);

  // Read the png info 
  png_read_info(png_ptr, info_ptr);

  // Extract image info 
  png_byte color_type = png_get_color_type(png_ptr, info_ptr);
  int width = png_get_image_width(png_ptr, info_ptr);
  int height = png_get_image_height(png_ptr, info_ptr);
  int rowsize = png_get_rowbytes(png_ptr, info_ptr);
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  png_byte depth = png_get_bit_depth(png_ptr, info_ptr);
  int bytes_per_pixel = depth / 8;
  assert(rowsize >= bytes_per_pixel * width);

  // Set ncomponents
  int ncomponents = 0;
  if (color_type == PNG_COLOR_TYPE_GRAY) ncomponents = 1;
  else if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ncomponents = 2;
  else if (color_type == PNG_COLOR_TYPE_RGB) ncomponents = 3;
  else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) ncomponents = 4;
  else { 
    fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    return 0;
  }

  // Allocate the pixels and row pointers 
  unsigned char *pixels = new unsigned char [ height * rowsize ]; 
  png_bytep *row_pointers = (png_bytep *) png_malloc(png_ptr, height * png_sizeof(png_bytep));
  for (int i = 0; i < height; i++) row_pointers[i] = &pixels[ (height - i - 1) * rowsize ];

  // Read the pixels 
  png_read_image(png_ptr, row_pointers);

  // Finish reading 
  png_read_end(png_ptr, info_ptr);

  // Fill in grid info
  grid_resolution[0] = width;
  grid_resolution[1] = height;
  grid_row_size = width;
  grid_size = width * height;
  world_to_grid_scale_factor = 1;
  world_to_grid_transform = R2identity_affine;
  grid_to_world_transform = R2identity_affine;
  if (grid_values) delete [] grid_values;
  grid_values = new RNScalar [ grid_size ];
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      int rgba[4];
      for (int k = 0; k < ncomponents; k++) {
        rgba[k] = 0;
        for (int d = 0; d < bytes_per_pixel; d++) {
          rgba[k] = rgba[k] << 8;
          rgba[k] += row_pointers[height - 1 - j][bytes_per_pixel*i+d];
        }
      }

      // Compute value
      RNScalar value = 0;
      if (ncomponents == 1) value = rgba[0];
      else if (ncomponents == 2) value = rgba[0];
      else if (ncomponents == 3) value = 0.3*rgba[0] + 0.59*rgba[1] + 0.11*rgba[2];
      else if (ncomponents == 4) value = 0.3*rgba[0] + 0.59*rgba[1] + 0.11*rgba[2];
      SetGridValue(i, j, value);
    }
  }

  // Free the row pointers 
  png_free(png_ptr, row_pointers);

  // Free the pixels
  delete [] pixels;

  // Clean up after the read, and free any memory allocated  
  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

  // Close the file 
  fclose(fp);

  // Return success 
  return 1;
#else
  RNFail("PNG not supported");
  return 0;
#endif
}



int R2Grid::
WritePNGFile(const char *filename) const
{
#ifdef RN_USE_PNG
  // Open the file 
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "Unable to open PNG file %s\n", filename);
    return 0;
  }
  
  // Create and initialize the png_struct 
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fclose(fp);
    return 0;
  }
  
  // Allocate/initialize the image information data. 
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fclose(fp);
    png_destroy_write_struct(&png_ptr,  NULL);
    return 0;
  }
  
  // Fill in the image data
  int width = grid_resolution[0];
  int height = grid_resolution[1];
  png_set_IHDR(png_ptr, info_ptr, width, height,
    16, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Copy the data into RNUInt16
  png_bytep data = new png_byte [ 2*grid_size ];
  png_bytep datap = data;
  for (int i = 0; i < grid_size; i++) {
    RNScalar value = grid_values[i];
    if (value == R2_GRID_UNKNOWN_VALUE) {
      *datap++ = 0;
      *datap++ = 0;
    }
    else {
      if (value > 65535) value = 65535;
      unsigned int ivalue = (unsigned int) value;
      *datap++ = (ivalue >> 8) & 0xFF;
      *datap++ = ivalue & 0xFF;
    }
  }

  // Allocate the row pointers
  png_bytep *row_pointers = (png_bytep *) png_malloc(png_ptr, height * png_sizeof(png_bytep));
  for (int i = 0; i < height; i++) row_pointers[i] = (png_bytep) &data[(height - i - 1) * 2*width];
  
  // Set up the output control 
  png_init_io(png_ptr, fp);
  
  // Write the png info 
  png_write_info(png_ptr, info_ptr);
  
  // Write the pixels 
  png_write_image(png_ptr, row_pointers);
  
  // Finish writing 
  png_write_end(png_ptr, info_ptr);
  
  // Free the row pointers 
  png_free(png_ptr, row_pointers);

  // Clean up after the write, and free any memory allocated 
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // Delete copy of data
  delete [] data;

  // Close the file 
  fclose(fp);

  // Return success
  return 1;
#else
  RNFail("PNG not supported");
  return 0;
#endif
}



////////////////////////////////////////////////////////////////////////
// IMAGE READ/WRITE
////////////////////////////////////////////////////////////////////////

int R2Grid::
ReadImage(const char *filename)
{
  // Allocate image
  R2Image *image = new R2Image();
  if (!image) {
    fprintf(stderr, "Unable to allocate image for %s\n", filename);
    return 0;
  }

  // Read image
  if (!image->Read(filename)) { 
    // fprintf(stderr, "Unable to read image from %s\n", filename);
    delete image; 
    return 0; 
  }

  // Fill in grid info
  grid_resolution[0] = image->Width();
  grid_resolution[1] = image->Height();
  grid_row_size = image->Width();
  grid_size = image->Width() * image->Height();
  world_to_grid_scale_factor = 1;
  world_to_grid_transform = R2identity_affine;
  grid_to_world_transform = R2identity_affine;
  if (grid_values) delete [] grid_values;
  grid_values = new RNScalar [ grid_size ];
  for (int j = 0; j < image->Height(); j++) {
    for (int i = 0; i < image->Width(); i++) {
      SetGridValue(i, j, image->PixelRGB(i, j).Luminance());
    }
  }

  // Delete image
  delete image;

  // Return success
  return 1;
}



int R2Grid::
WriteImage(const char *filename) const
{
  // Allocate image
  R2Image *image = new R2Image(grid_resolution[0], grid_resolution[1], 3);
  if (!image) {
    fprintf(stderr, "Unable to allocate image for %s\n", filename);
    return 0;
  }

  // Fill image pixel values
  RNInterval range = Range();
  if (range.Diameter() > 0) {
    for (int j = 0; j < image->Height(); j++) {
      for (int i = 0; i < image->Width(); i++) {
        RNScalar value = (GridValue(i, j) - range.Min()) / (range.Max() - range.Min());
        image->SetPixelRGB(i, j, RNRgb(value, value, value));
      }
    }
  }

  // Mark unknown values
  RNRgb unknown_rgb(1, 0.5, 0);
  for (int j = 0; j < image->Height(); j++) {
    for (int i = 0; i < image->Width(); i++) {
      if (GridValue(i, j) != R2_GRID_UNKNOWN_VALUE) continue;
      image->SetPixelRGB(i, j, unknown_rgb);
    }
  }

  // Write image
  if (!image->Write(filename)) { 
    // fprintf(stderr, "Unable to write image to %s\n", filename);
    delete image; 
    return 0; 
  }

  // Delete image
  delete image;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// OTHER I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R2Grid::
Print(FILE *fp) const
{
  // Check file
  if (!fp) fp = stdout;

  // Print values
  for (int j = 0; j < YResolution(); j++) {
    for (int i = 0; i < XResolution(); i++) {
      fprintf(fp, "%g ", GridValue(i, j));
    }
    fprintf(fp, "\n");
  }

  // Return number of grid values written
  return grid_size;
}




