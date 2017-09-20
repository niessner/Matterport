////////////////////////////////////////////////////////////////////////
// Source file for RGBD utility functions
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RGBD.h"



////////////////////////////////////////////////////////////////////////
// High-Level Channel Creation Functions (assumes relevant channels have been read)
////////////////////////////////////////////////////////////////////////

int RGBDCreatePositionChannels(RGBDImage *image,
  R2Grid& output_px_image, R2Grid& output_py_image, R2Grid& output_pz_image,
  RNBoolean world_coordinates)
{
  R4Matrix camera_to_world = R4identity_matrix;
  if (world_coordinates) camera_to_world = image->CameraToWorld().Matrix();
  return RGBDCreatePositionChannels(*(image->DepthChannel()),
    output_px_image, output_py_image, output_pz_image,
    image->Intrinsics(), camera_to_world);
}



int RGBDCreateNormalChannels(RGBDImage *image,
  R2Grid& output_nx_image, R2Grid& output_ny_image, R2Grid& output_nz_image,
  RNLength neighborhood_world_radius, int neighborhood_pixel_radius, RNBoolean neighborhood_search,
  RNBoolean world_coordinates)
{
  R2Grid boundary_image, px_image, py_image, pz_image, radius_image;
  if (!RGBDCreateBoundaryChannel(image, boundary_image)) return 0;
  if (!RGBDCreatePositionChannels(image, px_image, py_image, pz_image, world_coordinates)) return 0;
  return RGBDCreateNormalChannels(*(image->DepthChannel()),
    px_image, py_image, pz_image, boundary_image,
    output_nx_image, output_ny_image, output_nz_image, radius_image,
    image->WorldViewpoint(), image->WorldTowards(), image->WorldUp(),
    neighborhood_world_radius, neighborhood_pixel_radius, neighborhood_search);
}



int RGBDCreateBoundaryChannel(RGBDImage *image,
  R2Grid& output_boundary_image,
  RNScalar depth_threshold)
{
  return RGBDCreateBoundaryChannel(*(image->DepthChannel()),
    output_boundary_image, depth_threshold);
}




////////////////////////////////////////////////////////////////////////
// Mid-Level Channel Creation Functions
////////////////////////////////////////////////////////////////////////

int RGBDCreatePositionChannels(
  const R2Grid& input_undistorted_depth_image,
  R2Grid& output_px_image, R2Grid& output_py_image, R2Grid& output_pz_image,
  const R3Matrix& intrinsics_matrix, const R4Matrix& camera_to_world_matrix)
{
  // Initialize position images
  output_px_image = input_undistorted_depth_image; 
  output_px_image.Clear(R2_GRID_UNKNOWN_VALUE);
  output_py_image = output_px_image;
  output_pz_image = output_px_image;

  // Fill position images
  for (int i = 0; i < input_undistorted_depth_image.XResolution(); i++) {
    for (int j = 0; j < input_undistorted_depth_image.YResolution(); j++) {
      // Get depth
      RNScalar depth = input_undistorted_depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;

      // Get position in camera coordinate system 
      RNScalar x = ((i+0.5) - intrinsics_matrix[0][2]) * depth / intrinsics_matrix[0][0];
      RNScalar y = ((j+0.5) - intrinsics_matrix[1][2]) * depth / intrinsics_matrix[1][1];
      RNScalar z = -depth;
      R3Point position(x, y, z);

      // Transform by extrinsics matrix
      position = camera_to_world_matrix * position;

      // Fill position images
      output_px_image.SetGridValue(i, j, position.X());
      output_py_image.SetGridValue(i, j, position.Y());
      output_pz_image.SetGridValue(i, j, position.Z());
    }
  }

  // Return success 
  return 1;
}



int RGBDCreateBoundaryChannel(
  const R2Grid& input_depth_image, R2Grid& output_boundary_image,
  RNScalar depth_threshold)
{
  // Initialize position images
  output_boundary_image = input_depth_image; 
  output_boundary_image.Clear(R2_GRID_UNKNOWN_VALUE);

  // Mark bottom and top border boundaries
  for (int i = 0; i < output_boundary_image.XResolution(); i++) {
    for (int j = 0; j < output_boundary_image.YResolution(); j++) {
      output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
      if (RNIsPositive(input_depth_image.GridValue(i, j))) break;
    }
    for (int j = output_boundary_image.YResolution()-1; j >= 0; j--) {
      output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
      if (RNIsPositive(input_depth_image.GridValue(i, j))) break;
    }
  }

  // Mark left and right border boundaries
  for (int j = 0; j < output_boundary_image.YResolution(); j++) {
    for (int i = 0; i < output_boundary_image.XResolution(); i++) {
      output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
      if (RNIsPositive(input_depth_image.GridValue(i, j))) break;
    }
    for (int i = output_boundary_image.XResolution()-1; i >= 0; i--) {
      output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
      if (RNIsPositive(input_depth_image.GridValue(i, j))) break;
    }
  }

  // Create copy of depth image with holes filled
  R2Grid filled_input_depth_image(input_depth_image);
  filled_input_depth_image.Substitute(0, R2_GRID_UNKNOWN_VALUE);
  filled_input_depth_image.FillHoles();

  // Mark interior boundaries, silhouettes, and shadows
  for (int i = 1; i < output_boundary_image.XResolution()-1; i++) {
    for (int j = 1; j < output_boundary_image.YResolution()-1; j++) {
      // Get original depth
      RNScalar depth = input_depth_image.GridValue(i, j);

      // Check if in hole
      if (RNIsNegativeOrZero(depth)) {
        output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
        continue;
      }

      // Get filled depth
      depth = filled_input_depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;

      // Check depth relative to horizontal neighbors
      for (int k = 0; k < 4; k++) {
        int s = (k < 3) ? -1 : 0;
        int t = (k < 3) ? k-1 : -1;

        // Get depth on one side
        RNScalar depthA = filled_input_depth_image.GridValue(i-s, j-t);
        if (RNIsNegativeOrZero(depthA)) {
          output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
          break;
        }

        // Get depth on other side
        RNScalar depthB = filled_input_depth_image.GridValue(i+s, j+t);
        if (RNIsNegativeOrZero(depthB)) {
          output_boundary_image.SetGridValue(i, j, RGBD_BORDER_BOUNDARY);
          break;
        }

        // Check differences of depth for shadow/silhouette
        RNScalar deltaA = depth - depthA;
        RNScalar deltaB = depthB - depth;
        RNScalar threshold = depth * depth_threshold;
        if (threshold < 0.1) threshold = 0.1;
        if (deltaA < -threshold) {
          if (deltaA < 4*deltaB) {
            output_boundary_image.SetGridValue(i-s, j-t, RGBD_SHADOW_BOUNDARY);
            output_boundary_image.SetGridValue(i, j, RGBD_SILHOUETTE_BOUNDARY);
          }
        }
        else if (deltaA > threshold) {
          if (deltaA > 4*deltaB) {
            output_boundary_image.SetGridValue(i-s, j-t, RGBD_SILHOUETTE_BOUNDARY);
            output_boundary_image.SetGridValue(i, j, RGBD_SHADOW_BOUNDARY);
          }
        }
        if (deltaB < -threshold) {
          if (deltaB < 4*deltaA) {
            output_boundary_image.SetGridValue(i+s, j+t, RGBD_SILHOUETTE_BOUNDARY);
            output_boundary_image.SetGridValue(i, j, RGBD_SHADOW_BOUNDARY);
          }
        }
        else if (deltaB > threshold) {
          if (deltaB > 4*deltaA) {
            output_boundary_image.SetGridValue(i+s, j+t, RGBD_SHADOW_BOUNDARY);
            output_boundary_image.SetGridValue(i, j, RGBD_SILHOUETTE_BOUNDARY);
          }
        }
      }
    }
  }

  // Return success 
  return 1;
}



int RGBDCreateNormalChannels(const R2Grid& input_depth_image, 
  const R2Grid& input_px_image, const R2Grid& input_py_image, const R2Grid& input_pz_image, const R2Grid& boundary_image,
  R2Grid& output_nx_image, R2Grid& output_ny_image, R2Grid& output_nz_image, R2Grid& output_radius_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  RNLength neighborhood_world_radius, int neighborhood_pixel_radius, RNBoolean neighborhood_search)
{
  // Initialize normal images
  output_nx_image = input_depth_image; 
  output_nx_image.Clear(R2_GRID_UNKNOWN_VALUE);
  output_ny_image = output_nx_image;
  output_nz_image = output_nx_image;
  output_radius_image = output_nx_image;

  // Allocate array of neighborhood points
  int pixel_radius = neighborhood_pixel_radius * input_depth_image.XResolution() / 640.0 + 0.5;
  if (pixel_radius == 0) pixel_radius = 1;
  int neighborhood_pixel_radius_squared = pixel_radius * pixel_radius;
  RNScalar neighborhood_world_radius_squared = neighborhood_world_radius * neighborhood_world_radius;
  R3Point *neighborhood_points = new R3Point [ input_depth_image.NEntries() ];

  // Allocate search grid
  R2Grid search_grid(input_depth_image.XResolution(), input_depth_image.YResolution());
  const RNScalar *search_valuesp = search_grid.GridValues();
  search_grid.Clear(-1);

  // Fill normal images
  for (int i = 0; i < input_depth_image.XResolution(); i++) {
    for (int j = 0; j < input_depth_image.YResolution(); j++) {
      int neighborhood_npoints = 0;

      // Check depth
      RNScalar depth = input_depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;

      // Get image position
      R2Point image_position(i, j);

      // Get world position
      RNScalar x = input_px_image.GridValue(i, j);
      RNScalar y = input_py_image.GridValue(i, j);
      RNScalar z = input_pz_image.GridValue(i, j);
      R3Point world_position(x, y, z);

      // Create array of neighborhood points
      if (neighborhood_search) {
        // Initialize stack
        int seed_index;
        RNArray<const RNScalar *> stack;
        search_grid.IndicesToIndex(i, j, seed_index);
        const RNScalar *seed_valuep = &search_valuesp[seed_index];
        search_grid.SetGridValue(seed_index, seed_index);
        stack.Insert(seed_valuep);

        // Depth first search, not extending beyond boundaries
        while (!stack.IsEmpty()) {
          // Pop current point off stack
          int ix, iy;
          const RNScalar *current_valuep = stack.Tail(); stack.RemoveTail();
          int current_index = current_valuep - search_valuesp;
          search_grid.IndexToIndices(current_index, ix, iy);

          // Get current info
          RNScalar current_x = input_px_image.GridValue(ix, iy);
          RNScalar current_y = input_py_image.GridValue(ix, iy);
          RNScalar current_z = input_pz_image.GridValue(ix, iy);
          R3Point current_position(current_x, current_y, current_z);
          RNScalar current_boundary_value = boundary_image.GridValue(ix, iy);
          int current_boundary = (int) (current_boundary_value + 0.5);
        
          // Add point to array
          neighborhood_points[neighborhood_npoints++] = current_position;

          // Add adjacent points to stack
          for (int inx = ix-1; inx <= ix+1; inx++) {
            if ((inx < 0) || (inx >= search_grid.XResolution())) continue;
            for (int iny = iy-1; iny <= iy+1; iny++) {
              if ((iny < 0) || (iny >= search_grid.YResolution())) continue;

              // Check if neighbor has already been visited
              RNScalar neighbor_value = search_grid.GridValue(inx, iny);
              if (RNIsEqual(neighbor_value, seed_index)) continue;
              search_grid.SetGridValue(inx, iny, seed_index);

              // Check neighbor depth
              RNScalar neighbor_depth = input_depth_image.GridValue(inx, iny);
              if (RNIsNegativeOrZero(neighbor_depth)) continue;

              // Check neighbor image distance
              R2Point neighbor_image_position(inx, iny);
              if (R2SquaredDistance(image_position, neighbor_image_position) > neighborhood_pixel_radius_squared) continue;
            
              // Check neighbor world distance
              RNScalar neighbor_x = input_px_image.GridValue(inx, iny);
              RNScalar neighbor_y = input_py_image.GridValue(inx, iny);
              RNScalar neighbor_z = input_pz_image.GridValue(inx, iny);
              R3Point neighbor_world_position(neighbor_x, neighbor_y, neighbor_z);
              if (neighborhood_world_radius_squared > 0) {
                if (R3SquaredDistance(world_position, neighbor_world_position) > neighborhood_world_radius_squared) continue;
              }

              // Check if neighbor is across boundary
              int neighbor_boundary = (int) (boundary_image.GridValue(inx, iny) + 0.5);
              if ((current_boundary == RGBD_SHADOW_BOUNDARY) && (neighbor_boundary == RGBD_SILHOUETTE_BOUNDARY)) continue;
              if ((current_boundary == RGBD_SILHOUETTE_BOUNDARY) && (neighbor_boundary == RGBD_SHADOW_BOUNDARY)) continue;
        
              // Add neighbor to search
              int neighbor_index;
              search_grid.IndicesToIndex(inx, iny, neighbor_index);
              const RNScalar *neighbor_valuep = &search_valuesp[neighbor_index];
              stack.Insert(neighbor_valuep);
            }
          }
        }
      }
      else {
        // Add neighbor points to array
        for (int inx = i-neighborhood_pixel_radius; inx <= i+neighborhood_pixel_radius; inx++) {
          if ((inx < 0) || (inx >= input_depth_image.XResolution())) continue;
          for (int iny = j-neighborhood_pixel_radius; iny <= j+neighborhood_pixel_radius; iny++) {
            if ((iny < 0) || (iny >= input_depth_image.YResolution())) continue;
            
            // Check neighbor depth
            RNScalar neighbor_depth = input_depth_image.GridValue(inx, iny);
            if (RNIsNegativeOrZero(neighbor_depth)) continue;

            // Check neighbor world position
            RNScalar neighbor_x = input_px_image.GridValue(inx, iny);
            RNScalar neighbor_y = input_py_image.GridValue(inx, iny);
            RNScalar neighbor_z = input_pz_image.GridValue(inx, iny);
            R3Point neighbor_world_position(neighbor_x, neighbor_y, neighbor_z);
            if (neighborhood_world_radius_squared > 0) {
              if (R3SquaredDistance(world_position, neighbor_world_position) > neighborhood_world_radius_squared) continue;
            }

            // Add neighbor point to array
            neighborhood_points[neighborhood_npoints++] = neighbor_world_position;
          }
        }
      }

      // Check number of neighbor points
      if (neighborhood_npoints < 3) continue;

      // Solve for normal
      RNScalar variances[3];
      R3Point centroid = R3Centroid(neighborhood_npoints, neighborhood_points);
      R3Triad triad = R3PrincipleAxes(centroid, neighborhood_npoints, neighborhood_points, NULL, variances);
      R3Vector normal = triad[2];

      // Flip normal to point towards camera
      R3Plane plane(world_position, normal);
      if (R3SignedDistance(plane, viewpoint) < 0) normal.Flip(); 

      // Fill normal images
      output_nx_image.SetGridValue(i, j, normal.X());
      output_ny_image.SetGridValue(i, j, normal.Y());
      output_nz_image.SetGridValue(i, j, normal.Z());

      // Fill radius image
      RNScalar r = sqrt(variances[0]);
      if (neighborhood_pixel_radius > 1) r /= neighborhood_pixel_radius;
      output_radius_image.SetGridValue(i, j, r);
    }
  }

  // Delete array of neighborhood points
  delete [] neighborhood_points;

  // Return success 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// LOW-LEVEL UNDISTORTION UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////////////

int RGBDCreateUndistortedColorImage(
  const R2Grid& input_distorted_color_image,
  R2Grid& output_undistorted_color_image,
  const RGBDCamera& distorted_camera, const RGBDCamera& undistorted_camera)
{
  // Check output camera
  assert((undistorted_camera.k1_color == 0) && (undistorted_camera.k2_color == 0) &&
         (undistorted_camera.k3_color == 0) && (undistorted_camera.k4_color == 0) &&
         (undistorted_camera.p1_color == 0) && (undistorted_camera.p2_color == 0));

  // Get convenient variables
  if ((distorted_camera.colorWidth == 0) || (distorted_camera.colorHeight == 0)) return 0;
  if ((undistorted_camera.colorWidth == 0) || (undistorted_camera.colorHeight == 0)) return 0;
  double sx = (double) distorted_camera.colorWidth / (double) undistorted_camera.colorWidth;
  double sy = (double) distorted_camera.colorHeight / (double) undistorted_camera.colorHeight;

  // Check if don't need to scale or undistort
  if ((RNIsEqual(sx, 1.0) && RNIsEqual(sy, 1.0)) &&
      (distorted_camera.mx_color == undistorted_camera.mx_color) && (distorted_camera.my_color == undistorted_camera.my_color) &&
      (distorted_camera.k1_color == undistorted_camera.k1_color) && (distorted_camera.k2_color == undistorted_camera.k2_color) &&
      (distorted_camera.k3_color == undistorted_camera.k3_color) && (distorted_camera.k4_color == undistorted_camera.k4_color) &&
      (distorted_camera.p1_color == undistorted_camera.p1_color) && (distorted_camera.p2_color == undistorted_camera.p2_color)) {
    output_undistorted_color_image = input_distorted_color_image;
    return 1;
  }

  // Save copy of input image (so can pass same image both in and out)
  R2Grid source_color_image(input_distorted_color_image);

  // Initialize the output image
  output_undistorted_color_image = R2Grid(undistorted_camera.colorWidth, undistorted_camera.colorHeight);
  output_undistorted_color_image.Clear(R2_GRID_UNKNOWN_VALUE);
  
  // Compute output image via reverse mapping
  for (int output_iy = 0; output_iy < output_undistorted_color_image.YResolution(); output_iy++) {
    for (int output_ix = 0; output_ix < output_undistorted_color_image.XResolution(); output_ix++) {
      // Compute coordinates in undistorted output image
      R2Point undistorted_output_position(output_ix + 0.5, output_iy + 0.5);

      // Compute coordinates in undistorted input image
      // Note that this ignores p1 and p2
      R2Point undistorted_input_position(undistorted_output_position);
      RNScalar x = sx * (undistorted_output_position.X() - undistorted_camera.mx_color) + distorted_camera.mx_color;
      RNScalar y = sy * (undistorted_output_position.Y() - undistorted_camera.my_color) + distorted_camera.my_color;
      undistorted_input_position.Reset(x, y);
      
      // Compute coordinates in distorted input image
      R2Point distorted_input_position(undistorted_input_position);
      if ((distorted_camera.k1_color != 0.0) || (distorted_camera.k2_color != 0.0) ||
          (distorted_camera.k3_color != 0.0) || (distorted_camera.k4_color != 0.0) ||
          (distorted_camera.p1_color != 0.0) || (distorted_camera.k2_color != 0.0)) {
        double nx = (undistorted_input_position.X() - distorted_camera.mx_color) / distorted_camera.fx_color;
        double ny = (undistorted_input_position.Y() - distorted_camera.my_color) / distorted_camera.fy_color;
        double rr = nx*nx + ny*ny; double rrrr = rr*rr;
        double s = 1.0 + rr*distorted_camera.k1_color + rrrr*distorted_camera.k2_color + rrrr*rr*distorted_camera.k3_color + rrrr*rrrr*distorted_camera.k4_color;
        nx = s*nx + distorted_camera.p2_color*(rr + 2*nx*nx) + 2*distorted_camera.p1_color*nx*ny;
        ny = s*ny + distorted_camera.p1_color*(rr + 2*ny*ny) + 2*distorted_camera.p2_color*nx*ny;
        double x = nx*distorted_camera.fx_color + distorted_camera.mx_color;
        double y = ny*distorted_camera.fy_color + distorted_camera.my_color;
        distorted_input_position.Reset(x, y);
      }

      // Check if outside input image
      if (distorted_input_position.X() < 0) continue;
      if (distorted_input_position.X() > source_color_image.XResolution()-1) continue;
      if (distorted_input_position.Y() < 0) continue;
      if (distorted_input_position.Y() > source_color_image.YResolution()-1) continue;

#if 0
      // Check if center goes to center
      if (RNIsEqual(undistorted_output_position.X(), distorted_camera.mx_color, 0.5) &&
          RNIsEqual(undistorted_output_position.Y(), distorted_camera.my_color, 0.5) &&
          !R2Contains(undistorted_output_position, distorted_input_position)) RNAbort("HERE");
#endif
      
      // Copy value from input image to output image (bilinear interpolation)
      RNScalar value = source_color_image.GridValue(distorted_input_position);
      output_undistorted_color_image.SetGridValue(output_ix, output_iy, value);
    }
  }

  // Return success
  return 1;
}



int RGBDCreateUndistortedDepthImage(
   const R2Grid& input_distorted_depth_image,
   R2Grid& output_undistorted_depth_image,
   const RGBDCamera& distorted_camera, RGBDCamera& undistorted_camera)
{
  // Check output camera
  assert((undistorted_camera.k1_depth == 0) && (undistorted_camera.k2_depth == 0) &&
         (undistorted_camera.k3_depth == 0) && (undistorted_camera.k4_depth == 0) &&
         (undistorted_camera.p1_depth == 0) && (undistorted_camera.p2_depth == 0) &&
         (undistorted_camera.depthToColorExtrinsics.IsIdentity()));

  // Check inputs
  if ((distorted_camera.depthWidth == 0) || (distorted_camera.depthHeight == 0)) return 0;
  if ((undistorted_camera.depthWidth == 0) || (undistorted_camera.depthHeight == 0)) return 0;
  if ((input_distorted_depth_image.XResolution() == 0) || (input_distorted_depth_image.YResolution() == 0)) return 0;
  
  // Check if don't need to undistort or warp
  if ((undistorted_camera.depthWidth == distorted_camera.depthWidth) &&
      (undistorted_camera.depthHeight == distorted_camera.depthHeight) && 
      (distorted_camera.mx_depth == 0) && (distorted_camera.my_depth == 0) &&
      (distorted_camera.k1_depth == 0) && (distorted_camera.k2_depth == 0) &&
      (distorted_camera.k3_depth == 0) && (distorted_camera.k4_depth == 0) &&
      (distorted_camera.p1_depth == 0) && (distorted_camera.p2_depth == 0) &&
      (distorted_camera.depthToColorExtrinsics.IsIdentity())) {
    output_undistorted_depth_image = input_distorted_depth_image;
    return 1;
  }

  // Save copy of input image (so can pass same image both in and out)
  R2Grid source_depth_image(input_distorted_depth_image);

  // Initialize the output image
  output_undistorted_depth_image = R2Grid(undistorted_camera.depthWidth, undistorted_camera.depthHeight);
  output_undistorted_depth_image.Clear(R2_GRID_UNKNOWN_VALUE);

#ifdef MASK_HOLES
  // Create image with no holes 
  source_depth_image.Substitute(0, R2_GRID_UNKNOWN_VALUE);
  source_depth_image.FillHoles();
#endif
  
  // Compute output image via a combination of reverse and forward mapping
  double step_factor = 0.5; // make <1 to avoid pin-holes in output image (or use InterpolateDepthImage, which is faster but not as good)
  double x_step = step_factor * source_depth_image.XResolution() / (double) output_undistorted_depth_image.XResolution();
  double y_step = step_factor * source_depth_image.YResolution() / (double) output_undistorted_depth_image.YResolution();
  for (double undistorted_input_y = 0.5*y_step; undistorted_input_y < source_depth_image.YResolution(); undistorted_input_y += y_step) {
    for (double undistorted_input_x = 0.5*x_step; undistorted_input_x < source_depth_image.XResolution(); undistorted_input_x += x_step) {
      // Compute coordinates in distorted input image
      // Note that this ignores p1 and p2
      R2Point distorted_input_position(undistorted_input_x, undistorted_input_y);
      if ((distorted_camera.k1_depth != 0.0) || (distorted_camera.k2_depth != 0.0) ||
          (distorted_camera.k3_depth != 0.0) || (distorted_camera.k3_depth != 0.0) ||
          (distorted_camera.p1_depth != 0.0) || (distorted_camera.p2_depth != 0.0)) {
        double nx = (undistorted_input_x - distorted_camera.mx_depth) / distorted_camera.fx_depth;
        double ny = (undistorted_input_y - distorted_camera.my_depth) / distorted_camera.fy_depth;
        double rr = nx*nx + ny*ny; double rrrr = rr*rr;
        double s = 1.0 + rr*distorted_camera.k1_depth + rrrr*distorted_camera.k2_depth + rrrr*rr*distorted_camera.k3_depth + rrrr*rrrr*distorted_camera.k4_depth;
        nx = s*nx + distorted_camera.p2_depth*(rr + 2*nx*nx) + 2*distorted_camera.p1_depth*nx*ny;
        ny = s*ny + distorted_camera.p1_depth*(rr + 2*ny*ny) + 2*distorted_camera.p2_depth*nx*ny;
        double x = nx*distorted_camera.fx_depth + distorted_camera.mx_depth;
        double y = ny*distorted_camera.fy_depth + distorted_camera.my_depth;
        distorted_input_position.Reset(x, y);
      }

      // Get closest distorted input pixel
      int distorted_input_ix = (int) (distorted_input_position.X() + 0.5);
      int distorted_input_iy = (int) (distorted_input_position.Y() + 0.5);
      if (distorted_input_ix < 0) continue;
      if (distorted_input_ix >= source_depth_image.XResolution()) continue;
      if (distorted_input_iy < 0) continue;
      if (distorted_input_iy >= source_depth_image.YResolution()) continue;

#ifdef MASK_HOLES
      // Get depth, even if there was a hole there originally
      RNBoolean depth_was_missing = FALSE;
      RNScalar depth = source_depth_image.GridValue(distorted_input_ix, distorted_input_iy);
      if (depth <= 0) {
        depth = interpolated_depth_image.GridValue(distorted_input_ix, distorted_input_iy);
        if (depth <= 0) continue;
        depth_was_missing = TRUE;
      }
#else
      // Get/check depth
      RNScalar depth = source_depth_image.GridValue(distorted_input_ix, distorted_input_iy);
      if (depth <= 0) continue;
#endif
      
      // Find 3D position in depth camera coordinates
      R3Point world_position;
      world_position[0] = (undistorted_input_x - distorted_camera.mx_depth) * depth / distorted_camera.fx_depth;
      world_position[1] = (undistorted_input_y - distorted_camera.my_depth) * depth / distorted_camera.fy_depth;
      world_position[2] = -depth;

      // Transform into output (color) camera coordinates
      world_position = distorted_camera.depthToColorExtrinsics * world_position;

      // Project into output pixel coordinates
      R3Point undistorted_output_position;
      undistorted_output_position[0] = undistorted_camera.mx_depth + world_position[0] * undistorted_camera.fx_depth / -world_position[2];
      undistorted_output_position[1] = undistorted_camera.my_depth + world_position[1] * undistorted_camera.fy_depth / -world_position[2];
      undistorted_output_position[2] = world_position[2];

      // Snap to closest pixel in output image
      int undistorted_output_ix = (int) (undistorted_output_position[0] + 0.5);
      if (undistorted_output_ix < 0) continue;
      if (undistorted_output_ix >= output_undistorted_depth_image.XResolution()) continue;
      int undistorted_output_iy = (int) (undistorted_output_position[1] + 0.5);
      if (undistorted_output_iy < 0) continue;
      if (undistorted_output_iy >= output_undistorted_depth_image.YResolution()) continue;

#ifdef MASK_HOLES
      // Assign depth of zero in pixels where a hole is mapped
      if (depth_was_missing) undistorted_output_position[2] = 0.0;
#endif
      
      // Store only closest grid value at every ix,iy in output image
      RNScalar new_depth = -undistorted_output_position[2];
      RNScalar old_depth = output_undistorted_depth_image.GridValue(undistorted_output_ix, undistorted_output_iy);
      if ((old_depth == 0.0) || (old_depth == R2_GRID_UNKNOWN_VALUE) || (new_depth < old_depth)) {
        output_undistorted_depth_image.SetGridValue(undistorted_output_ix, undistorted_output_iy, new_depth);
      }

#if 0
      // Check if center goes to center
      if (distorted_camera.depthToColorExtrinsics.IsIdentity() &&
          undistorted_camera.depthToColorExtrinsics.IsIdentity() &&
          RNIsEqual(undistorted_output_position.X(), undistorted_camera.mx_depth, 0.5) &&
          RNIsEqual(undistorted_output_position.Y(), undistorted_camera.my_depth, 0.5)) {
        if (RNIsNotEqual(distorted_input_position.X(), distorted_camera.mx_depth, 0.5) ||
            RNIsNotEqual(distorted_input_position.Y(), distorted_camera.my_depth, 0.5)) {
          printf("%g %g | %g %g\n",
                 distorted_input_position.X(), distorted_input_position.Y(),
                 undistorted_output_position.X(), undistorted_output_position.Y());
                 
        }
      }
#endif
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// LOW-LEVEL IMAGE RESAMPLING UTILITIES
////////////////////////////////////////////////////////////////////////

int RGBDResampleDepthImage(R2Grid& image, R3Matrix& intrinsics_matrix, int xres, int yres)
{
  // Just checking
  if ((xres == 0) || (yres == 0)) return 0;
  if ((image.XResolution() == 0) || (image.YResolution() == 0)) return 0;

  // Compute scale factors
  double xscale = (double) image.XResolution() / (double) xres;
  double yscale = (double) image.YResolution() / (double) yres;

  // Copy original image and then initialize to zeroes
  R2Grid copy_image(image);
  image = R2Grid(xres, yres);

  // Resample images
  for (int iy = 0; iy < yres; iy++) {
    double cy = yscale * iy;
    int min_y = cy - 0.5*yscale + 0.5;
    if (min_y >= copy_image.YResolution()) continue;
    int max_y = cy + 0.5*yscale + 0.5;
    if (max_y < 0) continue;
    if (min_y < 0) min_y = 0;
    if (max_y >= copy_image.YResolution()) max_y = copy_image.YResolution()-1;
    if (max_y < min_y) max_y = min_y;
    
    for (int ix = 0; ix < xres; ix++) {
      double cx = xscale * ix;
      int min_x = cx - 0.5*xscale + 0.5;
      if (min_x >= copy_image.XResolution()) continue;
      int max_x = cx + 0.5*xscale + 0.5;
      if (max_x < 0) continue;
      if (min_x < 0) min_x = 0;
      if (max_x >= copy_image.XResolution()) max_x = copy_image.XResolution()-1;
      if (max_x < min_x) max_x = min_x;

      // Find minimum depth in neighborhood
      RNScalar min_depth = FLT_MAX;
      for (int j = min_y; j <= max_y; j++) {
        for (int i = min_x; i <= max_x; i++) {
          RNScalar depth = copy_image.GridValue(i, j);
          if ((depth <= 0) || (depth == R2_GRID_UNKNOWN_VALUE)) continue;
          if (depth < min_depth) min_depth = depth;
        }
      }

      // Set depth
      if (min_depth < FLT_MAX) {
        image.SetGridValue(ix, iy, min_depth);
      }
    }
  }

  // Update depth intrinsics matrix
  intrinsics_matrix[0][0] /= xscale;
  intrinsics_matrix[0][2] /= xscale;
  intrinsics_matrix[1][1] /= yscale;
  intrinsics_matrix[1][2] /= yscale;

  // Return success
  return 1;
}



int RGBDResampleColorImage(R2Image& image, R3Matrix& intrinsics_matrix, int xres, int yres)
{
  // Just checking
  if ((xres == 0) || (yres == 0)) return 0;
  if ((image.Width() == 0) || (image.Height() == 0)) return 0;

  // Compute scale factors
  double xscale = (double) image.Width() / (double) xres;
  double yscale = (double) image.Height() / (double) yres;
  if (RNIsNotEqual(xscale, yscale, 0.01)) {
    fprintf(stderr, "Warning: anisotropic scaling of color image by factor %g\n", yscale/xscale);
  }
  
  // Copy original image and then initialize to zeroes
  R2Image copy_image(image);
  image = R2Image(xres, yres, 3);

  // Resample images
  for (int iy = 0; iy < yres; iy++) {
    double cy = yscale * iy;
    int min_y = cy - 0.5*yscale + 0.5;
    if (min_y >= copy_image.Height()) continue;
    int max_y = cy + 0.5*yscale + 0.5;
    if (max_y < 0) continue;
    if (min_y < 0) min_y = 0;
    if (max_y >= copy_image.Height()) max_y = copy_image.Height()-1;
    if (max_y < min_y) max_y = min_y;

    for (int ix = 0; ix < xres; ix++) {
      double cx = xscale * ix;
      int min_x = cx - 0.5*xscale + 0.5;
      if (min_x >= copy_image.Width()) continue;
      int max_x = cx + 0.5*xscale + 0.5;
      if (max_x < 0) continue;
      if (min_x < 0) min_x = 0;
      if (max_x >= copy_image.Width()) max_x = copy_image.Width()-1;
      if (max_x < min_x) max_x = min_x;

      // Compute weighted sum of colors in neighborhood
      RNRgb sum(0,0,0);
      RNScalar weight = 0;
      for (int j = min_y; j <= max_y; j++) {
        for (int i = min_x; i <= max_x; i++) {
          RNRgb color = copy_image.PixelRGB(i, j);
          if (color == RNblack_rgb) continue; // ??? for borders
          RNScalar w = 1.0;
          sum += w * color;
          weight += w;
        }
      }

      // Compute average
      if (weight > 0) {
        RNRgb color = sum / weight;
        image.SetPixelRGB(ix, iy, color);
      }
    }
  }

  // Update depth intrinsics matrix
  intrinsics_matrix[0][0] /= xscale;
  intrinsics_matrix[0][2] /= xscale;
  intrinsics_matrix[1][1] /= yscale;
  intrinsics_matrix[1][2] /= yscale;

  // Return success
  return 1;
}



