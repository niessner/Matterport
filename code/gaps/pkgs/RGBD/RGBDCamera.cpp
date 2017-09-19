// Source file for RGBDCamera class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RGBD/RGBD.h"



////////////////////////////////////////////////////////////////////////
// Member functions
////////////////////////////////////////////////////////////////////////

RGBDCamera::
RGBDCamera(void)
  : colorWidth(0),    
    colorHeight(0), 
    depthWidth(0),  
    depthHeight(0), 
    fx_color(0),
    fy_color(0),
    mx_color(0), 
    my_color(0), 
    k1_color(0),
    k2_color(0),
    k3_color(0),
    k4_color(0),
    p1_color(0),
    p2_color(0),
    fx_depth(0), 
    fy_depth(0), 
    mx_depth(0), 
    my_depth(0),
    k1_depth(0),
    k2_depth(0),
    k3_depth(0),
    k4_depth(0),
    p1_depth(0),
    p2_depth(0),
    deviceId(NULL),
    deviceName(NULL),        
    sceneLabel(NULL),        
    sceneType(NULL),         
    numDepthFrames(0),      
    numColorFrames(0),      
    numIMUmeasurements(0)
{
  // Initialize everything
  depthToColorExtrinsics = R4identity_matrix;
}



RGBDCamera::
RGBDCamera(const RGBDCamera& camera)
  : colorWidth(camera.colorWidth),    
    colorHeight(camera.colorHeight), 
    depthWidth(camera.depthWidth),  
    depthHeight(camera.depthHeight), 
    fx_color(camera.fx_color),
    fy_color(camera.fy_color),
    mx_color(camera.mx_color), 
    my_color(camera.my_color), 
    k1_color(camera.k1_color),
    k2_color(camera.k2_color),
    k3_color(camera.k3_color),
    k4_color(camera.k4_color),
    p1_color(camera.p1_color),
    p2_color(camera.p2_color),
    fx_depth(camera.fx_depth), 
    fy_depth(camera.fy_depth), 
    mx_depth(camera.mx_depth), 
    my_depth(camera.my_depth),
    k1_depth(camera.k1_depth),
    k2_depth(camera.k2_depth),
    k3_depth(camera.k3_depth),
    k4_depth(camera.k4_depth),
    p1_depth(camera.p1_depth),
    p2_depth(camera.p2_depth),
    deviceId((camera.deviceId) ? strdup(camera.deviceId) : NULL),
    deviceName((camera.deviceName) ? strdup(camera.deviceName) : NULL),
    sceneLabel((camera.sceneLabel) ? strdup(camera.sceneLabel) : NULL),
    sceneType((camera.sceneType) ? strdup(camera.sceneType) : NULL),
    numDepthFrames(camera.numDepthFrames),      
    numColorFrames(camera.numColorFrames),      
    numIMUmeasurements(camera.numIMUmeasurements)
{
  // Initialize everything
  depthToColorExtrinsics = R4identity_matrix;
}



RGBDCamera::
~RGBDCamera(void)
{
  // Free everything
  if (deviceId) free(deviceId);
  if (deviceName) free(deviceName);
  if (sceneLabel) free(sceneLabel);
  if (sceneType) free(sceneType);
}



void RGBDCamera::
SetUndistortedParameters(void)
{
  // Do not change color center to minimize processing of color image
  // mx_color = 0.5 * colorWidth;
  // my_color = 0.5 * colorHeight;

  // Set color parameters
  k1_color = 0;
  k2_color = 0;
  k3_color = 0;
  k4_color = 0;
  p1_color = 0;
  p2_color = 0;

  // Scale focal length
  RNScalar xscaling = (depthWidth == 0) ? (double) colorWidth / (double) depthWidth : 1.0;
  RNScalar yscaling = (depthHeight == 0) ? (double) colorHeight / (double) depthHeight : 1.0;
  if (RNIsNotEqual(xscaling, 1.0)) fx_depth *= xscaling;
  if (RNIsNotEqual(yscaling, 1.0)) fy_depth *= yscaling;

  // Set depth parameters
  depthWidth = colorWidth;
  depthHeight = colorHeight;
  mx_depth = mx_color;
  my_depth = my_color;
  k1_depth = k1_color;
  k2_depth = k2_color;
  k3_depth = k3_color;
  k4_depth = k4_color;
  p1_depth = p1_color;
  p2_depth = p2_color;
}



////////////////////////////////////////////////////////////////////////
// ScanNet i/o functions
////////////////////////////////////////////////////////////////////////

int RGBDCamera::
ReadScanNetFile(const char *filename)
{
  // Open camera file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open camera file %s\n", filename);
    return 0;
  }

  // Read camera file
  double tmp_my_color = 0;
  double tmp_my_depth = 0;
  char buffer[4096], key[1024], value[1024];
  while (fgets(buffer, 4096, fp)) {
    if (sscanf(buffer, "%s = %s", key, value) != (unsigned int) 2) continue;
    if (!strcmp(key, "color_width")) colorWidth = atoi(value);    
    else if (!strcmp(key, "colorWidth")) colorWidth = atoi(value);    
    else if (!strcmp(key, "color_height")) colorHeight= atoi(value); 
    else if (!strcmp(key, "colorHeight")) colorHeight= atoi(value); 
    else if (!strcmp(key, "depth_width")) depthWidth = atoi(value);  
    else if (!strcmp(key, "depthWidth")) depthWidth = atoi(value);  
    else if (!strcmp(key, "depth_height")) depthHeight = atoi(value); 
    else if (!strcmp(key, "depthHeight")) depthHeight = atoi(value); 
    else if (!strcmp(key, "fx_color")) fx_color = atof(value);
    else if (!strcmp(key, "fy_color")) fy_color = atof(value);
    else if (!strcmp(key, "mx_color")) mx_color = atof(value); 
    else if (!strcmp(key, "my_color")) tmp_my_color = atof(value); 
    else if (!strcmp(key, "fx_depth")) fx_depth = atof(value); 
    else if (!strcmp(key, "fy_depth")) fy_depth = atof(value); 
    else if (!strcmp(key, "mx_depth")) mx_depth = atof(value); 
    else if (!strcmp(key, "my_depth")) tmp_my_depth = atof(value); 
    else if (!strcmp(key, "k1_color")) k1_color = atof(value); 
    else if (!strcmp(key, "k2_color")) k2_color = atof(value); 
    else if (!strcmp(key, "k1_depth")) k1_depth = atof(value); 
    else if (!strcmp(key, "k2_depth")) k2_depth = atof(value); 
    else if (!strcmp(key, "deviceId")) deviceId = strdup(value);
    else if (!strcmp(key, "deviceName")) deviceName = strdup(value);        
    else if (!strcmp(key, "sceneLabel")) sceneLabel = strdup(value);        
    else if (!strcmp(key, "sceneType")) sceneType = strdup(value);         
    else if (!strcmp(key, "numDepthFrames")) numDepthFrames = atoi(value);      
    else if (!strcmp(key, "numColorFrames")) numColorFrames = atoi(value);      
    else if (!strcmp(key, "numIMUmeasurements")) numIMUmeasurements = atoi(value);
    else if (!strcmp(key, "colorToDepthExtrinsics")) {
      // Note: variable name in file is inverse of what it is !?!?!?
      double m[16];
      sscanf(buffer, "%s = %lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", key,
        &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7],
        &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]);
      R4Matrix yzflip(1, 0, 0, 0,  0, -1, 0, 0,  0, 0, -1, 0,  0, 0, 0, 1);
      depthToColorExtrinsics = yzflip * R4Matrix(m) *yzflip;
    }      
    else if (!strcmp(key, "depthToColorExtrinsics")) {
      double m[16];
      sscanf(buffer, "%s = %lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", key,
        &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7],
        &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]);
      R4Matrix yzflip(1, 0, 0, 0,  0, -1, 0, 0,  0, 0, -1, 0,  0, 0, 0, 1);
      depthToColorExtrinsics = yzflip * R4Matrix(m) *yzflip;
    }      
  }

  // Close camera file
  fclose(fp);

  // Adjust my (specified with origin at top-left)
  if (tmp_my_color > 0) my_color = colorHeight - tmp_my_color - 1;
  if (tmp_my_depth > 0) my_depth = depthHeight - tmp_my_depth - 1;

  // Return success
  return 1;
}



int RGBDCamera::
WriteScanNetFile(const char *filename) const
{
  // Open camera file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open camera file %s\n", filename);
    return 0;
  }

  // Write camera file
  R4Matrix yzflip(1, 0, 0, 0,  0, -1, 0, 0,  0, 0, -1, 0,  0, 0, 0, 1);
  R4Matrix m = yzflip.Inverse() * depthToColorExtrinsics * yzflip.Inverse();
  if (depthToColorExtrinsics.IsIdentity()) m = R4identity_matrix;
  fprintf(fp, "colorWidth = %d\n", colorWidth);
  fprintf(fp, "colorHeight = %d\n", colorHeight); 
  fprintf(fp, "depthWidth = %d\n", depthWidth);  
  fprintf(fp, "depthHeight = %d\n", depthHeight); 
  fprintf(fp, "fx_color = %g\n", fx_color);
  fprintf(fp, "fy_color = %g\n", fy_color);
  fprintf(fp, "mx_color = %g\n", mx_color); 
  fprintf(fp, "my_color = %g\n", colorHeight - my_color - 1); 
  fprintf(fp, "fx_depth = %g\n", fx_depth); 
  fprintf(fp, "fy_depth = %g\n", fy_depth); 
  fprintf(fp, "mx_depth = %g\n", mx_depth); 
  fprintf(fp, "my_depth = %g\n", depthHeight - my_depth - 1); 
  fprintf(fp, "k1_color = %g\n", k1_color); 
  fprintf(fp, "k2_color = %g\n", k2_color); 
  fprintf(fp, "k1_depth = %g\n", k1_depth); 
  fprintf(fp, "k2_depth = %g\n", k2_depth); 
  fprintf(fp, "deviceId = %s\n", deviceId);
  fprintf(fp, "deviceName = %s\n", deviceName);
  fprintf(fp, "sceneLabel = %s\n", sceneLabel);
  fprintf(fp, "sceneType = %s\n", sceneType);
  fprintf(fp, "numDepthFrames = %d\n", numDepthFrames);
  fprintf(fp, "numColorFrames = %d\n", numColorFrames);
  fprintf(fp, "numIMUmeasurements = %d\n", numIMUmeasurements);
  fprintf(fp, "depthToColorExtrinsics = %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
     m[0][0], m[0][1], m[0][2], m[0][3],
     m[1][0], m[1][1], m[1][2], m[1][3],
     m[2][0], m[2][1], m[2][2], m[2][3],
     m[3][0], m[3][1], m[3][2], m[3][3]);

  // Close camera file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Matterport i/o functions
////////////////////////////////////////////////////////////////////////

int RGBDCamera::
ReadMatterportFile(const char *filename)
{
  // Open camera file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open camera file %s\n", filename);
    return 0;
  }

  // Read camera file
  double width, height, fx, fy, cx, cy, k1, k2, p1, p2, k3;
  if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &width, &height,
    &fx, &fy, &cx, &cy, &k1, &k2, &p1, &p2, &k3) != (unsigned int) 11) {
    fprintf(stderr, "Unable to read Matterport intrinsics matrix.\n");
    return 0;
  }

  // Close camera file
  fclose(fp);

  // Set values
  colorWidth = width;    
  colorHeight = height; 
  depthWidth = width;  
  depthHeight = height; 
  fx_color = fx;
  fy_color = fy;
  mx_color = cx; 
  my_color = height - cy - 1; 
  k1_color = k1;
  k2_color = k2;
  k3_color = k3;
  k4_color = 0;
  p1_color = p1;
  p2_color = p2;
  fx_depth = fx; 
  fy_depth = fy; 
  mx_depth = cx; 
  my_depth = height - cy - 1;
  k1_depth = k1;
  k2_depth = k2;
  k3_depth = k3;
  k4_depth = 0;
  p1_depth = p1;
  p2_depth = p2;
  depthToColorExtrinsics = R4identity_matrix;

  // Return success
  return 1;
}



int RGBDCamera::
WriteMatterportFile(const char *filename) const
{
  // Open camera file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open camera file %s\n", filename);
    return 0;
  }

  // Write intrinsics file
  fprintf(fp, "%d %d %g %g %g %g %g %g %g %g %g\n", 
    colorWidth, colorHeight, fx_color, fy_color, mx_color, my_color,
    k1_color, k2_color, p1_color, p2_color, k3_color);

  // Close camera file
  fclose(fp);

  // Return success
  return 1;
}



