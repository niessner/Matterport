# Matterport3D
The Matterport3D V1.0 dataset contains data captured throughout 90 properties
with a Matterport Pro Camera.   The camera consists of three Primesense Carmine
depth and RGB sensors, oriented in down, middle, and up
orientations, and placed on an eye-height tripod.  During capture, the camera around its vertical
axis and stops at 60-degree intervals in order to capture 18
multi-exposure HDR images per panoramic sweep.  

This repository includes the raw data for the dataset plus derived data,annotated 
data, and scripts/models for several scene understanding tasks.

File naming conventions:
---------------------
The files are organized by camera, panorama, and yaw with the following naming convention:

<panorama_uuid>_<prefix><camera_index>_<yaw_index>.<extension>

where <panorama_uuid> is a unique string, <camera_index> is [0-5], and <yaw_index> is [0-2].  <prefix> is 'j' for HDR images, 'i' for tone-mapped color images, 'd' for depth images, "skybox" for skybox images, "pose" for camera pose files, and "intrinsics" for camera intrinsics files.  


Data organization
=========================

The main dataset resides in the "data" directory.  There is 
a separate subdirectory for every property, which is named by a unique string (e.g., "1pXnuDYAj8r").   
Within each property directory, there are separate directories for different types 
of data as follows:

Data provided directly by Matterport:
- matterport_camera_intrinsics = camera intrinsics provided by matterport
- matterport_camera_poses = camera extrinsics provided by matterport
- matterport_color_images = color images provided by matterport
- matterport_depth_images = depth images provided by matterport
- matterport_hdr_images = HDR images provided by matterport
- matterport_mesh = textured mesh provided by matteport
- matterport_skybox_images = skybox images provided by matteport

Data computed for the convenience of future research:
- undistorted_camera_intrinsics = camera intrinsics after undistortion
- undistorted_color_images = color images after undistortion
- undistorted_depth_images = depth images after undistortion
- undistorted_normal_images = normal and boundary images aligned with undistorted depth images
- poisson_meshes = meshes resulting from poisson mesh reconstruction
- vh_meshes = meshes resulting from voxel hashing
- overlaps = statistics about how much each pair of images overlaps

Manually specified annotations:
- cameras = camera extrinsics for manually chosen good view(s)
- houses = manually drawn floorplans and room labels

Details about each of these data directories follow.


matterport_hdr_images
---------------------

Raw HDR images in jxr format.


matterport_color_images
---------------------

Tone-mapped color images in jpg format.


matterport_depth_images
---------------------

Reprojected depth images aligned with the
color images.  Every depth image a 16 bit
PNG containing the pixel's distance in the z-direction from the
camera center (_not_ the euclidean distance from the camera
center), 0.25 mm per value (divide by 4000 to get meters).
A zero value denotes 'no reading'.


matterport_camera_intrinsics
---------------------

Intrinsic parameters for every camera, stored as ASCII in the following format (using OpenCV's notation):

width height fx fy cx cy k1 k2 p1 p2 k3

The Matterport camera axis convention is:
    x-axis: right
    y-axis: down
    z-axis: "look"

The image origin is at the top-left of the image.

Sample code to map from coordinates without and with radial distortion follows:

    double nx = (undistorted_position_x - cx) / fx;
    double ny = (undistorted_position_y - cy) / fy;
    double rr = nx*nx + ny*ny;
    double rrrr = rr*rr;
    double s = 1.0 + rr*k1 + rrrr*k2 + rrrr*rr*k3 + rrrr*rrrr*k4;
    nx = s*nx + p2*(rr + 2*nx*nx) + 2*p1*nx*ny;
    ny = s*ny + p1*(rr + 2*ny*ny) + 2*p2*nx*ny;
    double distorted_position_x = nx*fx + cx;
    double distorted_position_y = ny*fy + cy;

*Note that camea intrinsics are the same across different yaw_indices, so yaw_index is omitted in their filenames*


matterport_camera_poses 
---------------------

Camera pose files contain a 4x4 matrix that transforms column
vectors from camera to global coordinates.  Translation units are
in meters.  Global coordinate system is arbitrary, but the z-axis
generally points upward if the camera was level during capture.


matterport_mesh
---------------------
Textured mesh for the entire property.  The subdirectory contains a single .obj file, a single .mtl file, and textures in .jpg and .png format as referenced by the mtl file.

undistorted_camera_intrinsics
---------------------

Camera intrinsics after undistortion.   Files are in the same format as the matterport camera intrinsics.   Note that the last five numbers (k1 k2 k3 p1 p2) are zero.

undistorted_color_images 
---------------------

Tone-mapped color images after undistortion.   Radial distortion is removed, but (cx,cy) is not shifted to the center of the image.  Files are in JPG format.


undistorted_depth_images = depth images after undistortion
---------------------

Depth images after undistortion.   Pixels of these depth images should (approximately) map to corresponding pixels in the undistorted_color_images.   The files are in 16-bit PNG format with the same scaling as matterport_depth_images (0.25mm per unit).


undistorted_normal_images
---------------------

Normal, boundary, and radius maps estimated for the undistorted_depth_images.  Normals are stored in three 16-bit PNG files ... 


poisson_meshes
---------------------
Surface meshes reconstructed from the depth images using [Screened Poisson Surface Reconstruction](http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version9.01/).
The files are binary PLY format. 

vh_meshes
---------------------
Surface meshes reconstructed from the depth images using [VoxelHashing](https://github.com/niessner/VoxelHashing).
The files are binary PLY format. 


Benchmark Task Data
=======================

ODO list files

**Keypoint Matching**  
TODO brief description and link to own readme for how to use code

**View Overlap Prediction**  
TODO brief description and link to own readme for how to use code

**Surface Normal Estimation**  
TODO brief description and link to own readme for how to use code
