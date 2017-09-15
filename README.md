# Matterport3D
The Matterport3D V1.0 dataset contains data captured throughout 90 properties
with a Matterport Pro Camera.   

This repository includes the raw data for the dataset plus derived data, annotated 
data, and scripts/models for several scene understanding tasks.

Downloading:
------------
You must indicate that you agree to the terms of use as specified in the dataset [license](http://dovahkiin.stanford.edu/matterport/public/MP_TOS.pdf) by sending an email to the mailing list at: `matterport3d@googlegroups.com`.  We will then provide download access to the dataset.

File naming conventions:
---------------------

The Matterport Pro camera consists of three Primesense Carmine
depth and RGB sensors placed on an eye-height tripod.   The
three cameras (camera index 0-2) are oriented diagonally up, flat, and diagonally down,
cumulativley covering most of the vertical field of view.
During capture, the camera rotates 
around its vertical axis and stops at 60-degree intervals (yaw index 0-5).
So, for each tripod placement (panorama_uuid), a total of 18 
images are captured.  Files associated with these images adhere to the 
following naming convention:

<panorama_uuid>_<prefix><camera_index>_<yaw_index>.<extension>

where <panorama_uuid> is a unique string, <camera_index> is [0-5], and <yaw_index> is [0-2].  <prefix> is 'j' for HDR images, 'i' for tone-mapped color images, 'd' for depth images, "skybox" for skybox images, "pose" for camera pose files, and "intrinsics" for camera intrinsics files.  The extension is ".jxr" for HDR images, ".jpg" for tone-mapped color images, and ".png" for depth and normal images.


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
- undistorted_camera_parameters = camera intrinsics and extrinsics after undistortion
- undistorted_color_images = color images after undistortion
- undistorted_depth_images = depth images after undistortion
- undistorted_normal_images = normal and boundary images aligned with undistorted depth images
- poisson_meshes = meshes resulting from poisson mesh reconstruction

Manually specified annotations:
- cameras = camera extrinsics for manually chosen good view(s)
- region_segmentations = manually specified floorplans with instance and semantic category labels
- object_segmentations = manually specified object instance and semantic category labels for surfaces of a 3D mesh

Task specific data:
- view_overlap_data = computed overlaps (loop closures) between pairs of views within the same house 

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


undistorted_camera_parameters
---------------------

An ascii file indicating the intrinsic and extrinsic camera parameters for every image.   Each line of the file is a separate command with the following format, where <n> is the number of scan commands (imaages) in the file, fx and fy are focal lengths in pixels, cx and cy specific the optical center of the image in pixels (not necessarily the center of the image), depth and color image filenames are relative to the depth and color directories, and camera-to-world-matrix is 16 values in row-major order providing a matrix that takes a canonical camera into a right-handed world coordinate system (the inverse of a typical extrinsics matrix).   The canonical camera has its eye at the origin, its view direction at (0,0-1), and its up direction at (0,1,0).   Note that all images have their origins in the bottom-left of the image.

    dataset matterport
    n_images <n>
    depth_directory undistorted_depth_images
    color_directory undistorted_color_images
    intrinsics_matrix <fx> 0 <cx>  0 <fy> <cy> 0 0 1
    scan <depth_image_filename> <color_image_filename> <camera-to-world-matrix>


undistorted_color_images 
---------------------

Tone-mapped color images after undistortion.   Radial distortion is removed, but (cx,cy) is not shifted to the center of the image.  Files are in JPG format.  Note that all images have their origins in the bottom-left of the image.


undistorted_depth_images = depth images after undistortion
---------------------

Depth images after undistortion.   Pixels of these depth images should (approximately) map to corresponding pixels in the undistorted_color_images.   The files are in 16-bit PNG format with the same scaling as matterport_depth_images (0.25mm per unit).  Note that all images have their origins in the bottom-left of the image.


undistorted_normal_images
---------------------

Normal, boundary, and radius maps estimated from the undistorted_depth_images.  

Normal maps are stored in three 16-bit PNG files (_nx.png, _ny.png, and _nz.png), where the integer values in the file are 32768*(1 + n), where n is a normal coordinate in range [-1,1].

Boundary maps are stored in a 16-bit PNG file (_boundary.png), where each pixel has one of the following values:
1 = on the border, 
2 = on a silhouette edge (boundary of an occluder), and
4 = on a shadow edge (boundary for occluded surface ... must be adjacent to a silhouette edge).
 
Radius maps are stored in a 16-bit PNG file (_radius.png), where integer values in the file are 4000 times the "radius" (average distance) to the pixel's neighbors in 3D.

Note that all images have their origins in the bottom-left of the image.


poisson_meshes
---------------------
Surface meshes reconstructed from the depth images using [Screened Poisson Surface Reconstruction](http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version9.01/).
The files are binary PLY format. 


cameras
---------------------
Camera extrinsics for manually chosen good view(s).

    exterior.cam - manually chosen camera viewpoints to view houses from a bird's eye view

Each .cam file has one line per camera with ascii numbers indicating the following camera parameters separated by spaces:

    vx vy vz  tx ty tz  ux uy uz  xfov yfov 1

where (vx, vy, vz) is the eye viewpoint of the camera, (tx, ty, tz) is the view direction, (ux, uy, uz) is the up direction, and xfov and yfov are the half-angles of the horizontal and vertical fields of view of the camera in radians (the angle from the central ray to the leftmost/bottommost ray in the field of view).


house_floorplans
---------------------
A list of manually specified floor and region boundaries along with semantic regon labels.

Each .house file has a sequence of ascii lines with fields separated by spaces in the following format:

    H name label #images #panoramas #vertices #surfaces #regions #levels  0 0 0 0 0 0 0 0
    L level_index #regions label  px py pz  xlo ylo zlo xhi yhi zhi  0 0 0 0 0
    R region_index level_index #panoramas #surfaces label  px py pz  xlo ylo zlo xhi yhi zhi  0 0 0 0 0
    S surface_index region_index #vertices #surfaces label  px py pz  nx ny nz  xlo ylo zlo xhi yhi zhi 0 0 0 0 0
    V vertex_index surface_index label  px py pz  nx ny nz  0 0 0
    P name panorama_index region_index #images  px py pz  0 0 0 0 0
    I name panorama_index panorama_index  px py pz  0 0 0 0 0
   
where xxx_index indicates the index of the xxx in the house file (starting at 0), #xxxs indicates how many xxxs will appear later in the file that back reference (associate) to this entry, (px,py,pz) is a representative position, (nx,ny,nz) is a normal direction, and (xlo, ylo, zlo, xhi, yhi, zhi) is an axis-aligned bounding box, and 0 is a value that can be ignored.   The extent of each region is defined by a prism with its vertical extent dictated by zlo and zhi as its horizontal cross-section dictated by the counter-clockwise set of polygon vertices associated with the first surface assocated with the region.  

The label of each region is a string with the following conventions:

    'a' = bathroom (should have a toilet and a sink)
    'b' = bedroom
    'c' = closet
    'd' = dining room (includes “breakfast rooms” other rooms people mainly eat in)
    'e' = entryway/foyer/lobby (should be the front door, not any door)
    'f' = familyroom (should be a room that a family hangs out in, not any area with couches)
    'g' = garage
    'h' = hallway
    'i' = library (should be room like a library at a university, not an individual study)
    'j' = laundryroom/mudroom (place where people do laundry, etc.)
    'k' = kitchen
    'l' = living room (should be the main “showcase” living room in a house, not any area with couches)
    'm' = meetingroom/conferenceroom
    'n' = lounge (any area where people relax in comfy chairs/couches that is not the family room or living room
    'o' = office (usually for an individual, or a small set of people)
    'p' = porch/terrace/deck/driveway (must be outdoors on ground level)
    'r' = rec/game (should have recreational objects, like pool table, etc.)
    's' = stairs
    't' = toilet (should be a small room with ONLY a toilet)
    'u' = utilityroom/toolroom 
    'v' = tv (must have theater-style seating)
    'w' = workout/gym/exercise
    'x' = outdoor areas containing grass, plants, bushes, trees, etc.
    'y' = balcony (must be outside and must not be on ground floor)
    'z' = other room (it is clearly a room, but the function is not clear)
    'B' = bar
    'C' = classroom
    'D' = dining booth
    'S' = spa/sauna
    'Z' = junk (reflections of mirrors, random points floating in space, etc.)
    '-' = no label 
    
    
object_segmentations
---------------------
A set of manually specified segment, object instance, and semantic category labels for walls, floors, ceilings, doors, windows, and "furniture-sized" objects. 

The meshes and annotations are split into regions for each of processing.  The labels are provided as annotations on 3D meshes.   A ply file provides the raw geometry for each region.   Json files indicate how each triangle of the mesh is associated with a "segment", how segments are associated with object instances, and, how object instances are associated with semantic categories as follows:

    regionX.ply = 3D mesh in ply format.   In addition to the usual fields, there are three additional fields for each face:
        face_material = unique id of segment containing this face
        face_segment = unique id of object instance containing this face
        face_category = unique id of the category label for the object instance containing this face 
            (i.e., mapping to the "index" column of the category.tsv file)
        
    regionX.fsegs.json = JSON file indicating which segment contains each face of the mesh in regionX.ply, where
        segIndices = an array of unique segment IDs, one per face in the order it appears in the mesh 
               (i.e., the Kth entry provides the unique segment ID for the Kth face of the mesh)   
               
    regionX.semseg.json = JSON file containing an array of object instances with category labels, where
        segGroups = an array of object instances, each with the following fields
            label = string indicating the raw label provided by a human annotator for that object instance
                (this label maps to the "raw category" column of categories.tsv)
            segments = an array containing the unique ids for all segments in this object instance
                (the unique ids for segments map to ones found in regionX.fsegs.json)


view_overlap_data
---------------------
Statistics about how images "overlap" with one another.  There are three main files. The first (xxx_iis.txt) has a line for every pair of images whose visible surfaces overlap.   The second (xxx_vi.txt) has a line for a sampling of vertices on the specified surface mesh indicating which images see it. The third (xxx_iv.txt) has a line for each image indicating a sampling of vertices of the mesh are seen by it.   The other files (xxx_iis_iou.png and (xxx_iis_count.png) have the same information in more convenient matrix forms.  Details of the file formats follow.

xxx_iis.txt

    C matterport_camera_poses/xxx.conf
    II iid1 iid2 iou isect union, count1 count2
    followed by many lines starting with 'II' in the same format

        where II is a keyword (it is capital i twice)
        iid1 and iid2 are zero-based integer indices of images in the conf file,
        count1 is the number of pixels in image1 that see surface points visible in image2,
        count2 is the number of pixels in image2 that see surface points visible in image1,
        isect is min(count1, count2),
        union is count0+count1-isect, and
        iou is isect/union

xxx_vi.txt

    C matterport_camera_poses/xxx.conf
    M poisson_meshes/xxx.ply
    VI vid px py pz nx ny nz n iid1 iid2 ...
    followed by many lines starting with 'VI' in the same format

        where VI is a keyword
        where vid is a zero-based integer index of a vertex in the mesh file
        (px,py,pz) is the position of the vertex,
        (nx,ny,nz) is the normal of the vertex,
        n is the number of images that see it, and
        iid1, iid2, ... is a list of n zero-based image indices that see the vertex

xxx_iv.txt

    C matterport_camera_poses/xxx.conf
    M poisson_meshes/xxx.ply
    IV iid n vid1 vid2 ...
    followed by many lines starting with 'IV' in the same format

        where IV is a keyword
        iid is a zero-based integer index of an image in the conf file
        n is the number of vertices seen by it, and
        vid1, vid2, ... is a list of n zero-based vertex indices seen in the image

xxx_iis_iou.png

    An NxN 16-bit png file

        where N is the number of images in the conf file, and
        pixel (i,j) indicates "intersection over union" (iou) of pixels
        in images i, j multiplied by 1000.

xxx_iis_count.pfm

    An NxN 16-bit pfm file

        where N is the number of images in the conf file, and pixel (i,j) indicates the
        number of pixels in image i that see surface points visible in image j.
        The PFM file format has a header line with "Pf" in ASCII, a second line with
        "<width> <height>" in ASCII, and a third line with "-1.0" in ASCII, followed by
        <width>*<height> 4-byte floats in row-major order (i varies more quickly).

Actually, there are four versions of the image-image overlap files (iis, iip, iig, and iiv),
representing different ways in which overlaps were computed.   In each of the first three files,
overlaps between images i and j are counted by looping over pixels of image i and counting
the number of times the ...

  iis - 3D backprojection of pixel in image i is within 5cm of
        the 3D backprojection of an any pixel in image j.
  iig - 3D backprojection of pixel in image i i falls into the
        same 5-cm-wide grid cell as the backprojection of any pixel in image j.
  iip - projection from pixel in image i to 3D and the back into image j falls on
        a pixel in j with depth +/- 10% of observed depth at that pixel.

In the fourth file (iiv), overlaps between image i and j are counted by looping over
vertices of a mesh and counting the number of vertices visible in both images

  iiv - a vertex of a mesh visible at pixel i is also visible in image j


Benchmark Task Data
=======================

**Image Keypoint Matching**  
The image keypoint matching task aims to establish correspondences between keypoints in RGB image data.   It leverages the wide variety of camera baselines in the Matterport3D dataset as training and test data. Please see [keypointmatch] (keypointmatch).

**View Overlap Prediction**  
The view overlap prediction task aims to predict how much the views of two images overlap (what fraction of the visible surfaces are shared between the views).  It leverages the wide variety of camera baselines in the Matterport3D dataset as training and test data.  Please check [`https://github.com/niessner/Matterport/tree/master/view_overlap`](./view_overlap) for train/test codes, pretrained models, and auxiliary data for the experiments.  Please see [`here`](./view_overlap/readme.md) for how to download  the data and run the scripts for this task.

**Surface Normal Estimation**  
The surface normal estimation task aims toi predict pixelwise surface normals from RGB images.   It leverages normals estimated from the vast number of RGB-D image pairs in the Matterport3D dataset as training and testing data.  Please check [`https://github.com/niessner/Matterport/tree/master/surface_normal`](./surface_normal) for train/test codes, pretrained models, and auxiliary data for the experiments. Note that to run the experiments for surface normal estimation, you don't need to download the whole dataset. Please see [`here`](./surface_normal/readme.md) for how to download the data and run the scripts for this task.

**Semantic Voxel Labeling**  
The semantic voxel labeling task predicts per-voxel class labels for a scan. Please see [semantic_voxel_label](semantic_voxel_label).

**Room-Type Categorization**  
The room type categorization task aims to predict the semantic category of the region (e.g., bedroom, kitchen, patio, etc.) containing the camera viewpoint of an RGB image or panorama.   It leverages semantic boundaries and labels for manually-specified regions in the Matterport3D dataset.  Please check [`https://github.com/niessner/Matterport/tree/master/room_categorization`](./room_categorization) for train/test codes, pretrained models, and auxiliary data for the experiments.  Please see [`here`](./room_categorization/readme.md) for how to download the data and run the scripts for this task.

