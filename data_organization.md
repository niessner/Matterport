Data organization
=========================

The main dataset resides in the "data" directory.  There is a separate subdirectory for every house, which is named by a unique string (e.g., "1pXnuDYAj8r").  Within each house directory, there are separate directories for different types of data as follows:

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
- poisson_meshes = low-resolution meshes resulting from poisson mesh reconstruction

Manually specified annotations:
- cameras = camera extrinsics for manually chosen good view(s)
- house_segmentations = manually specified semantic segmentations of houses into regions
- region_segmentations = manually specified semantic segmentations of regions into objects 

Details about each of the data directories follow.



General file naming conventions for matterport images:
---------------------

The Matterport Pro camera consists of three Primesense Carmine depth and RGB sensors placed on an eye-height tripod.   The
three cameras (camera index 0-2) are oriented diagonally up, flat, and diagonally down, cumulatively covering most of the vertical field of view. During capture, the camera rotates  around its vertical axis and stops at 60-degree intervals (yaw index 0-5). So, for each tripod placement (panorama_uuid), a total of 18 images are captured.  Files associated with these images adhere to the following naming convention:

    <panorama_uuid>_<imgtype><camera_index>_<yaw_index>.<extension>

where <panorama_uuid> is a unique string, <camera_index> is [0-5], and <yaw_index> is [0-2].  <imgtype> is 'j' for HDR images, 'i' for tone-mapped color images, 'd' for depth images, "skybox" for skybox images, "pose" for camera pose files, and "intrinsics" for camera intrinsics files.  The extension is ".jxr" for HDR images, ".jpg" for tone-mapped color images, and ".png" for depth and normal images.


matterport_hdr_images
---------------------

Raw HDR images in jxr format.   

Each RGB channel provides a 16-bit value (jxr_val[i]) that can be mapped to intensity (col_val[i]) with:

    for (int i = 0; i < 3; i++) {
      if (jxr_val[i] <= 3000) col_val[i] = jxr_val[i] * 8e-8
      else col_val[i] = 0.00024 * 1.0002 ** (jxr_val[i] - 3000)
    }
    
    col_val[0] *= 0.8
    col_val[1] *= 1.0
    col_val[2] *= 1.6
    
    
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
Textured mesh for the entire house.  The subdirectory contains a single .obj file, a single .mtl file, and textures in .jpg and .png format as referenced by the mtl file.


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


undistorted_depth_images 
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

The files are binary PLY format.  xxx.ply provides a mesh that has been simplified to a reasonable polygon count (e.g., 500K) and thus provides a plausible trade-off between accuracy and size.  Additionally, three meshes are provided corresponding to different depth levels of the poisson surface reconstruction (e.g., xxx_9.ply corresponds to depth 9, xxx_10.ply corresponds to depth 10, etc.).   Higher depths correspond to finer resolution meshes.    

These meshes are generally lower resolution than the ones in house_segmentation/xxx.ply, which were generated by splitting the house into regions, reconstructing a mesh for each region, and then merging into one mesh.


cameras
---------------------
Camera extrinsics for manually chosen good view(s).

    exterior.cam - manually chosen camera viewpoints to view houses from a bird's eye view

Each .cam file has one line per camera with ascii numbers indicating the following camera parameters separated by spaces:

    vx vy vz  tx ty tz  ux uy uz  xfov yfov 1

where (vx, vy, vz) is the eye viewpoint of the camera, (tx, ty, tz) is the view direction, (ux, uy, uz) is the up direction, and xfov and yfov are the half-angles of the horizontal and vertical fields of view of the camera in radians (the angle from the central ray to the leftmost/bottommost ray in the field of view).


house_segmentations
---------------------
A manually specified decomposition of a house into levels, room-like regions, and objects with semantic labels.

The data for each house xxx is represented in four files:

    xxx.house = ascii file containing a list of regions, objects, etc.
        The .house file has a sequence of ascii lines with fields separated by spaces in the following format:

        H name label #images #panoramas #vertices #surfaces #segments #objects #categories #regions #portals #levels  0 0 0 0 0  xlo ylo zlo xhi yhi zhi  0 0 0 0 0
        L level_index #regions label  px py pz  xlo ylo zlo xhi yhi zhi  0 0 0 0 0
        R region_index level_index 0 0 label  px py pz  xlo ylo zlo xhi yhi zhi  height  0 0 0 0
        P portal_index region0_index region1_index label  xlo ylo zlo xhi yhi zhi  0 0 0 0
        S surface_index region_index 0 label px py pz  nx ny nz  xlo ylo zlo xhi yhi zhi  0 0 0 0 0
        V vertex_index surface_index label  px py pz  nx ny nz  0 0 0
        P name  panorama_index region_index 0  px py pz  0 0 0 0 0
        I image_index panorama_index  name camera_index yaw_index e00 e01 e02 e03 e10 e11 e12 e13 e20 e21 e22 e23 e30 e31 e32 e33  i00 i01 i02  i10 i11 i12 i20 i21 i22  width height  px py pz  0 0 0 0 0
        C category_index category_mapping_index category_mapping_name mpcat40_index mpcat40_name 0 0 0 0 0
        O object_index region_index category_index px py pz  a0x a0y a0z  a1x a1y a1z  r0 r1 r2 0 0 0 0 0 0 0 0 
        E segment_index object_index id area px py pz xlo ylo zlo xhi yhi zhi  0 0 0 0 0
   
        where xxx_index indicates the index of the xxx in the house file (starting at 0),
        #xxxs indicates how many xxxs will appear later in the file that back reference (associate) to this entry,
        (px,py,pz) is a representative position, (nx,ny,nz) is a normal direction,
        (xlo, ylo, zlo, xhi, yhi, zhi) is an axis-aligned bounding box,
        camera_index is in [0-5], yaw_index is in [0-2],
        (e00 e01 e02 e03 e10 e11 e12 e13 e20 e21 e22 e23 e30 e31 e32 e33) are the extrinsic matrix of a camera,
        (i00 i01 i02  i10 i11 i12 i20 i21 i22) are the intrinsic matrix for a camera,
        (px, py, pz, a0x, a0y, a0z, a1x, a1y, a1z, r0, r1, r2) define the center, axis directions, and radii of an oriented bounding box, 
        height is the distance from the floor, and 
        0 is a value that can be ignored.

        The extent of each region is defined by a prism with its vertical extent dictated by its height and
        its horizontal cross-section dictated by the counter-clockwise set of polygon vertices associated
        with each surface assocated with the region.

        The extent of each object is defined by the oriented bounding box of the 'O' command.
        The set of faces associated with each segment are ones whose 'face_material' field
        in the xxx.ply file (described next) matches the segment 'id' in the 'S' command.

    xxx.ply = 3D mesh in ply format.   In addition to the usual fields, there are three additional fields for each face:
        face_material = unique id of segment containing this face
        face_segment = unique id of object instance containing this face
        face_category = unique id of the category label for the object instance containing this face 
            (i.e., mapping to the "index" column of the category.tsv file)
        
    xxx.fsegs.json = JSON file indicating which segment contains each face of the mesh in regionX.ply, where
        segIndices = an array of unique segment IDs, one per face in the order it appears in the mesh 
               (i.e., the Kth entry provides the unique segment ID for the Kth face of the mesh)   
               
    xxx.semseg.json = JSON file containing an array of object instances with category labels, where
        segGroups = an array of object instances, each with the following fields
        label = string indicating the raw label provided by a human annotator for that object instance
                (this label maps to the "raw category" column of categories.tsv)
        segments = an array containing the unique ids for all segments in this object instance
                (the unique ids for segments map to ones found in regionX.fsegs.json)

The 'label' of each region is a string with the following conventions:

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

The label of each object is defined by its category index.   For each category, the .house file provides a category_mapping_index, category_mapping_name, mcat40_index, and mcat40_name.   The category_mapping_index maps into the first column of [metadata/category_mapping.tsv](metadata/category_mapping.tsv).  Further information these object categories (including how they map to WordNet synsets) can be
found in [metadata/category_mapping.tsv](metadata/category_mapping.tsv).  The mpcat40_index maps into the first column of [metadata/mpcat40.tsv](metadata/mpcat40.tsv), which provides further information about them (including a standard color to display for each one).


region_segmentations
---------------------
A set of manually specified segment, object instance, and semantic category labels for walls, floors, ceilings, doors, windows, and "furniture-sized" objects for each region of each house. 

The region_segmentations and labels are provided as annotations on 3D meshes.   A ply file provides the raw geometry for each region.   Json files indicate how each triangle of the mesh is associated with a "segment", how segments are associated with object instances, and, how object instances are associated with semantic categories as follows (this is the same as for house_segmentations):

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


