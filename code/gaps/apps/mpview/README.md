# mpview

mpview is a C++ application for parsing and viewing houses in the Matterport3D dataset.   The code for parsing input files can be found in mp.cpp.   The code for interactive visualization is mostly in mpview.cpp.  Please feel free to copy/use any of this code (including gaps/pkgs) in your own applications.

## Usage

The mpview application can be invoked by simplying typing "mpview -input_house house_segmentations/xxx.house -v".   However, it also supports control over how input files are parsed and how data is displayed on startup as follows:   

    Usage: mpview [options]
      -input_house <filename> : input region segmentations file (e.g., xxx/region_segmentations/xxx.house)
      -input_scene <filename> : input textured mesh (e.g., xxx/matterport_mesh/*/*.obj)
      -input_mesh <filename> : input poisson mesh (e.g., xxx/poisson_meshes/xxx_11.ply)
      -input_categories <filename> : input categories tsv file (e.g., metadata/*.tsv)
      -input_segments <filename> : input json file with face segments (e.g., xxx/object_segmentations/*.fsegs.json)
      -input_objects <filename> : input json file with objects and labels (e.g., xxx/object_segmentations/*.semseg.json)
      -input_configuration <filename> : input file with images and panorama (e.g., xxx/undistorted_camera_parameters/xxx.conf)
      -output_image <filename> : save an image to <filename> and exit
      -background <r> <g> <b> : background color (with each component in [0.0-1.0])
      -window <width> <height> : window size in pixels
      -camera <ex> <ey> <ez> <tx> <ty> <tz> <ux> <uy> <uz> : initial camera extrinsics
      -batch : exit without starting interactive viewer
      -v : print verbose (recommended)
    
    Typical usage for viewing house segmentations:
      cd scans/17DRP5sb8fy (or any other house)
      mpview -input_house house_segmentations/*.house -v   OR 
      mpview -input_house house_segmentations/*.house -input_scene matterport_mesh/*/*.obj -v   OR
      mpview -input_house house_segmentations/*.house -input_mesh house_segmentations/*.ply -v
    
    Typical usage for viewing region segmentations:
      cd scans/17DRP5sb8fy (or any other house)
      mpview -input_mesh region_segmentations/region0.ply -input_segments region_segmentations/region0.fsegs.json -input_objects region_segmentations/region0.semseg.json -input_categories ../../metadata/category_mapping.tsv -v


## Command interface

Upon invoking the program, a window will popup (unless the -batch command is given on the command line).   By default, mpview will show oriented bounding boxes for objects colored by their semantic labels.  You may interactively change the view by dragging the mouse, change what is shown by keyboard strokes, change the coloring scheme by hitting the space bar, etc.   Here is an overview of the available commands:

    Camera control:
       Left-button-drag = rotate the camera view 
       Middle-button-drag = zoom the camera view
       Right-button-drag = pan the camera view
       Scroll-wheel = zoom the camera view
       Left-click = set the center of rotation and zooming
       PAGE UP/DOWN = match the camera parameters to the next/prev image
    
    Query commands:
       Left-button-click = select a region, object, or image
       Left-button-double-click = print information about what is under cursor to stdout
    
    Display options:
       A = toggle display of Cartesian axes
       B = toggle display of bounding boxes
       C = toggle display of image viewpoints
       E = toggle display of edges (wireframe)
       F = toggle display of faces
       H = print this list of commands
       I = toggle display of images
       L = toggle display of text labels
       M = toggle display of mesh
       O = toggle display of objects
       P = toggle display of panorama centers
       R = toggle display of regions
       S = toggle display of scene
       V = toggle display of vertices
       X = toggle display of clip box (only clips if visible)
       SPACE = cycles through different colors schemes
    
    Clipbox control:
       ARROW UP/DOWN = move top of clipbox 
       ARROW LEFT/RIGHT = move bottom of clipbox
       ESC = reset the clipbox (and remove all selections)
    
    Quit:
       Ctrl-Q = quit 
    
