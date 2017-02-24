#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "opencv_core" for configuration "Debug"
set_property(TARGET opencv_core APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_core PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_core249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG ""
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_core249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_core )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_core "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_core249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_core249d.dll" )

# Import target "opencv_flann" for configuration "Debug"
set_property(TARGET opencv_flann APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_flann PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_flann249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_flann249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_flann )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_flann "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_flann249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_flann249d.dll" )

# Import target "opencv_imgproc" for configuration "Debug"
set_property(TARGET opencv_imgproc APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_imgproc PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_imgproc249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_imgproc249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_imgproc )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_imgproc "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_imgproc249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_imgproc249d.dll" )

# Import target "opencv_highgui" for configuration "Debug"
set_property(TARGET opencv_highgui APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_highgui PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_highgui249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_imgproc"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_highgui249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_highgui )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_highgui "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_highgui249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_highgui249d.dll" )

# Import target "opencv_features2d" for configuration "Debug"
set_property(TARGET opencv_features2d APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_features2d PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_features2d249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_features2d249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_features2d )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_features2d "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_features2d249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_features2d249d.dll" )

# Import target "opencv_calib3d" for configuration "Debug"
set_property(TARGET opencv_calib3d APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_calib3d PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_calib3d249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_calib3d249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_calib3d )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_calib3d "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_calib3d249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_calib3d249d.dll" )

# Import target "opencv_ml" for configuration "Debug"
set_property(TARGET opencv_ml APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_ml PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_ml249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_ml249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_ml )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_ml "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_ml249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_ml249d.dll" )

# Import target "opencv_video" for configuration "Debug"
set_property(TARGET opencv_video APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_video PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_video249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_imgproc"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_video249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_video )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_video "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_video249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_video249d.dll" )

# Import target "opencv_legacy" for configuration "Debug"
set_property(TARGET opencv_legacy APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_legacy PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_legacy249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_legacy249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_legacy )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_legacy "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_legacy249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_legacy249d.dll" )

# Import target "opencv_objdetect" for configuration "Debug"
set_property(TARGET opencv_objdetect APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_objdetect PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_objdetect249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_imgproc;opencv_highgui"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_objdetect249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_objdetect )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_objdetect "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_objdetect249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_objdetect249d.dll" )

# Import target "opencv_photo" for configuration "Debug"
set_property(TARGET opencv_photo APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_photo PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_photo249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_imgproc"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_photo249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_photo )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_photo "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_photo249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_photo249d.dll" )

# Import target "opencv_gpu" for configuration "Debug"
set_property(TARGET opencv_gpu APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_gpu PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_gpu249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video;opencv_legacy;opencv_objdetect;opencv_photo"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_gpu249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_gpu )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_gpu "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_gpu249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_gpu249d.dll" )

# Import target "opencv_ocl" for configuration "Debug"
set_property(TARGET opencv_ocl APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_ocl PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_ocl249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_objdetect;opencv_video"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_ocl249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_ocl )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_ocl "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_ocl249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_ocl249d.dll" )

# Import target "opencv_nonfree" for configuration "Debug"
set_property(TARGET opencv_nonfree APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_nonfree PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_nonfree249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video;opencv_legacy;opencv_objdetect;opencv_photo;opencv_gpu;opencv_ocl"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_nonfree249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_nonfree )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_nonfree "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_nonfree249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_nonfree249d.dll" )

# Import target "opencv_contrib" for configuration "Debug"
set_property(TARGET opencv_contrib APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_contrib PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_contrib249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video;opencv_legacy;opencv_objdetect;opencv_photo;opencv_gpu;opencv_ocl;opencv_nonfree"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_contrib249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_contrib )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_contrib "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_contrib249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_contrib249d.dll" )

# Import target "opencv_stitching" for configuration "Debug"
set_property(TARGET opencv_stitching APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_stitching PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_stitching249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video;opencv_legacy;opencv_objdetect;opencv_photo;opencv_gpu;opencv_ocl;opencv_nonfree"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_stitching249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_stitching )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_stitching "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_stitching249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_stitching249d.dll" )

# Import target "opencv_superres" for configuration "Debug"
set_property(TARGET opencv_superres APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_superres PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_superres249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video;opencv_legacy;opencv_objdetect;opencv_photo;opencv_gpu;opencv_ocl"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_superres249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_superres )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_superres "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_superres249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_superres249d.dll" )

# Import target "opencv_ts" for configuration "Debug"
set_property(TARGET opencv_ts APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_ts PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_ts249d.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_ts )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_ts "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_ts249d.lib" )

# Import target "opencv_videostab" for configuration "Debug"
set_property(TARGET opencv_videostab APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(opencv_videostab PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_videostab249d.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "opencv_core;opencv_flann;opencv_imgproc;opencv_highgui;opencv_features2d;opencv_calib3d;opencv_ml;opencv_video;opencv_legacy;opencv_objdetect;opencv_photo;opencv_gpu"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_videostab249d.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS opencv_videostab )
list(APPEND _IMPORT_CHECK_FILES_FOR_opencv_videostab "${_IMPORT_PREFIX}/x64/vc12/lib/opencv_videostab249d.lib" "${_IMPORT_PREFIX}/x64/vc12/bin/opencv_videostab249d.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
