# Benchmark Tasks

Using the Matterport3D dataset, we present several benchmarking tasks.  For each task, all train/test code, pretrained models, and auxiliary data for the experiments are provided.

## Tasks

### Image Keypoint Matching
The image keypoint matching task aims to establish correspondences between keypoints in RGB image data.   It leverages the wide variety of camera baselines in the Matterport3D dataset as training and test data. Please see [keypoint_match](keypoint_match) for details.

### View Overlap Prediction
The view overlap prediction task aims to predict how much the views of two images overlap (what fraction of the visible surfaces are shared between the views).  It leverages the wide variety of camera baselines in the Matterport3D dataset as training and test data.  Please see [view_overlap](view_overlap) for details.

### Surface Normal Estimation
The surface normal estimation task aims to predict pixelwise surface normals from RGB images.   It leverages normals estimated from the vast number of RGB-D image pairs in the Matterport3D dataset as training and testing data.  Please see [surface_normal](surface_normal) for details.

### Region Type Classification  
The region classification task aims to predict the semantic category of the region (e.g., bedroom, kitchen, patio, etc.) containing the camera viewpoint of an RGB image or panorama.   It leverages semantic boundaries and labels for manually-specified regions in the Matterport3D dataset.  Please see [region_classification](region_classification) for details.

### Semantic Voxel Labeling
The semantic voxel labeling task predicts per-voxel class labels for a scan. Please see [semantic_voxel_label](semantic_voxel_label) for details.


