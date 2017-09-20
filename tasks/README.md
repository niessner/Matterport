# Benchmark Tasks

## Image Keypoint Matching
The image keypoint matching task aims to establish correspondences between keypoints in RGB image data.   It leverages the wide variety of camera baselines in the Matterport3D dataset as training and test data. Please check [keypoint_match](keypoint_match).

## View Overlap Prediction
The view overlap prediction task aims to predict how much the views of two images overlap (what fraction of the visible surfaces are shared between the views).  It leverages the wide variety of camera baselines in the Matterport3D dataset as training and test data.  Please check [view_overlap](view_overlap) for train/test codes, pretrained models, and auxiliary data for the experiments.  Please see [here](view_overlap/readme.md) for how to download  the data and run the scripts for this task.

## Surface Normal Estimation
The surface normal estimation task aims toi predict pixelwise surface normals from RGB images.   It leverages normals estimated from the vast number of RGB-D image pairs in the Matterport3D dataset as training and testing data.  Please check [surface_normal](surface_normal) for train/test codes, pretrained models, and auxiliary data for the experiments. Note that to run the experiments for surface normal estimation, you don't need to download the whole dataset. Please see [here](surface_normal/readme.md) for how to download the data and run the scripts for this task.

## Semantic Voxel Labeling
The semantic voxel labeling task predicts per-voxel class labels for a scan. Please see [semantic_voxel_label](semantic_voxel_label).

## Region Classification  
The region classification task aims to predict the semantic category of the region (e.g., bedroom, kitchen, patio, etc.) containing the camera viewpoint of an RGB image or panorama.   It leverages semantic boundaries and labels for manually-specified regions in the Matterport3D dataset.  Please check [region_classification](region_classification) for how to download the data and run the script for the task (all train/test code, pretrained models, and auxiliary data for the experiments are provided).
