# Semantic Voxel Labeling

Predicts per-voxel class labels for a scan, which are labeled column-by-column. 
Code is adapted from the [semantic voxel labeling task of ScanNet](https://github.com/ScanNet/ScanNet/tree/master/Tasks/), but uses only occupancy information in a 31x31x62 neighborhood to labeling the 1x1x62 center column.

## Installation
Training uses [Torch7](http://torch.ch/docs/getting-started.html), with torch packages `cudnn`, `cunn`, `hdf5`, `xlua`.

## Training
```
th train.lua --train_data [path to train h5 file list] --test_data [path to test h5 file list] --save [output path]
```
with the appropriate paths to the train/test data (use '--help' to see more options)