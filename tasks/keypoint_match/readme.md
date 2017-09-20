# Keypoint Matching

Learning feature descriptors for matching keypoints to establish correspondences between image data..

## Installation
Training uses [Torch7](http://torch.ch/docs/getting-started.html), with torch packages `cudnn`, `cunn`.

## Training
Code is in [2dmatch](2dmatch).

```
th train.lua --basePath [path to keypointmatch data (images, matches files)] --train_data [list of train scene ids] --save [output path]
```
(use '--help' to see more options)
