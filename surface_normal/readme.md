# 3D surface normal estimation from a single color image on Matterport Dataset

This is an extension of the 3D surface normal estimation code from: https://github.com/yindaz/surface_normal, with the support of training, testing, and evaluating on the matterport dataset.
	
Please contact me (yindaz AT cs DOT princeton DOT edu) if you have any problem.

## Testing
You need to specify the `-root_path` in [`config.lua`](./config.lua) to the root path of the matterport dataset.

We provide several pretrained models for testing. You can download them in [`./model/`](./model/).

Run
```
th main_test_matterport.lua -test_model ./model/sync_mp.t7
```
and the result will be in [`./result/`](./result/). The estimated normal will be saved in 16-bit PNG format, where 0-65535 in R,G,B channel correspond to [-1, 1] for the X,Y,Z component of the normal vector. We use the camera coordinates defined as - X points to the camera right, Y points to the camera forward, and Z points to the camera upward. For example, right facing wall are very red, floor are very blue, and you rarely see green as it's parallel to the camera viewing direction. 

## Training

Run
```
th main_train_matterport.lua -ps ./model/train_example
```
The training process will start, and snapshots will be saved as `./model/train_example_iter_xxx.t7`. To finetune a model, run
```
th main_train_matterport.lua -finetune -finetune_model ./model/train_example.t7 -ps ./model/train_finetune_example
```

## Evaluate
In MATLAB, run 
```
evaluate_normal_mp('./result/sync_mp_matterport_test/','Path to root of matterport dataset')
```
This function returns the pixelwise auglar loss on each testing image.

## Data
To train on synthetic image, you can find the training data from http://pbrs.cs.princeton.edu. Specifically,
- `Color image`: http://pbrs.cs.princeton.edu/pbrs_release/data/mlt_v2.zip (278GB)
- `Surface normal ground truth`: http://pbrs.cs.princeton.edu/pbrs_release/data/normal_v2.zip (27GB)
- `Data list`: http://pbrs.cs.princeton.edu/pbrs_release/data/data_goodlist_v2.txt

To experiment on NYUv2 data,
- `Color image and ground truth`: http://pbrs.cs.princeton.edu/pbrs_release/nyu/nyu_data.zip. This file is converted using data from http://www.cs.nyu.edu/~deigen/dnl/ and http://cs.nyu.edu/~silberman/datasets/nyu_depth_v2.html
- `Training data list`: http://pbrs.cs.princeton.edu/pbrs_release/nyu/trainNdxs.txt
- `Testing data list`: http://pbrs.cs.princeton.edu/pbrs_release/nyu/testNdxs.txt

