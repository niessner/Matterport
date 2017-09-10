# 3D surface normal estimation from a single color image on Matterport Dataset

This is an extension of the 3D surface normal estimation code from: https://github.com/yindaz/surface_normal, with the support of training, testing, and evaluating on the matterport dataset.
	
Please contact me (yindaz AT cs DOT princeton DOT edu) if you have any problem.

## Data
To run surface normal experiment, you need to download the undistorted color images and surface normal images. If you have signed the agreement and have the download script, you may call:
```
python download_mp.py --type undistorted_color_images undistorted_normal_images --task_data surface_normal_data [surface_normal_pretrain]
```
This gives you necessary training and testing image together with the ground truth. It also downloads a zip of files containing the list of training/testing data. You should unzip the list files here. Optionally, we provide several pretrained models for quick testing, which can be downloaded with `surface_normal_pretrain`. The pretrained models are trained as following:

| Training                               | Model name     |
| -------------------------------------- | -------------- |
| Synthetic data                         | sync.t7        |
| Synthetic data + Matterport3D          | sync_mp.t7     |
| Synthetic data + Matterport3D + NYUv2  | sync_mp_nyu.t7 |

## Testing
You need to specify the `-root_path` in [`config.lua`](./config.lua) to the root path of the matterport dataset. 

Assuming that pretrained models are saved in [`model`](./model), run
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

## External Data
To train on synthetic image, you can find the training data from http://pbrs.cs.princeton.edu. Specifically,
- `Color image`: http://pbrs.cs.princeton.edu/pbrs_release/data/mlt_v2.zip (278GB)
- `Surface normal ground truth`: http://pbrs.cs.princeton.edu/pbrs_release/data/normal_v2.zip (27GB)
- `Data list`: http://pbrs.cs.princeton.edu/pbrs_release/data/data_goodlist_v2.txt

To experiment on NYUv2 data,
- `Color image and ground truth`: http://pbrs.cs.princeton.edu/pbrs_release/nyu/nyu_data.zip. This file is converted using data from http://www.cs.nyu.edu/~deigen/dnl/ and http://cs.nyu.edu/~silberman/datasets/nyu_depth_v2.html
- `Training data list`: http://pbrs.cs.princeton.edu/pbrs_release/nyu/trainNdxs.txt
- `Testing data list`: http://pbrs.cs.princeton.edu/pbrs_release/nyu/testNdxs.txt

