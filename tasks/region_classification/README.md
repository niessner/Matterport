# Region Classification

In this task we predict the semantic category for each region within a house.

## Requirements
- Install [Torch](http://torch.ch/docs/getting-started.html) on a machine with CUDA GPU
- Install [cuDNN v5](https://developer.nvidia.com/cudnn)
- download pretrained resnet model
    ```
    wget https://d2j0dndfm35trm.cloudfront.net/resnet-50.t7
    ```


## pre-trained models
| Training                               | Model name            |
| -------------------------------------- | ----------------------|
| Single Image                           | snapshot-single.t7    |
| Panorama Image                         | snapshot-pano.t7      |

To download the models:
```
python download_mp.py  --task_data regoin_classification_data 
```

## Training 

```
imgList=[path training image list] labelList=[path training image list] gpu=1 snapshotsFolder=[output path] th train.lua
```

Example:
```
imgList=./data/train_room_pano_image.txt labelList=./data/train_room_pano_label.txt gpu=1 snapshotsFolder=./pano th train.lua
```

The snapshot models will be saved in [output path].

## Testing

```
imgList=[path training image list]  labelList=[path training image list] gpu=1 modelPath=[path to model] th test.lua
```
