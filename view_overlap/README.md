# Code for view overlap prediction
## Requirements
- Install [Torch](http://torch.ch/docs/getting-started.html) on a machine with CUDA GPU
- Install [cuDNN v5](https://developer.nvidia.com/cudnn)
- download pretrained resnet model
    ```
    cd trainingcode
    wget https://d2j0dndfm35trm.cloudfront.net/resnet-50.t7
    ```
- matterport data path: ```$path2matterportdata```
- training testing split:  ```../data/trainval```
## Training
```
cd trainingcode
dataRoot=$path2matterportdata loseType=diff th train.lua 
```
This program trains a view overlap prediction model with following options:
```loseType= diff ```  train with regrassion lose for overlap esitmation.
```loseType= binary ```  train without regrassion lose.
The snapshot models will be saved in
```trainingcode/checkpoints/<$name>/<number_of_iteration>.net```
## Testing
```
name=train1diff  th testFeat.lua
```
This program will load the trained model in folder ```$name```, and compute the image feature for all images in the test scene, the result is saved in ```checkpoints/<$name>/result/....h5```

## Evaluation
The following code will load the computed feature vectors in "testing" step and compute normalized discounted cumulative gain for each testing sequence. 

```
cd trainingcode
matlab &
mpOverlapEval('train1diff')
```