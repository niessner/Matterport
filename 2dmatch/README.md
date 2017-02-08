# 2DMatch for MP Data

## Requirements
- Install [Torch](http://torch.ch/docs/getting-started.html) on a machine with CUDA GPU
- Install [cuDNN v5](https://developer.nvidia.com/cudnn)

## Instructions

```bash
git clone https://github.com/andyzeng/2dmatch-mp.git 2dmatch-mp
cd 2dmatch-mp
wget https://d2j0dndfm35trm.cloudfront.net/resnet-101.t7
th train.lua 

th test.lua # test on a saved snapshot model
```

For visualizations, enable `qtwidget` calls and run `qlua -lenv train.lua`.