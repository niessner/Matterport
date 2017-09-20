require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'

-- Custom requires
require 'utils'
require 'DataLoader'

-- imgList=/n/fs/modelnet/matterportMultiView//data/test_room_pano_image.txt labelList=/n/fs/modelnet/matterportMultiView//data/test_room_pano_label.txt gpu=7 modelPath=./pano/snapshot-15000.t7 th test.lua 2>&1 | tee pano_result.txt
-- imgList=/n/fs/modelnet/matterportMultiView//data/test_room_single_image.txt labelList=/n/fs/modelnet/matterportMultiView//data/test_room_single_label.txt gpu=7 modelPath=./single/snapshot-15000.t7 th test.lua 2>&1 | tee single_result.txt


-- User options for model and data loading
options = {
  imgList = './data/train-imgs.txt',
  labelList = './data/train-labels.txt',
  batchSize = 8,
  doShuffle = false,
  doFlip = false,
  modelPath = './snapshots/snapshot-1000.t7',
  gpu = 1
}
for k,v in pairs(options) do options[k] = tonumber(os.getenv(k)) or os.getenv(k) or options[k] end

-- Set RNG seed
math.randomseed(os.time())

-- Set GPU device (0 - CPU)
cutorch.setDevice(options.gpu)

-- Create data sampler
print('Loading dataset...')
local dataLoader = DataLoader(options)

-- Load classification model
local model = torch.load(options.modelPath)
model:add(nn.LogSoftMax())
model:evaluate()
model = model:cuda()

-- Print image paths and classification results
for imgIdx = 1,dataLoader.numImgs do
    -- print('Classification ('..imgIdx..'/'..dataLoader.numImgs..'): '..dataLoader.imgPaths[imgIdx])
    local img = image.load(dataLoader.imgPaths[imgIdx])
    img = image.scale(img,224,224)
    local mean = {0.485,0.456,0.406}
    local std = {0.229,0.224,0.225}
    for c=1,3 do
        img[c]:add(-mean[c])
        img[c]:div(std[c])
    end
    local output = model:forward(img:cuda()):float()
    maxConf,predClass = torch.max(output,2)
    print(predClass[1][1]..','..dataLoader.imgLabels[imgIdx])
end

print('Finished.')



