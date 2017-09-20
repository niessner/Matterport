require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'
local threads = require 'threads'
threads.Threads.serialization('threads.sharedserialize')

local DataLoader = torch.class('DataLoader')

-- Initialize multi-thread data loader
function DataLoader:__init(options)

    -- Get user options
    self.imgList = options.imgList
    self.labelList = options.labelList
    self.batchSize = options.batchSize
    self.doShuffle = options.doShuffle
    self.doFlip = options.doFlip

    -- Read paths for training images and labels
    self.imgPaths = getLinesFromFile(self.imgList)
    self.imgLabels = getLinesFromFile(self.labelList)
    self.numImgs = #self.imgPaths

    -- Initialize shuffle indices and count number of classes
    self.shuffleIdx = {}
    self.numClasses = 0;
    for imgIdx = 1,self.numImgs do
        self.imgLabels[imgIdx] = tonumber(self.imgLabels[imgIdx])
        self.shuffleIdx[imgIdx] = imgIdx
        self.numClasses = math.max(self.numClasses,self.imgLabels[imgIdx])
    end

    -- Shuffle training data
    self.shuffleIdx = shuffleTable(self.shuffleIdx,self.numImgs)
    self.trainIdx = 1
    self.trainEpochIdx = 1
    self.trainEpochSize = self.numImgs

    -- Define multi-thread pool
    self.nthread = self.batchSize
    self.pool = threads.Threads(self.nthread,
        function(threadid)
            pcall(require, 'image')
            math.randomseed(threadid)
        end)
end

-- Load mini-batch of data with multi-thread
function DataLoader:getMiniBatch()
    local input = torch.zeros(self.batchSize,3,224,224)
    local label = torch.zeros(self.batchSize,1)

    -- Get paths and labels of images to be loaded
    local batchImgPaths = {}
    for sampleIdx = 1,self.batchSize do
        batchImgPaths[sampleIdx] = self.imgPaths[self.shuffleIdx[self.trainIdx]]
        label[sampleIdx] = self.imgLabels[self.shuffleIdx[self.trainIdx]]

        -- Re-shuffle data if at end of training epoch
        if self.trainIdx == self.trainEpochSize then
            self.shuffleIdx = shuffleTable(self.shuffleIdx,self.numImgs)
            self.trainIdx = 1
            self.trainEpochIdx = self.trainEpochIdx+1
        else
            self.trainIdx = self.trainIdx+1
        end
    end

    -- Load images with multi-thread
    local doFlip = self.doFlip
    for jobIdx = 1,self.nthread do
        self.pool:addjob(
            function()

                -- Load and pre-process image
                local img = image.load(batchImgPaths[jobIdx])
                if doFlip and torch.uniform()>0.7 then
                   img = image.hflip(img)
                end
                img = image.scale(img,224,224)
                local mean = {0.485,0.456,0.406}
                local std = {0.229,0.224,0.225}
                for c=1,3 do
                    img[c]:add(-mean[c])
                    img[c]:div(std[c])
                end
                input[jobIdx] = img:reshape(1,3,224,224);

                return __threadid
            end)
    end

    self.pool:synchronize()
    collectgarbage()

    return input,label
end
