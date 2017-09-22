require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'
require 'optim'

-- Custom requires
require 'utils'
require 'model'
require 'DataLoader'

-- imgList=/n/fs/modelnet/matterportMultiView//data/train_room_pano_image.txt labelList=/n/fs/modelnet/matterportMultiView//data/train_room_pano_label.txt gpu=8 snapshotsFolder=./pano th train.lua
-- imgList=/n/fs/modelnet/matterportMultiView//data/train_room_single_image.txt labelList=/n/fs/modelnet/matterportMultiView//data/train_room_single_label.txt gpu=8 snapshotsFolder=./single th train.lua

-- User options for data loading
options = {
  imgList = './data/train-imgs.txt',
  labelList = './data/train-labels.txt',
  batchSize = 8,
  doShuffle = true,
  doFlip = true,
  snapshotsFolder = './snapshots',
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

-- Build deep model
print('Building deep model...')
print('Number of training classes: '..dataLoader.numClasses)
local model,criterion = getModel(dataLoader.numClasses)
model:training()
model = model:cuda()
criterion = criterion:cuda()

-- Construct window for visualizing training examples
-- local visWindow = qtwidget.newwindow(1344,448)

-- Get model parameters
local params,gradParams = model:getParameters()

for trainIter = 1,200000 do

	-- Create a mini-batch of training examples
	local input,label = dataLoader:getMiniBatch()

	
	-- Convert input and label to GPU memory
	input = input:cuda()
	label = label:cuda()

	local feval = function(x)

	    -- Update model parameters
	    if x ~= params then
	        params:copy(x) 
	    end

	    -- Reset gradients
	    gradParams:zero() 

		local output = model:forward(input)

	    -- Compute classification loss (object class)
	    local loss = criterion:forward(output,label)
	    dloss = criterion:backward(output,label)

	    -- Backward pass
        model:backward(input,dloss)

	    if trainIter%10 == 0 then
	    	print('Training iteration '..trainIter..': '..loss)
	    end

	    return loss,gradParams
	end

	-- Update model parameters (SGD)
	local config = {learningRate = 0.001,momentum = 0.99}
	optim.sgd(feval,params,config)

	-- Save training snapshot of model
	if trainIter%5000 == 0 then
		local filename = paths.concat(options.snapshotsFolder,'snapshot-'..trainIter..'.t7')
		os.execute('mkdir -p '..sys.dirname(filename))
		torch.save(filename, model:clearState()) 
	end
end

       
