require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'
require 'optim'

-- Custom files
require 'model'
require 'sys'
-- require 'qtwidget' -- for visualizing images
dofile('util.lua')

opt_string = [[
	-h,--help									   	print help
	-s,--save			(default "logs")			subdirectory to save logs
	-b,--batchSize		(default 32)				batch size
	-g,--gpu_index		(default 0)			 		GPU index (start from 0)
	--max_epoch			(default 10)				maximum number of epochs
	--basePath			(default "/mnt/raid/datasets/Matterport/Matching1/")	base path for train/test data
	--train_data		(default "scenes_train_small.txt")	 						txt file containing train
	--test_data			(default "scenes_test_small.txt")	  							txt file containing test
	--patchSize			(default 64)				patch size to extract (resized to 224)
	--matchFileSkip		(default 80)				only use every skip^th keypoint match in file
	--use_bottleneck	(default true)
	--learn_metric		(default false)
]]

opt = lapp(opt_string)

-- print help or chosen options
if opt.help == true then
	print('Usage: th train.lua')
	print('Options:')
	print(opt_string)
	os.exit()
else
	print(opt)
end

-- set gpu
cutorch.setDevice(opt.gpu_index+1)

-- Set RNG seed
--math.randomseed(os.time())
math.randomseed(0)
torch.manualSeed(0)


-- Load model and criterion
model,criterion = getModel(opt.use_bottleneck, opt.learn_metric)
model = model:cuda()
critrerion = criterion:cuda()
parameters, gradParameters = model:getParameters()
--print(model)

-- Construct window for visualizing training examples
-- visWindow = qtwidget.newwindow(1792,672)


-- load training and testing files
train_files = getDataFiles(paths.concat(opt.basePath,opt.train_data), opt.basePath) --filter out non-existent scenes
test_files = getDataFiles(paths.concat(opt.basePath,opt.test_data), opt.basePath)   --filter out non-existent scenes
print('#train files = ' .. #train_files)
print('#test files = ' .. #test_files)


-- config logging
testLogger = optim.Logger(paths.concat(opt.save, 'test.log'))
testLogger:setNames{'epoch', 'iteration', 'current train loss', 'avg train loss'}
do
	local optfile = assert(io.open(paths.concat(opt.save, 'options.txt'), 'w'))
	local cur = io.output()
	io.output(optfile)
	serialize(opt)
	serialize(train_files)
	serialize(test_files)
	io.output(cur)
	optfile:close()
end

local patchSize = opt.patchSize
local saveInterval = 1000

------------------------------------
-- Training routine
--
function train()
	model:training()
	epoch = epoch or 1 -- if epoch not defined, assign it as 1
	print('epoch ' .. epoch)
	--if epoch % opt.epoch_step == 0 then optimState.learningRate = optimState.learningRate/2 end

	--load in the train data (positive and negative matches) 
	local poss, negs = loadMatchFiles(opt.basePath, train_files, patchSize/2, opt.matchFileSkip)
	print(#poss)
	--print(poss)
	--print(negs)
	collectgarbage()

	--pre-allocate memory
	local inputs0 = torch.zeros(opt.batchSize, 1, 64, 64):cuda()
	local inputs1 = torch.zeros(opt.batchSize, 1, 64, 64):cuda()
	local targets = torch.ones(opt.batchSize):cuda()
	targets[{{opt.batchSize/2+1,opt.batchSize}}]:fill(-1) -- first half is positive, second half is negative

	local totalloss = 0	
	local indicesPos = torch.randperm(#poss)
	local indicesNeg = torch.randperm(#negs)
	local numSamples = 2*#poss
	local numIters = math.floor(numSamples/opt.batchSize)

	for iter = 1,numSamples,opt.batchSize do
		-- print progress bar :D		
		local trainIter = (iter-1)/opt.batchSize+1
		xlua.progress(trainIter, numIters)
		if iter + opt.batchSize > numSamples then break end --don't use last batch 
		collectgarbage()
		local offset = (iter - 1)/2
		for k = 1,opt.batchSize/2 do --create a mini batch
			local idxPos = indicesPos[offset+k]
			local idxNeg = indicesNeg[offset+k]
			local scenePos = poss[idxPos][1]
			local imgPathPos = paths.concat(opt.basePath,scenePos,'images')
			local sceneNeg = negs[idxNeg][1]
			local imgPathNeg = paths.concat(opt.basePath,sceneNeg,'images')
			
			local p0,p1 = getTrainingExample(imgPathPos, poss[idxPos][2], poss[idxPos][3], patchSize)
			local n0,n1 = getTrainingExample(imgPathNeg, negs[idxNeg][2], negs[idxNeg][3], patchSize)
			inputs0[{k,1,{},{}}]:copy(p0)
			inputs1[{k,1,{},{}}]:copy(p1)
			inputs0[{opt.batchSize/2+k,{},{},{}}]:copy(n0)
			inputs1[{opt.batchSize/2+k,{},{},{}}]:copy(n1)
			
			--[[print(k .. '\t' .. opt.batchSize/2+k)
			image.save('p0.png', p0)
			image.save('p1.png', p1)
			image.save('n0.png', n0)
			image.save('n1.png', n1)
			io.read()--]]
		end

		-- a function that takes single input and return f(x) and df/dx
		local curLoss = -1
		local feval = function(x)
			if x ~= parameters then parameters:copy(x) end
			gradParameters:zero()
			-- Forward and backward pass
			local inputs = {inputs0, inputs1}
			local output = model:forward(inputs)
			local loss = criterion:forward(output, targets)
			--print('Training iteration '..trainIter..': '..loss)
			local dLoss = criterion:backward(output, targets)
			model:backward(inputs,dLoss)
			curLoss = loss
			totalloss = totalloss + loss
			return loss,gradParameters
		end

		config = {learningRate = 0.001,momentum = 0.99}
		optim.sgd(feval, parameters, config)

		if testLogger then
			paths.mkdir(opt.save)
			testLogger:add{tostring(epoch), tostring(trainIter), curLoss, totalloss / trainIter}
			testLogger:style{'-','-','-','-'}
		end
  

		if trainIter > 0 and (trainIter % saveInterval == 0 or trainIter == numIters) then
			local filename = paths.concat(opt.save, 'model_' .. tostring(epoch) .. '-' .. tostring(trainIter) .. '.net')
			print('==> saving model to '..filename)
			torch.save(filename, model)--:clearState())
		end	   
	end
	
	epoch = epoch + 1
end

-----------------------------------------
-- Start training
--
for i = 1,opt.max_epoch do
	train()
end

