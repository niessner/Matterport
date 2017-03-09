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
    -h,--help                                       print help
    -s,--save               (default "logs")        subdirectory to save logs
    -b,--batchSize          (default 8)             batch size
    -g,--gpu_index          (default 0)             GPU index (start from 0)
    --max_epoch             (default 10)	        maximum number of epochs
    --basePath              (default "/mnt/raid/datasets/Matterport/Matching1/")        base path for train/test data
    --train_data            (default "scenes_trainval.txt")     txt file containing train
    --test_data             (default "scenes_test.txt")      txt file containing test
    --patchSize             (default 64)            patch size to extract (resized to 224)
    --matchFileSkip         (default 10)            only use every skip^th keypoint match in file
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
model,criterion = getModel()
model = model:cuda()
critrerion = criterion:cuda()
model:zeroGradParameters()
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
    local inputs_anc = torch.zeros(opt.batchSize, 3, 224, 224):cuda()
    local inputs_pos = torch.zeros(opt.batchSize, 3, 224, 224):cuda()
    local inputs_neg = torch.zeros(opt.batchSize, 3, 224, 224):cuda()

    local totalloss = 0	
    local indices = torch.randperm(#poss)
    local numIters = math.floor(#poss/opt.batchSize)

    for iter = 1,#poss,opt.batchSize do
        -- print progress bar :D		
        local trainIter = (iter-1)/opt.batchSize+1
        xlua.progress(trainIter, numIters)
        if iter + opt.batchSize > #poss then break end --don't use last batch 
        collectgarbage()
        for k = iter,iter+opt.batchSize-1 do --create a mini batch
            local idx = indices[k]
            local sceneName = poss[idx][1]
            local imgPath = paths.concat(opt.basePath,sceneName,'images')
			
            local anc,pos,neg = getTrainingExampleTriplet(imgPath, poss[idx][2], poss[idx][3], negs[idx][3], patchSize)
	    --pos = torch.add(anc, torch.rand(3, 224, 224):float() * 0.1)
            --neg = torch.rand(3, 224, 224):float()
            inputs_anc[{k-iter+1,{},{},{}}]:copy(pos) --match
            inputs_pos[{k-iter+1,{},{},{}}]:copy(anc) --anchor
            inputs_neg[{k-iter+1,{},{},{}}]:copy(neg) --non-match
        end

        -- a function that takes single input and return f(x) and df/dx
        local curLoss = -1
        local feval = function(x)
            if x ~= parameters then parameters:copy(x) end
            gradParameters:zero()
            -- Forward and backward pass
            local inputs = {inputs_pos, inputs_anc, inputs_neg}
            local output = model:forward(inputs)
            local loss = criterion:forward(output)
            --print('Training iteration '..trainIter..': '..loss)
            local dLoss = criterion:backward(output)
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
  

        if trainIter > 0 and (trainIter % saveInterval == 0 or trainIter == #indices) then
            local filename = paths.concat(opt.save, 'model_' .. tostring(epoch) .. '-' .. tostring(trainIter) .. '.net')
            print('==> saving model to '..filename)
            torch.save(filename, model)--:clearState())
        end	   
    end
	
    epoch = epoch + 1
end

-------------------------------------
-- Test routine
--

function test()
    --[[model:evaluate()
    for fn = 1, #test_files do
        local pos,anc,neg = loadDataFile(train_files[train_file_indices[fn])

        local filesize = (#current_data)[1]
        local indices = torch.randperm(filesize):long():split(opt.batchSize)
        for t, v in ipairs(indices) do
            local inputs = { pos:index(1,v):cuda(), anc:index(1,v):cuda(), neg:index(1,v):cuda() }       
            local outputs = model:forward(inputs)
            local loss = criterion:forward(output)
        end
    end
    --TODO PRINT TEST LOSS/ACCURACY 
    if testLogger then
        paths.mkdir(opt.save)
        testLogger:add{train_acc, confusion.totalValid * 100}
        testLogger:style{'-','-'}
    end
    --]]

    -- save model every 10 epochs
    if epoch % 10 == 0 then
      local filename = paths.concat(opt.save, 'model_' ..tostring(epoch) .. '.net')
      print('==> saving model to '..filename)
      torch.save(filename, model:clearState())
    end 
end

-----------------------------------------
-- Start training
--
for i = 1,opt.max_epoch do
    train()
    --test()
end

