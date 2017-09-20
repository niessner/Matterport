require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'
require 'optim'
require 'DataLoader'
-- Custom files
require 'model'
require 'util'
-- require 'qtwidget' -- for visualizing images

opt = {
   batchSize = 8,         -- number of samples to produce
   loadSize = 224,         -- resize the loaded image to loadsize maintaining aspect ratio. 0 means don't resize. -1 means scale randomly between [0.5,2] -- see donkey_folder.lua
   fineSize = 224,         -- size of random crops
   nc = 1,                 -- # of channels in input
   nThreads = 4,           -- #  of data loading threads to use
   niter = 25,             -- #  of iter at starting learning rate
   lr = 0.0002,            -- initial learning rate for adam
   beta1 = 0.5,            -- momentum term of adam
   ntrain = math.huge,     -- #  of examples per epoch. math.huge for full dataset
   display = 1,            -- display samples while training. 0 = false
   display_id = 10,        -- display window id.
   display_iter = 50,      -- # number of iterations after which display is updated
   gpu = 3,                -- gpu = 0 is CPU mode. gpu=X is GPU mode on GPU X
   name = 'train1',        -- name of the experiment you are running
   manualSeed = 0,         -- 0 means random seed
   dataRoot = '/n/fs/rgbd/data/matterport/v1/',
   checkpoints_dir = './checkpoints/',
   continue_train = 0,
   snapshotname = 'lastest', 
   dataset ='mp',
   net='',
   use_color = 0,
   save_latest_freq = 1000,
   phase = 'train',
   seqlist = '../data/trainval',
   loseType = 'binary' --{'binary' or 'diff'}
}


for k,v in pairs(opt) do opt[k] = tonumber(os.getenv(k)) or os.getenv(k) or opt[k] end
if opt.dataset == 'sun3d' then 
   opt.dataRoot = '/n/fs/rgbd/data/sun3d/'
   opt.seqlist = '../data/sun3d_train'
   opt.name = 'train_sun3d'
end
if opt.display == 0 then opt.display = false end
if opt.use_color == 0 then opt.use_color = false end
if opt.continue_train == 0 then opt.continue_train = false end


if opt.gpu>0 then 
   cutorch.setDevice(opt.gpu)
end

-- Set RNG seed
math.randomseed(os.time())

-- Load model and criterion
if opt.loseType =='diff' then 
   opt.name = opt.name..opt.loseType
end

local model
if opt.continue_train then
   print('loading previously trained model ...')
   model = torch.load(opt.net)
   opt.name = opt.name..'_ft'
else
	if opt.loseType =='binary' then
	   model = getModel()
	else
	   model = getModelRatio()
	end
end


local dataset = DataLoader(opt)
-- Defines triplet loss from "Deep Metric Learning using Triplet Network" http://arxiv.org/abs/1412.6622
local criterionTP = nn.DistanceRatioCriterion(true)
-- local criterionMSE = nn.MSECriterion()
local criterionAbs = nn.AbsCriterion()

local input_matchPatch = torch.zeros(opt.batchSize,3,opt.fineSize,opt.fineSize)
local input_anchorPatch   = torch.zeros(opt.batchSize,3,opt.fineSize,opt.fineSize)
local input_nonMatchPatch   = torch.zeros(opt.batchSize,3,opt.fineSize,opt.fineSize)
local matchRatio =  torch.zeros(opt.batchSize,1)
if opt.gpu >0 then 
	model = model:cuda()
	criterionTP = criterionTP:cuda()
	-- criterionMSE = criterionMSE:cuda()
	criterionAbs = criterionAbs:cuda()
	input_nonMatchPatch = input_nonMatchPatch:cuda()
	input_anchorPatch = input_anchorPatch:cuda()
	input_matchPatch  = input_matchPatch:cuda()
	matchRatio = matchRatio:cuda()
end

paths.mkdir(paths.concat(opt.checkpoints_dir, opt.name))
local params,gradParams = model:getParameters()
for trainIter = 1,15000 do

	local matchPatch,anchorPatch,nonMatchPatch,matchRatio_c = dataset:getBatch()
	input_matchPatch:copy(matchPatch)
	input_anchorPatch:copy(anchorPatch)
	input_nonMatchPatch:copy(nonMatchPatch)
	matchRatio:copy(matchRatio_c)

	-- Visualize training examples (console command: qlua -lenv train.lua)
	-- image.display{image=mosaic,offscreen=false,win=visWindow}

	local feval = function(x)

	    -- Update model parameters
	    if x ~= params then
	        params:copy(x) 
	    end
	    -- Reset gradients
	    gradParams:zero() 

	    -- Forward and backward pass
	    local input = {input_matchPatch,input_anchorPatch,input_nonMatchPatch}
	    local output = model:forward(input)
	    local loss,lossTP,lossR = 0,0,0
		local dlossR1,dLossTP
	    
	    lossTP = criterionTP:forward({output[1],output[2]})
	    dLossTP = criterionTP:backward({output[1],output[2]})
	    

	    if opt.loseType =='binary' then 
	       model:backward(input,dLossTP)
	       loss = lossTP
	    else
	    	lossR = criterionAbs:forward(output[3], matchRatio)
	    	dlossR = criterionAbs:backward(output[3],matchRatio)
	    	model:backward(input,{dLossTP[1],dLossTP[2],dlossR})
	    	loss = lossTP+lossR
	    end
	    
	    print('Training iteration '..trainIter..'lossTP: '..lossTP .. 'dlossR: ' .. lossR)
	    return loss,gradParams
	end

	-- Update model parameters (SGD)
	config = {learningRate = 0.001,momentum = 0.99}
	optim.sgd(feval,params,config)

	-- Save training snapshot of model
	if trainIter%50 == 1 then
	   image.save(paths.concat(opt.checkpoints_dir, opt.name, 'matchPatch.png'), image.toDisplayTensor(matchPatch))
	   image.save(paths.concat(opt.checkpoints_dir, opt.name, 'anchorPatch.png'), image.toDisplayTensor(anchorPatch))
	   image.save(paths.concat(opt.checkpoints_dir, opt.name, 'nonMatchPatch.png'), image.toDisplayTensor(nonMatchPatch))
    end

    if trainIter%1000 == 1 then
       local filename = paths.concat(opt.checkpoints_dir, opt.name, "lastest.net")
       print('saving:' .. filename)
       torch.save(filename, model:clearState()) 
    end

	if trainIter%1000 == 0 then
		local filename = paths.concat(opt.checkpoints_dir, opt.name, trainIter..".net")
		torch.save(filename, model:clearState()) 
	end
end
