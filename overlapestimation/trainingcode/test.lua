require 'cutorch';
require 'cunn';
require 'cudnn'
require 'image'
require 'util'
require 'hdf5'
require 'DistanceRatioCriterion'
require 'DataLoader'

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
   continue_train = 1,
   use_color = 0,
   numtest = 1600,
   snapshotname = 'lastest', 
   phase = 'test',
   listfilename = 'ovTestpairs',
   loseType = 'binary' --{'binary' or 'diff'}
}
for k,v in pairs(opt) do opt[k] = tonumber(os.getenv(k)) or os.getenv(k) or opt[k] end
if opt.display == 0 then opt.display = false end
if opt.use_color == 0 then opt.use_color = false end

        
if opt.loseType =='diff' then 
   opt.name = opt.name..opt.loseType
end
print(opt.name)


local dataset = DataLoader(opt)
local lines,numoflines,maxPathLength = lines_from('../data/'..opt.listfilename)
opt.numtest = numoflines

-- Load training snapshot model
local model = torch.load(paths.concat(opt.checkpoints_dir, opt.name, opt.snapshotname..".net"))
-- Extract one tower from the training model (3 cloned towers)
local patchEncoder = model:get(1):get(1):clone()
patchEncoder:evaluate()

local posRatioModel,distanceAll
if opt.loseType == 'diff' then 
   posRatioModel = model:get(2):get(3):clone()
   distanceAll = torch.Tensor(opt.numtest,2):fill(0)
else
   distanceAll = torch.Tensor(opt.numtest,1):fill(0)
end



for testId = 1,opt.numtest do 
	local idx1 = string.find(lines[testId],' ')
	local sceneId =  string.sub(lines[testId],1,idx1-1)
	local idx2 = string.find(lines[testId],' ', idx1+1)
	local I =  tonumber(string.sub(lines[testId],idx1+1,idx2-1))
	local idx3 = string.find(lines[testId],' ', idx2+1)
	local J =  tonumber(string.sub(lines[testId],idx2+1,idx3-1))
	-- Get example triplet patches
	local matchPatch,anchorPatch = dataset:getBatch({sceneId, I,J})

	-- Compute features per patch
	matchFeat = patchEncoder:forward(matchPatch:cuda()):clone()
	anchorFeat = patchEncoder:forward(anchorPatch:cuda()):clone()
	-- Descriptor distance between matching patches
	
	if opt.loseType == 'diff' then 
	   local overlapPred = posRatioModel:forward({matchFeat,anchorFeat})
	   distanceAll[testId][1] = anchorFeat:dist(matchFeat)
	   distanceAll[testId][2] = overlapPred:float()
	else
	   distanceAll[testId] = anchorFeat:dist(matchFeat);
	end
	
	-- image.save(paths.concat(opt.checkpoints_dir, opt.name, 'testPatch'.. testId .. distanceAll[testId] ..'.png'), image.toDisplayTensor(torch.cat(matchPatch,anchorPatch,2)))
	if testId%100 ==1å then 
	   print('Test: '.. testId..'total '..opt.numtest )
	end
end

local myFile = hdf5.open(paths.concat(opt.checkpoints_dir, opt.name, opt.listfilename..'matchresult.h5'),  'w')
myFile:write('distanceAll', distanceAll)
myFile:close()