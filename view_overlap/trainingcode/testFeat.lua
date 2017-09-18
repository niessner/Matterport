require 'cutorch';
require 'cunn';
require 'cudnn'
require 'image'
require 'util'
require 'hdf5'
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
   name = 'train1diff',        -- name of the experiment you are running
   manualSeed = 0,         -- 0 means random seed
   dataRoot = '/n/fs/rgbd/data/matterport/v1/',
   checkpoints_dir = './checkpoints/',
   continue_train = 1,
   use_color = 0,
   snapshotname = 'lastest', 
   phase = 'test',
   sceneId = '',
   dataset = 'mp',
   loseType = 'binary' --{'binary' or 'diff'}
}
for k,v in pairs(opt) do opt[k] = tonumber(os.getenv(k)) or os.getenv(k) or opt[k] end
if opt.display == 0 then opt.display = false end
if opt.use_color == 0 then opt.use_color = false end

if opt.dataset == 'sun3d' then 
   opt.dataRoot = '/n/fs/rgbd/data/sun3d/'
   opt.seqlist = '../data/sun3d_test'
else
   opt.seqlist = '../data/test'
end     
print(opt.name)


local dataset = DataLoader(opt)
-- local lines,numoflines,maxPathLength = lines_from('../data/'..opt.listfilename)


-- Load training snapshot model
local model = torch.load(paths.concat(opt.checkpoints_dir, opt.name, opt.snapshotname..".net"))
-- Extract one tower from the training model (3 cloned towers)
local patchEncoder = model:get(1):get(1):clone()
patchEncoder:evaluate()


for i = 1,dataset.numofScene do 
   opt.sceneId = dataset.sceneIds[i]
   print(opt.sceneId)
   opt.numtest = dataset:SceneSize(opt.sceneId)
   local featAll = torch.Tensor(opt.numtest,512):fill(0)

   for testId = 1,opt.numtest do 
   	
   	-- Get example triplet patches
   	local matchPatch,anchorPatch = dataset:getBatch({opt.sceneId, testId,testId})

   	-- Compute features per patch
   	local matchFeat = patchEncoder:forward(matchPatch:cuda()):clone()
   	-- Descriptor distance between matching patches
   	featAll[testId]:copy(matchFeat)
   	if testId%100 ==1 then 
   	   print('Test: '.. testId..'total '..opt.numtest)
   	end
   end
   local filename = paths.concat(opt.checkpoints_dir, opt.name,'result', string.gsub(opt.sceneId,'/','--')..'featAll.h5')
   paths.mkdir(paths.concat(opt.checkpoints_dir, opt.name,'result'))
   print('write to: '..filename)
   local myFile = hdf5.open(filename,  'w')
   myFile:write('featAll', featAll)
   myFile:close()
end