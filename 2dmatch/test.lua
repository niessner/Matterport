require 'cutorch';
require 'cunn';
require 'cudnn'
require 'image'
require 'util'
-- require 'qtwidget'

local basePath = '/mnt/raid/datasets/Matterport/Matching/'
local patchSize = 64

opt_string = [[
    -h,--help                                      			 print help
    --gpu_index     		(default 0)                      GPU index
	--model					(default "model.net")		     torch model file path
    --train_data            (default "scenes_train.txt")     txt file containing train
    --test_data             (default "scenes_test.txt")      txt file containing test
	--output_dir			(default "output")			 	 output dir for match distances
]]

opt = lapp(opt_string)

-- print help or chosen options
if opt.help == true then
    print('Usage: th test.lua')
    print('Options:')
    print(opt_string)
    os.exit()
else
    print(opt)
end

-- Set RNG seed
--math.randomseed(os.time())
math.randomseed(0)

-- set gpu
cutorch.setDevice(opt.gpu_index+1)

print('mem1: ')
print(cutorch.getMemoryUsage(1))

-- Load training snapshot model
assert(paths.filep(opt.model))
local patchEncoder
do 
	print('loading model: ' .. opt.model .. '...')
	local model = torch.load(opt.model)
	print('mem2: ')
	print(cutorch.getMemoryUsage(1))
	model:clearState()
	collectgarbage()
	model:float() --model uses crazy amount of memory...
	--print(model)
	
	print('mem3: ')
	print(cutorch.getMemoryUsage(1))
	
	-- Extract one tower from the training model (3 cloned towers)
	patchEncoder = model:get(1):get(1):clone()
	print('done')
	
	print('mem4: ')
	print(cutorch.getMemoryUsage(1))
end
collectgarbage()

patchEncoder:cuda()
--print(patchEncoder)

print('mem5: ')
print(cutorch.getMemoryUsage(1))


-- load training and testing files
local train_files = {}
if opt.train_data ~= '' and paths.filep(paths.concat(basePath,opt.train_data)) then
	train_files = getDataFiles(paths.concat(basePath,opt.train_data))
	print(train_files)
end
local test_files = getDataFiles(paths.concat(basePath,opt.test_data))
print(test_files)

function test(data_files)
	-- Set to testing mode (disable dropout and batch normalization)
	patchEncoder:evaluate()
	
	--load in the data (positive and negative matches) 
	local poss, negs = loadMatchFiles(basePath, data_files, patchSize/2)
	
	local inputs = {}
	
	for k = 1,#poss do
		-- print progress bar :D		
		xlua.progress(k, #poss)
		
		local sceneName = poss[k][1]
		local imgPath = paths.concat(basePath,sceneName,'images')
		
		-- Get example triplet patches
		local anc,pos,neg = getTrainingExampleTriplet(imgPath, poss[k][2], poss[k][3], negs[k][3], patchSize)
		
		-- Compute features per patch
		local ancFeat = patchEncoder:forward(anc:resize(1,3,224,224):cuda()):clone()
		local posFeat = patchEncoder:forward(pos:resize(1,3,224,224):cuda()):clone()
		local negFeat = patchEncoder:forward(neg:resize(1,3,224,224):cuda()):clone()
		
		-- descriptor distance between patches
		print('Distance between matching patches: ' ..       ancFeat:dist(posFeat))
		print('Distance between non-matching patches 1: ' .. posFeat:dist(negFeat))
		print('Distance between non-matching patches 2: ' .. ancFeat:dist(negFeat))
		
		io.read() --TODO save this to a file
	end
end

-- Test!
if #train_files > 0 then test(train_files) end
test(test_files)

