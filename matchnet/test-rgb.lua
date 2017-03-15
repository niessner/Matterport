require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'

-- Custom files
require 'model-rgb'
require 'sys'
-- require 'qtwidget' -- for visualizing images
dofile('util-rgb.lua')

opt_string = [[
	-h,--help									   	print help
	-b,--batchSize		(default 32)				batch size
	-g,--gpu_index		(default 0)			 		GPU index (start from 0)
	--model				(default "")				model
	--basePath			(default "/mnt/raid/datasets/Matterport/Matching1/")	base path for train/test data
	--train_data		(default "scenes_train_small.txt")	 						txt file containing train
	--test_data			(default "scenes_test_small.txt")	  							txt file containing test
	--patchSize			(default 64)				patch size to extract (resized to 224)
	--matchFileSkip		(default 80)				only use every skip^th keypoint match in file
	--use_bottleneck	(default true)
	--learn_metric		(default false)
    --output_dir        (default "output")
	--imWidth           (default 640)           image dimensions in data folder
	--imHeight          (default 512)           image dimensions in data folder
	--detectImWidth     (default 1280)          image dimensions for key detection
	--detectImHeight    (default 1024)          image dimensions for key detection
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


-- Load model
local model
if paths.filep(opt.model) then 
	print('loading model: ' .. opt.model)
	model = torch.load(opt.model)
else
        print('using untrained model: ')
	model = getModel(opt.use_bottleneck, opt.learn_metric)
	model = model:cuda()
end
local criterion = nn.HingeEmbeddingCriterion(1)

-- load training and testing files
local train_files = {}
if paths.filep(paths.concat(opt.basePath,opt.train_data)) then
    train_files = getDataFiles(paths.concat(opt.basePath,opt.train_data), opt.basePath) --filter out non-existent scenes
end
local test_files = getDataFiles(paths.concat(opt.basePath,opt.test_data), opt.basePath)   --filter out non-existent scenes
print('#train files = ' .. #train_files)
print('#test files = ' .. #test_files)
local scaleX = opt.imWidth / opt.detectImWidth
local scaleY = opt.imHeight / opt.detectImHeight

------------------------------------
-- Training routine
--
function test(data_files, outpath)
	model:evaluate()
	
	--load in the train data (positive and negative matches) 
	local poss, negs = loadMatchFiles(opt.basePath, data_files, opt.patchSize/2, opt.matchFileSkip, opt.imWidth, opt.imHeight, scaleX, scaleY)
	print(#poss)
	print(#negs)
	collectgarbage()

	--pre-allocate memory
	local inputs0 = torch.zeros(1, 3, 64, 64):cuda()
	local inputs1 = torch.zeros(1, 3, 64, 64):cuda()
	local tgtPos = torch.ones(1)
	local tgtNeg = -torch.ones(1)
	
        -- output files
        local splitter = ',' --output as csv
        local outfile_anc = assert(io.open(paths.concat(outpath, 'anc_feat.txt'), "w"))
        local outfile_pos = assert(io.open(paths.concat(outpath, 'pos_feat.txt'), "w"))
        local outfile_neg = assert(io.open(paths.concat(outpath, 'neg_feat.txt'), "w"))
	
	local loss = 0
        local count = 0
        local avgPosDist = 0
        local avgNegDist = 0
        local num = #poss
	for k = 1,num do
		-- print progress bar :D		
		xlua.progress(k, num)
		collectgarbage()
		local scenePos = poss[k][1]
		local imgPathPos = paths.concat(opt.basePath,scenePos,'images')
		local sceneNeg = negs[k][1]
		local imgPathNeg = paths.concat(opt.basePath,sceneNeg,'images')
		
		local p0,p1 = getTrainingExample(imgPathPos, poss[k][2], poss[k][3], opt.patchSize)
		local n0,n1 = getTrainingExample(imgPathNeg, negs[k][2], negs[k][3], opt.patchSize)
		inputs0[{1,{},{},{}}]:copy(p0)
		inputs1[{1,{},{},{}}]:copy(p1)
		
		local dpos = model:forward({inputs0, inputs1}):clone()
                --print(output)
		local p0_feat = model:get(1):get(1).output:float()
		local p1_feat = model:get(1):get(2).output:float()
		local f = criterion:forward(dpos, tgtPos)
                --print(f)
                avgPosDist = avgPosDist + dpos[1]
		loss = loss + f
		
		inputs0[{1,{},{},{}}]:copy(n0)
		inputs1[{1,{},{},{}}]:copy(n1)
		local dneg = model:forward({inputs0, inputs1}):clone()
                --print(output)
		--local n0_feat = model:get(1):get(1).output:float()
		local n1_feat = model:get(1):get(2).output:float()
		f = criterion:forward(dneg, tgtNeg)
                --print(f)
                avgNegDist = avgNegDist + dneg[1]
		loss = loss + f
                --io.read()

                if dpos[1] < dneg[1] then count = count + 1 end
		for i=1,p0_feat:size(1) do
                    outfile_anc:write(string.format("%.6f", p0_feat[i]))
                    outfile_pos:write(string.format("%.6f", p1_feat[i]))
                    outfile_neg:write(string.format("%.6f", n1_feat[i]))
                    if i < p0_feat:size(1) then 
                        outfile_anc:write(splitter)
                        outfile_pos:write(splitter)
                        outfile_neg:write(splitter)
                    end
                end
                outfile_anc:write('\n')
                outfile_pos:write('\n')
                outfile_neg:write('\n')
        end
	print('loss: ' .. 0.5*loss/num)
        print('count/num = ' .. count/num .. ' (' .. count .. '/' .. num ..')')
        print('avg pos dist = ' .. avgPosDist/num)
        print('avg neg dist = ' .. avgNegDist/num)
end

if #train_files > 0 then
    local outpath = paths.concat(opt.output_dir, 'train')
    if not paths.dirp(outpath) then paths.mkdir(outpath) end
    test(train_files, outpath)
end
do 
    local outpath = paths.concat(opt.output_dir, 'test')
    if not paths.dirp(outpath) then paths.mkdir(outpath) end
    test(test_files, outpath)
end

