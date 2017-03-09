require 'cutorch';
require 'cunn';
require 'cudnn'
require 'image'
require 'model'
dofile('util.lua')
-- require 'qtwidget'

-- example usage:  th test.lua --gpu_index 0 --model ~/data/matthias/Matterport/2dmatch/logs/model_20000.net --train_data none test_data scenes_test.txt

local basePath = '/mnt/raid/datasets/Matterport/Matching1/'

opt_string = [[
    -h,--help                                                print help
    --gpu_index             (default 0)                      GPU index
    --batchSize             (default 8)
    --model                 (default "model.net")            torch model file path
    --test_data             (default "scenes_test.txt")      txt file containing test
    --train_data            (default "")                     txt file containing train (optional)
    --output_dir            (default "output")               output dir for match distances
    --patchSize             (default 64)            patch size to extract (resized to 224)
    --matchFileSkip         (default 10)            only use every skip^th keypoint match in file
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
torch.manualSeed(0)

if (paths.dirp(opt.output_dir)) then
    print(sys.COLORS.red .. 'warning: output dir ' .. opt.output_dir .. ' already exists, press key to continue' .. sys.COLORS.none)
    io.read()
end

-- set gpu
cutorch.setDevice(opt.gpu_index+1)

-- Load training snapshot model
local model
if paths.filep(opt.model) then
    model = torch.load(opt.model)
else
    model = getModel()
    model:cuda()
end
-- Extract one tower from the training model (3 cloned towers)
--local patchEncoder = model:get(1):get(1):clone()

local test_files = getDataFiles(paths.concat(basePath,opt.test_data), basePath)
print(test_files)

function test(data_files, outpath)
    assert(paths.dirp(outpath))
    -- Set to testing mode (disable dropout and batch normalization)
    --patchEncoder:evaluate()
    
    --load in the data (positive and negative matches) 
    local poss, negs = loadMatchFiles(basePath, data_files, opt.patchSize/2, opt.matchFileSkip)
    print(#poss)
    print(#negs)
    local count = 0
    local err = 0
    local avgPosDist = 0
    local avgNegDist = 0
    
    -- output files
    splitter = ',' --output as csv
    local outfile_anc = assert(io.open(paths.concat(outpath, 'anc_feat.txt'), "w"))
    local outfile_pos = assert(io.open(paths.concat(outpath, 'pos_feat.txt'), "w"))
    local outfile_neg = assert(io.open(paths.concat(outpath, 'neg_feat.txt'), "w"))
    
    --pre-allocate memory
    local inputs_anc = torch.zeros(opt.batchSize, 3, 224, 224):cuda()
    local inputs_pos = torch.zeros(opt.batchSize, 3, 224, 224):cuda()
    local inputs_neg = torch.zeros(opt.batchSize, 3, 224, 224):cuda()

    local indices = torch.randperm(#poss)
    local numIters = math.ceil(#poss/opt.batchSize)
    for iter = 1,#poss,opt.batchSize do
        local trainIter = (iter-1)/opt.batchSize+1
        xlua.progress(trainIter, numIters)
        local curBatchSize = math.min(opt.batchSize, #poss - iter)
        for k = iter,iter+curBatchSize-1 do
            local idx = indices[k]
            --print('idx = ' .. idx)
            local sceneName = poss[idx][1]
            local imgPath = paths.concat(basePath,sceneName,'images')
            local anc,pos,neg = getTrainingExampleTriplet(imgPath, poss[idx][2], poss[idx][3], negs[idx][3], opt.patchSize)
            inputs_pos[{k-iter+1,{},{},{}}]:copy(pos)
            inputs_anc[{k-iter+1,{},{},{}}]:copy(anc)
            inputs_neg[{k-iter+1,{},{},{}}]:copy(neg)
        end
        
        local output = model:forward({inputs_pos, inputs_anc, inputs_neg})
        local ancft = model:get(1):get(2).output:float()
        local posft = model:get(1):get(1).output:float()
        local negft = model:get(1):get(3).output:float()
 
        -- descriptor distance between patches
        local dists_pos = torch.zeros(curBatchSize)
        local dists_neg = torch.zeros(curBatchSize)
        for k = 1,curBatchSize do
            dists_pos[k] = torch.dist(ancft[k], posft[k])
            dists_neg[k] = torch.dist(ancft[k], negft[k])
        end
        local loss = -torch.log( torch.cdiv( torch.exp(-dists_pos), torch.exp(-dists_pos)+torch.exp(-dists_neg) ) )
        err = err + torch.sum(loss)
        avgPosDist = avgPosDist + torch.sum(dists_pos)
        avgNegDist = avgNegDist + torch.sum(dists_neg)

        --[[print('dist anc pos = ' .. torch.sum(dists_pos)/opt.batchSize); print(dists_pos)
        print('dist anc neg = ' .. torch.sum(dists_neg)/opt.batchSize); print(dists_neg)
        print('loss = ' .. torch.sum(loss)/opt.batchSize); print(loss)
        print(output[1])
        print(output[2])
        io.read()--]]

        for k = 1,curBatchSize do
            if dists_pos[k] < dists_neg[k] then count = count + 1 end
            for i=1,ancft:size(2) do
                outfile_anc:write(string.format("%.6f", ancft[k][i]))
                outfile_pos:write(string.format("%.6f", posft[k][i]))
                outfile_neg:write(string.format("%.6f", negft[k][i]))
                if i < ancft:size(2) then 
                    outfile_anc:write(splitter)
                    outfile_pos:write(splitter)
                    outfile_neg:write(splitter)
                end
            end
            outfile_anc:write('\n')
            outfile_pos:write('\n')
            outfile_neg:write('\n')
        end
    end
    print('count = ' .. count .. ' of ' .. #poss)
    print(count/#poss)
    print('loss = ' .. err .. ' / ' .. #poss .. ' = ' .. err/#poss)
    print('avg pos dist = ' .. avgPosDist .. ' / ' .. #poss .. ' = ' .. avgPosDist/#poss)
    print('avg neg dist = ' .. avgNegDist .. ' / ' .. #poss .. ' = ' .. avgNegDist/#poss)
end

do 
    local outpath = paths.concat(opt.output_dir, 'test')
    if not paths.dirp(outpath) then paths.mkdir(outpath) end
    test(test_files, outpath)
end
