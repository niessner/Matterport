require 'cutorch';
require 'cunn';
require 'cudnn'
require 'image'
require 'model'
dofile('util.lua')

-- example usage:  th test.lua --gpu_index 0 --model model.net --test_data scenes_test.txt

local basePath = 'data/keypointmatch'

opt_string = [[
    -h,--help                                                print help
    --gpu_index             (default 0)                      GPU index
    --model                 (default "model.net")            torch model file path
    --test_data             (default "scenes_test.txt")      txt file containing test
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
local patchSize = opt.patchSize
print('loading model: ' .. opt.model .. '...')
local model = torch.load(opt.model)
local patchEncoder = model:get(1):get(1):clone()
patchEncoder:evaluate()

-- load test files
local test_files = getDataFiles(paths.concat(basePath,opt.test_data), basePath)
print(test_files)

function test(data_files, outpath)
    assert(paths.dirp(outpath))
    
    --load in the data (positive and negative matches) 
    local poss, negs = loadMatchFiles(basePath, data_files, patchSize/2, opt.matchFileSkip)
    
    -- output files
    splitter = ',' --output as csv
    local outfile_anc = assert(io.open(paths.concat(outpath, 'anc_feat.txt'), "w"))
    local outfile_pos = assert(io.open(paths.concat(outpath, 'pos_feat.txt'), "w"))
    local outfile_neg = assert(io.open(paths.concat(outpath, 'neg_feat.txt'), "w"))
    
    for k = 1,#poss do
        -- print progress bar :D        
        xlua.progress(k, #poss)
        
        local sceneName = poss[k][1]
        local imgPath = paths.concat(basePath,sceneName,'images')
        
        -- Get example triplet patches
        local anc,pos,neg = getTrainingExampleTriplet(imgPath, poss[k][2], poss[k][3], negs[k][3], patchSize)
        
        -- Compute features per patch
        local ancFeat = torch.squeeze(patchEncoder:forward(anc:resize(1,3,224,224):cuda()):float())
        local posFeat = torch.squeeze(patchEncoder:forward(pos:resize(1,3,224,224):cuda()):float())
        local negFeat = torch.squeeze(patchEncoder:forward(neg:resize(1,3,224,224):cuda()):float())
        --print(ancFeat:size())
        assert(ancFeat:size(1) == posFeat:size(1) and ancFeat:size(1) == negFeat:size(1))

        -- -- descriptor distance between patches
        --print('Distance between matching patches: ' ..       ancFeat:dist(posFeat))
        --print('Distance between non-matching patches 1: ' .. posFeat:dist(negFeat))
        --print('Distance between non-matching patches 2: ' .. ancFeat:dist(negFeat))
        --io.read()
        local distAncPos = ancFeat:dist(posFeat)
        local distAncNeg = ancFeat:dist(negFeat)        
        -- log to file
        for i=1,ancFeat:size(1) do
            outfile_anc:write(string.format("%.6f", ancFeat[i]))
            outfile_pos:write(string.format("%.6f", posFeat[i]))
            outfile_neg:write(string.format("%.6f", negFeat[i]))
            if i < ancFeat:size(1) then 
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

-- Test!
if not paths.dirp(opt.output_dir) then paths.mkdir(opt.output_dir) end
do 
    local outpath = paths.concat(opt.output_dir, 'test')
    if not paths.dirp(outpath) then paths.mkdir(outpath) end
    test(test_files, outpath)
end