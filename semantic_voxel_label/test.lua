require 'nn'
require 'cunn'
require 'cudnn'
require 'xlua'
require 'torch'
require 'hdf5'
require 'paths'
require 'util'
dofile './provider.lua'

opt_string = [[
    --model         (default "logs/model.net")              torch model file path
    --h5_list_path  (default "")                            file of h5 data paths
    --gpu_index     (default 0)                             GPU index
    --output_dir    (default "output/")
    --num_classes   (default 42)
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

save_output_files = true
if save_output_files and not paths.dirp(opt.output_dir) then
    paths.mkdir(opt.output_dir)
end

-- set gpu
cutorch.setDevice(opt.gpu_index+1)

num_classes = opt.num_classes --42
local excluded_classes = {} --exclude empty/unannotated
excluded_classes[1] = true 
excluded_classes[num_classes] = true
print('using #classes = ' .. num_classes)
voxel_dimensions = torch.LongTensor{62, 31, 31}
surround_x = math.floor(voxel_dimensions[3]/2)
surround_y = math.floor(voxel_dimensions[2]/2)

print('Loading model...')
model = torch.load(opt.model):cuda()
model:evaluate()
print(model)


test_files = getDataFiles(opt.h5_list_path)
print('#test files = ' .. #test_files)

local confusion = torch.Tensor(num_classes, num_classes):fill(0)
count = 0
total = 0
countPerClass = torch.zeros(num_classes)
totalPerClass = torch.zeros(num_classes)
for i = 1,#test_files do
    --print('Loading scene')
    local scene_occ, scene_label = loadDataFile(test_files[i], num_classes, class_map)
    local output_scene = torch.zeros(voxel_dimensions[1], scene_label:size(2), scene_label:size(3))

    local outfile = paths.basename(test_files[i], '.h5') ..  '-dim_' .. tostring(output_scene:size(3)) .. '-' .. tostring(output_scene:size(2)) .. '-' .. tostring(output_scene:size(1)) .. '.bin'
    outfile = paths.concat(opt.output_dir, outfile)

    local scene_height = math.min(voxel_dimensions[1], scene_occ:size(1))
    local ystart = surround_y + 1
    local yend = scene_occ:size(2) - surround_y
    local xstart = surround_x + 1
    local xend = scene_occ:size(3) - surround_x

    local inputs = torch.CudaTensor(1, 1, voxel_dimensions[1], voxel_dimensions[2], voxel_dimensions[3])
    for y = ystart,yend do
        xlua.progress(y-ystart+1, yend-ystart+1)
        for x = xstart,xend do
            local target = torch.squeeze(scene_label[{{},{y},{x}}]:clone())
            local middle = scene_occ[{{},{y},{x}}]
            if torch.sum(middle) > 0 then --only predict for non-empty
                --collect sample
                local sample = torch.zeros(voxel_dimensions[1], voxel_dimensions[2], voxel_dimensions[3])
                sample[{{1,scene_height},{},{}}] = scene_occ[{{1,scene_height},{y-surround_y,y+surround_y},{x-surround_x,x+surround_x}}]
                inputs:copy(sample)
                local outputs = torch.squeeze(model:forward(inputs))
                local val, idx = torch.max(outputs:sub(1,num_classes-1):double(), 1)
                idx = torch.squeeze(idx)
                for k = 1,scene_height do
                    local tgt = target[k]
                    if excluded_classes[tgt] ~= true then
                        total = total + 1
                        totalPerClass[tgt] = totalPerClass[tgt] + 1
                        if idx[k] == tgt then 
                            count = count + 1
                            countPerClass[tgt] = countPerClass[tgt] + 1
                        end
                        confusion[tgt][idx[k]] = confusion[tgt][idx[k]] + 1
                        output_scene[k][y][x] = idx[k]-1
                    elseif tgt == num_classes then
                        output_scene[k][y][x] = idx[k]-1
                    end
                end
            end --non-empty
        end --x
    end --y
    if save_output_files then writeTensorToFile(output_scene, outfile) end
end --scenes
for r = 1,num_classes do
    for c = 1,num_classes do
        if totalPerClass[r] > 0 then
            confusion[r][c] = confusion[r][c] / totalPerClass[r]
        end
    end
end
saveConfusionToFile(confusion, paths.concat(opt.output_dir, 'confusion.txt'))
print('accuracy = ' .. count/total .. ' = ' .. count .. ' / ' .. total)
print('accuracy per class:')
for k = 1,num_classes do
    print('class ' .. k .. ':\t' .. countPerClass[k]/totalPerClass[k] .. '\t' .. countPerClass[k] .. '\t' .. totalPerClass[k])
end



