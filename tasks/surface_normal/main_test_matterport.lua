require 'nn'
require 'cutorch'
require 'cunn'
require 'cudnn'
require 'nngraph'
require 'optim'
require 'image'

require 'BatchIterator'
require 'utils'
require 'utils_matterport'
require 'BatchIterator_matterport'
-- require 'hdf5'

local config = dofile('config.lua')
config = config.parse(arg)
print(config)
cutorch.setDevice(config.gpuid)

local tmp1 = split(config.test_model, "/")
config.result_path = config.result_path .. "/" .. string.sub(tmp1[#tmp1],1,-4)


config.result_path = config.result_path .. "_matterport_test/"

os.execute("mkdir " .. config.result_path)

-- local model = dofile('model_multi_task.lua')(config.do_normal, config.do_semantic, config.do_boundary, config.do_room)
local model = dofile(config.model)(config) 

parameters, gradParameters = model:getParameters()
model:cuda()
parameters, gradParameters = model:getParameters()
parameters:copy(torch.load(config.test_model))

-- dataset
local train_data = {}
local test_data  = loadMatterport(config.test_file, config.root_path)
local batch_iterator = BatchIterator(config, train_data, test_data)
batch_iterator:setBatchSize(1)

local test_count = 0

while batch_iterator.epoch==0 and test_count<=config.max_count do
	local batch = batch_iterator:nextBatchMatterport('test', config)
	local currName = batch_iterator:currentName('test')
    print(currName)
	local k = split(currName, "/")
    saveName = k[#k-2] .. "_" .. k[#k]
	print(string.format("Testing %s", saveName))
	

	local inputs = batch.pr_color
    inputs = inputs:contiguous():cuda()
    local outputs = model:forward(inputs)

    local ch, h, w = 0, 0, 0
    local normal_est, normal_mask, normal_gnd, f_normal, df_do_normal, normal_outputs = nil,nil,nil,nil,nil,nil

    normal_est = outputs
    ch, h, w = normal_est:size(2), normal_est:size(3), normal_est:size(4)
    normal_est = normal_est:permute(1, 3, 4, 2):contiguous()
    normal_est = normal_est:view(-1, ch)
    local normalize_layer = nn.Normalize(2):cuda()
    normal_outputs = normalize_layer:forward(normal_est)
	normal_outputs = normal_outputs:view(1, h, w, ch)
	normal_outputs = normal_outputs:permute(1, 4, 2, 3):contiguous()
	normal_outputs = normal_outputs:view( ch, h, w)
	normal_outputs = normal_outputs:float()

	image.save(string.format("%s%s_normal_est.png", config.result_path, saveName), normal_outputs:add(1):mul(0.5))


    test_count = test_count + 1
end

print("Finish!")


