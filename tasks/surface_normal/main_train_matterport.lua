require 'nn'
require 'cutorch'
require 'cunn'
require 'cudnn'
require 'nngraph'
require 'optim'

require 'BatchIterator'
require 'utils'
require 'utils_matterport'
require 'BatchIterator_matterport'


-- config
local config = dofile('config.lua')
-- print(arg)
config = config.parse(arg)
-- print(config)
cutorch.setDevice(config.gpuid)
print("Start: " .. config.ps)

-- model
local model = dofile(config.model)(config) 
parameters, gradParameters = model:getParameters()
model:cuda()
parameters, gradParameters = model:getParameters()

-- resume training
if config.resume_training then
    print('loading saved model weight...')
    parameters:copy(torch.load(config.saved_model_weights))
    config.optim_state = torch.load(config.saved_optim_state)
end

if config.finetune then
    print('finetune from saved model weight...')
    parameters:copy(torch.load(config.finetune_model))
    print('set up learning rate...')
    config.optim_state.learningRate = config.finetune_init_lr
end

-- criterion
local criterion_n = nn.CosineEmbeddingCriterion():cuda()

-- dataset
local train_data = loadMatterport(config.train_file, config.root_path)
local test_data  = loadMatterport(config.test_file, config.root_path)
local batch_iterator = BatchIterator(config, train_data, test_data)

-- logger
local logger = optim.Logger(config.log_path .. 'log', true)

-- main training
for it_batch = 1, math.floor(config.nb_epoch * #batch_iterator.train.data / config.batch_size) do

    -- print(">>")
    local batch = batch_iterator:nextBatchMatterport('train', config)
    -- print(batch.norm_valid:size())
    -- print("batch done")
    -- print('4')
    -- print(batch_iterator:currentName('train'))
    -- image.save('color.png', batch.pr_color[1]:add(0.5))
    -- image.save('norml.png', batch.cam_normal[1]:add(1):mul(0.5))
    -- image.save('valid.png', batch.norm_valid[1])
    -- image.load()
    -- print(batch.input:size())
    -- local temp = batch.input[{1,{},{},{}}]:view(3,240,320)
    -- image.save('color1.png', temp:add(0.5))
    -- temp = batch.input[{3,{},{},{}}]:view(3,240,320)
    -- image.save('color3.png', temp:add(0.5))
    -- temp = batch.input[{5,{},{},{}}]:view(3,240,320)
    -- image.save('color5.png', temp:add(0.5))

    -- inputs and targets
    local inputs = batch.pr_color
    inputs = inputs:contiguous():cuda()
    -- print(inputs:size())
    
    local feval = function(x)
        -- prepare
        collectgarbage()
        if x ~= parameters then
            parameters:copy(x)
        end
        
        -- forward propagation
        -- print('111')
        local est = model:forward(inputs)
        -- print(est:size())
        local valid = batch.norm_valid
        valid = valid:cuda()
        local gnd = batch.cam_normal
        gnd = gnd:cuda()

        bz, ch, h, w = est:size(1), est:size(2), est:size(3), est:size(4)
        est = est:permute(1,3,4,2):contiguous():view(-1,ch)
        local normalize_layer = nn.Normalize(2):cuda()
        est_n = normalize_layer:forward(est)
        gnd = gnd:permute(1,3,4,2):contiguous():view(-1,ch)

        -- print(est_n:size())
        -- print(gnd:size())

        f = criterion_n:forward({est_n, gnd}, torch.Tensor(est_n:size(1)):cuda():fill(1))
        df = criterion_n:backward({est_n, gnd}, torch.Tensor(est_n:size(1)):cuda():fill(1))
        -- print(df)
        df = df[1]
        -- print(df:size())
        df = normalize_layer:backward(est, df)


        -- print(valid:size())
        valid = valid:view(-1,1):expandAs(df)

        -- print(valid:size())
        df[torch.eq(valid,0)] = 0

        df = df:view(-1, h, w, ch)
        df = df:permute(1, 4, 2, 3):contiguous()

        gradParameters:zero()
        model:backward(inputs, df)

        -- print
        if it_batch % config.print_iters == 0 then
            print( it_batch, f)
        end

        -- log
        if it_batch % config.log_iters == 0 then
            -- logger:add{['normal_loss'] = f}
            -- logger:add{['segmentation_loss'] = fs}
            -- logger:add{ f_normal, f_semantic, f_boundary, f_room}
            logger:add{ f }
        end

        -- return
        -- return f_normal + f_semantic + f_boundary + f_room, gradParameters
        return f, gradParameters

    end

    -- optimizer
    optim.rmsprop(feval, parameters, config.optim_state)

    -- save
    if it_batch % config.snapshot_iters == 0 then
        print('saving model weight...')
        local filename
        filename = config.ps .. 'iter_' .. it_batch .. '.t7'
        torch.save(filename, parameters)
        filename = config.ps .. 'iter_' .. it_batch .. '_state.t7'
        torch.save(filename, config.optim_state)
    end

    -- lr
    if it_batch % config.lr_decay == 0 then
        config.optim_state.learningRate = config.optim_state.learningRate / config.lr_decay_t
        config.optim_state.learningRate = math.max(config.optim_state.learningRate, config.optim_state.learningRateMin)
        print('decresing lr... new lr:', config.optim_state.learningRate)
    end
end

print('saving model weight...')
local filename
filename = config.ps .. 'final' .. '.t7'
torch.save(filename, parameters)
filename = config.ps .. 'final' .. '_state.t7'
torch.save(filename, config.optim_state)
