require 'cutorch'
require 'cunn'
require 'cudnn'

-- Build model (triplet ResNet-101)
function getModel()

    -- Load ResNet-101 pre-trained on ImageNet
    local resnetModel = torch.load('resnet-50.t7')
    resnetModel:remove(11) -- Remove last FC layer
    resnetModel:insert(nn.Linear(2048,512)) -- Reduce descriptor vector size
    resnetModel:insert(nn.Normalize(2)) -- Normalize descriptor vector
    resnetModel:training()

    -- Load 3 copies of ResNet-101 in a triple tower model:
    -- First tower: matching patch
    -- Second tower: anchor patch
    -- Third tower: non-matching patch
    local tripletModel = nn.ParallelTable()
    tripletModel:add(resnetModel)
    tripletModel:add(resnetModel:clone('weight','bias', 'gradWeight','gradBias')) 
    tripletModel:add(resnetModel:clone('weight','bias', 'gradWeight','gradBias'))

    -- Build pairwise sample distance model (L2)
    local posDistModel = nn.Sequential() -- Similar sample distance w.r.t anchor sample
    posDistModel:add(nn.NarrowTable(1,2)):add(nn.PairwiseDistance(2))
    local negDistModel = nn.Sequential() -- Different sample distance w.r.t anchor sample
    negDistModel:add(nn.NarrowTable(2,2)):add(nn.PairwiseDistance(2))
    local distanceModel = nn.ConcatTable():add(posDistModel):add(negDistModel)

    -- Build complete model
    local model = nn.Sequential():add(tripletModel):add(distanceModel)
    return model
end
function getModelRandom()
    
end
function getModelRatio()

    -- Load ResNet-101 pre-trained on ImageNet
    local resnetModel = torch.load('resnet-50.t7')
    resnetModel:remove(11) -- Remove last FC layer
    resnetModel:insert(nn.Linear(2048,512)) -- Reduce descriptor vector size
    resnetModel:insert(nn.Normalize(2)) -- Normalize descriptor vector
    resnetModel:training()

    -- Load 3 copies of ResNet-101 in a triple tower model:
    -- First tower: matching patch
    -- Second tower: anchor patch
    -- Third tower: non-matching patch
    local tripletModel = nn.ParallelTable()
    tripletModel:add(resnetModel)
    tripletModel:add(resnetModel:clone('weight','bias', 'gradWeight','gradBias')) 
    tripletModel:add(resnetModel:clone('weight','bias', 'gradWeight','gradBias'))

    -- Build pairwise sample distance model (L2)
    local posDistModel = nn.Sequential() -- Similar sample distance w.r.t anchor sample
    posDistModel:add(nn.NarrowTable(1,2)):add(nn.PairwiseDistance(2))
    local negDistModel = nn.Sequential() -- Different sample distance w.r.t anchor sample
    negDistModel:add(nn.NarrowTable(2,2)):add(nn.PairwiseDistance(2))

    local posRatioModel = nn.Sequential() -- Different sample distance w.r.t anchor sample
    posRatioModel:add(nn.NarrowTable(1,2)):add(nn.JoinTable(2))
    posRatioModel:add(nn.Linear(1024,256)):add(nn.LeakyReLU())
    posRatioModel:add(nn.Linear(256,1)):add(nn.Clamp(0,1))   
    -- local negRatioModel = nn.Sequential() -- Different sample distance w.r.t anchor sample
    -- negRatioModel:add(nn.NarrowTable(2,2)):add(nn.JoinTable(2)):add(nn.Linear(256,1))


    local distanceModel = nn.ConcatTable():add(posDistModel):add(negDistModel):add(posRatioModel)

    -- Build complete model
    local model = nn.Sequential():add(tripletModel):add(distanceModel)

    

    return model
end






