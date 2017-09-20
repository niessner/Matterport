require 'cutorch'
require 'cunn'
require 'cudnn'

-- Build model (triplet ResNet-101)
function getModel(useBottleneck)--, useMetricNetwork)
	
	local featureNet = nn.Sequential()
	featureNet:add(nn.SpatialConvolution(1, 24, 7, 7, 1, 1, 3, 3))		--conv0
	featureNet:add(nn.ReLU())
	featureNet:add(nn.SpatialMaxPooling(3, 3, 2, 2, 1, 1))				--pool0
	featureNet:add(nn.SpatialConvolution(24, 64, 5, 5, 1, 1, 2, 2))		--conv1
	featureNet:add(nn.ReLU())
	featureNet:add(nn.SpatialMaxPooling(3, 3, 2, 2, 1, 1))				--pool1
	featureNet:add(nn.SpatialConvolution(64, 96, 3, 3, 1, 1, 1, 1))		--conv2
	featureNet:add(nn.ReLU())
	featureNet:add(nn.SpatialConvolution(96, 96, 3, 3, 1, 1, 1, 1))		--conv3
	featureNet:add(nn.ReLU())
	featureNet:add(nn.SpatialConvolution(96, 64, 3, 3, 1, 1, 1, 1))		--conv4
	featureNet:add(nn.ReLU())
	featureNet:add(nn.SpatialMaxPooling(3, 3, 2, 2, 1, 1))				--pool4
	local fsize = 4096
	featureNet:add(nn.View(fsize))
	if useBottleneck then --bottleneck
		featureNet:add(nn.Linear(fsize, 512))
		featureNet:add(nn.ReLU())
		fsize = 512
	end
	featureNet:add(nn.Normalize(2))--normalize descriptor
	
	local model = nn.Sequential()
	local featureEncoder = nn.ParallelTable()
	featureEncoder:add(featureNet)
	featureEncoder:add(featureNet:clone('weight','bias', 'gradWeight','gradBias')) 
	model:add(featureEncoder)
	
	local criterion
	--[[if useMetricNetwork then--metric network
		model:add(nn.JoinTable(2))
		model:add(nn.Linear(fsize*2, 512))								--fc1
		model:add(nn.ReLU())
		model:add(nn.Linear(512, 512))									--fc2
		model:add(nn.ReLU())
		model:add(nn.Linear(512, 2))									--fc3
		
		criterion = nn.CrossEntropyCriterion()
	else--]] --L2 pariwise distance
		model:add(nn.PairwiseDistance(2))
		criterion = nn.HingeEmbeddingCriterion(1)
		--criterion = nn.CosineEmbeddingCriterion(0.8)
	--[[end--]]
	
	return model,criterion
end








