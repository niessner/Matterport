require 'cutorch'
require 'cunn'
require 'cudnn'

function getModel(numClasses)

    -- Load ResNet-50 pre-trained on ImageNet
    local model = torch.load('resnet-50.t7')
    model:remove(11)
    model:add(nn.Linear(2048,numClasses))
 
    -- Add classification loss layer
    local criterion = nn.CrossEntropyCriterion()

    return model,criterion
end