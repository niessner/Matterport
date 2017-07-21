require 'nn'
require 'cunn'
require 'cudnn'
require 'nngraph'

-- cnn: VGG-16 with batch norm
local function create_conv_unit(c1, c2, p, k, d)
    local conv = nn.Sequential()
    if d == 1 then
        conv:add(cudnn.SpatialConvolution(c1, c2, k, k, 1, 1, p, p, 1))
    else
        conv:add(nn.SpatialDilatedConvolution(c1, c2, k, k, p, p, d, d))
    end
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    return conv
end


local function create_conv_2(c1, c2)
    dilate = d or 1
    local conv = nn.Sequential()
    conv:add(cudnn.SpatialConvolution(c1, c2, 3, 3, 1, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    conv:add(cudnn.SpatialConvolution(c2, c2, 3, 3, 1, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    return conv
end

local function create_conv_3(c1, c2)
    dilate = d or 1
    local conv = nn.Sequential()
    conv:add(cudnn.SpatialConvolution(c1, c2, 3, 3, 1, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    conv:add(cudnn.SpatialConvolution(c2, c2, 3, 3, 1, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    conv:add(cudnn.SpatialConvolution(c2, c2, 3, 3, 1, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    return conv
end

local function create_deconv_3(c1, c2)
    local conv = nn.Sequential()
    conv:add(nn.SpatialFullConvolution(c1, c2, 3, 3, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    conv:add(nn.SpatialFullConvolution(c2, c2, 3, 3, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    conv:add(nn.SpatialFullConvolution(c2, c2, 3, 3, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    return conv
end

local function create_deconv_2(c1, c2)
    local conv = nn.Sequential()
    conv:add(nn.SpatialFullConvolution(c1, c2, 3, 3, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    conv:add(nn.SpatialFullConvolution(c2, c2, 3, 3, 1, 1, 1, 1))
    conv:add(cudnn.SpatialBatchNormalization(c2, nil, nil, nil))
    conv:add(cudnn.ReLU(true))
    return conv
end

-- local function create_simple_conv(c1,c2,)

local function create_model(config)
    -- local conv1_1 = create_conv_unit(3,64,3,1,1)
    -- local conv1_2 = create_conv_unit(64,64,3,1,1)

    -- local conv2_1 = create_conv_unit(64,128,3,1,1)
    -- local conv2_2 = create_conv_unit(128,128,3,1,1)

    -- local conv3_1 = create_conv_unit(128,256,3,1,1)
    -- local conv3_2 = create_conv_unit(256,256,3,1,1)
    -- local conv3_3 = create_conv_unit(256,256,3,1,1)

    -- local conv4_1 = create_conv_unit(256,512,3,1,1)
    -- local conv4_2 = create_conv_unit(512,512,3,1,1)
    -- local conv4_3 = create_conv_unit(512,512,3,1,1)

    -- local conv5_1 = create_conv_unit(512,512,3,1,1)
    -- local conv5_2 = create_conv_unit(512,512,3,1,1)
    -- local conv5_3 = create_conv_unit(512,512,3,1,1)

    local input_channel = config.input_channel    
    local output_channel = config.output_channel


    -- encoder
    local conv1 = create_conv_2(input_channel, 64)
    local conv2 = create_conv_2(64, 128)
    local conv3 = create_conv_3(128, 256)
    local conv4 = create_conv_3(256, 512)
    local conv5 = create_conv_3(512, 512)

    local pool1 = nn.SpatialMaxPooling(2, 2, 2, 2, 0, 0)
    local pool2 = nn.SpatialMaxPooling(2, 2, 2, 2, 0, 0)
    local pool3 = nn.SpatialMaxPooling(2, 2, 2, 2, 0, 0)
    local pool4 = nn.SpatialMaxPooling(2, 2, 2, 2, 0, 0)

    local input = nn.Identity()()
    local features1 = conv1(input) -- 3 -> 64
    local features2 = conv2(pool1(features1)) -- 64 -> 128
    local features3 = conv3(pool2(features2)) -- 128 ->256
    local features4 = conv4(pool3(features3)) -- 256 -> 512
    local features5 = conv5(pool4(features4)) -- 512 -> 512

    local deconv5 = create_deconv_3(512, 512)
    local deconv4 = create_deconv_3(512+512, 256)
    local deconv3 = create_deconv_3(256+256, 128)
    local deconv2 = create_deconv_2(128+128, 64)
    local deconv1 = create_deconv_2(64+64, output_channel)
    deconv1:remove(6)

    local defeature5 = deconv5(features5)
    local defeature4 = nn.JoinTable(1,3)({nn.SpatialMaxUnpooling(pool4)(defeature5),features4})
    local defeature3t = deconv4(defeature4)
    local defeature3 = nn.JoinTable(1,3)({nn.SpatialMaxUnpooling(pool3)(defeature3t),features3})
    local defeature2t = deconv3(defeature3)
    local defeature2 = nn.JoinTable(1,3)({nn.SpatialMaxUnpooling(pool2)(defeature2t),features2})
    local defeature1t = deconv2(defeature2)
    local defeature1 = nn.JoinTable(1,3)({nn.SpatialMaxUnpooling(pool1)(defeature1t),features1})
    local output = deconv1(defeature1)

    local model = nn.gModule({input}, {output})
    return model
end

return create_model
