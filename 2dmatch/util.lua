-- Some useful random util functions 

-- Lookup filenames in directory (with search query string)
function scanDir(directory,query)
    local i, t, popen = 0, {}, io.popen
    local pfile = popen('ls -a "'..directory..'"')
    for filename in pfile:lines() do
        if string.find(filename,query) then
            i = i+1
            t[i] = filename
        end
    end
    pfile:close()
    return t
end

-- Pre-process images for ResNet-101 pre-trained on ImageNet (224x224 RGB mean-subtracted std-divided)
function preprocessImg(img)
    img = image.scale(img,224,224)
    local mean = {0.485,0.456,0.406}
    local std = {0.229,0.224,0.225}
    for i=1,3 do
        img[i]:add(-mean[i])
        img[i]:div(std[i])
    end
    return img

    -- Old code (please double check)
    -- img = image.scale(img,224,224)
    -- local mean_pixel = torch.DoubleTensor({123.68, 116.779, 103.939})
    -- img = img:mul(255.0)+1
    -- mean_pixel = mean_pixel:view(3,1,1):expandAs(img)
    -- img:add(-1,mean_pixel)
    -- return img
end

-- Includes: matching tote image and product image, and non-matching tote image
function getTrainingExampleTriplet()
    local matchImg = image.load(paths.concat('data','sample','frame-000001.color.jpg'),3,'float')
    local anchorImg = image.load(paths.concat('data','sample','frame-000000.color.jpg'),3,'float')
    local nonMatchImg = image.load(paths.concat('data','sample','frame-000000.color.jpg'),3,'float')

    -- Pixel locations of patch centers (x,y)
    local matchPixelLoc = {821,109}
    local anchorPixelLoc = {808,897}
    local nonMatchPixelLoc = {661,794}

    -- Extract 64x64 patches
    local patchSize = 64
    local matchPatch = image.crop(matchImg,matchPixelLoc[1]-patchSize/2,matchPixelLoc[2]-patchSize/2,matchPixelLoc[1]+patchSize/2,matchPixelLoc[2]+patchSize/2)
    local anchorPatch = image.crop(anchorImg,anchorPixelLoc[1]-patchSize/2,anchorPixelLoc[2]-patchSize/2,anchorPixelLoc[1]+patchSize/2,anchorPixelLoc[2]+patchSize/2)
    local nonMatchPatch = image.crop(nonMatchImg,nonMatchPixelLoc[1]-patchSize/2,nonMatchPixelLoc[2]-patchSize/2,nonMatchPixelLoc[1]+patchSize/2,nonMatchPixelLoc[2]+patchSize/2)

    -- Preprocess image patches
    matchPatch = preprocessImg(matchPatch)
    anchorPatch = preprocessImg(anchorPatch)
    nonMatchPatch = preprocessImg(nonMatchPatch)

    return matchPatch,anchorPatch,nonMatchPatch
end

-- Get subset of a 1D table
function subrange(t,first,last)
  local sub = {}
  for i = first,last do
    sub[#sub+1] = t[i]
  end
  return sub
end

-- Recursively freeze layers of the model
function recursiveModelFreeze(model)
    for i = 1,model:size() do
        local tmpLayer = model:get(i)
        if torch.type(tmpLayer):find('Convolution') or torch.type(tmpLayer):find('Linear') then

            -- Set parameter update functions to empty functions
            tmpLayer.accGradParameters = function() end
            tmpLayer.updateParameters = function() end
        end
        if torch.type(tmpLayer):find('Sequential') or torch.type(tmpLayer):find('ConcatTable') then
            recursiveModelFreeze(tmpLayer)
        end
    end
end

-- Load depth file (saved as 16-bit PNG in centimeters)
function loadDepth(filename)
    depth = image.load(filename)*65536/10000
    depth = depth:clamp(0.2,1.2) -- Depth range of Intel RealSense F200
    depth = depth:csub(0.440931) -- Subtract average mean depth value from training data
    return depth
end