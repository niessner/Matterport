
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