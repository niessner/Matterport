-- Some useful random util functions 

function loadMatchFile(file, skip)
    assert(paths.filep(file))

    local keypoints = {}

    local i = 0
    for line in io.lines(file) do    
        if i >= 3 then
            if (i-3)%skip == 0 or (i-3)%skip == 1 then 
                local parts = {}
                for p in line:gmatch('([^\t]+)') do table.insert(parts, p) end
                table.insert(keypoints, parts)
                --print(parts)
                --io.read()
            end
            
            --if i >= 10+2 then break end --todo remove this part here
            --if i % 100000 == 0 then
            --    print('i ' .. i)
            --    print(gcinfo()/1024 .. ' mb')
            --end
            --if #keypoints >= 100000 then break end
        end
        i = i + 1
    end

    return keypoints

end


function strToVec2(str) 
    local parts = {}
    for p in str:gmatch('%w+') do table.insert(parts, p) end
    assert(#parts == 2) 
    return torch.Tensor{tonumber(parts[1]), tonumber(parts[2])};
end

function inBounds(p, bounds, padding)
    for i = 1, p:size(1), 1 do
        if p[i] < padding or p[i] + padding > bounds[i] then
            return false
        end
    end
    return true
end


function loadMatchFiles(basePath, files, padding, skip, imWidth, imHeight, scaleX, scaleY)
    --load in the train data (positive and negative matches) 
    local poss = {}
    local negs = {}
    for fn = 1, #files do
        local sceneName = files[fn] 
        local file_pos = paths.concat(basePath, sceneName, 'matches.txt')
        local file_neg = paths.concat(basePath, sceneName, 'negatives.txt')

        local _pos = loadMatchFile(file_pos, skip)
        local _neg = loadMatchFile(file_neg, skip)

        assert(#_pos == #_neg)

        for i = 1, #_pos do 
            --[[local scale = 2.0    --because our images are only half the size
            _pos[i][4] = strToVec2(_pos[i][4]) / scale
            _neg[i][4] = strToVec2(_neg[i][4]) / scale--]]
            _pos[i][4] = strToVec2(_pos[i][4])
            _pos[i][4][1] = _pos[i][4][1] * scaleX; _pos[i][4][2] = _pos[i][4][2] * scaleY
            _neg[i][4] = strToVec2(_neg[i][4])
            _neg[i][4][1] = _neg[i][4][1] * scaleX; _neg[i][4][2] = _neg[i][4][2] * scaleY

            assert(_pos[i][1] == _neg[i][1])    --make sure the match index is the same
        end


        local bounds = torch.Tensor{imWidth, imHeight}
        for i = 1, #_pos, 2 do
            assert(_pos[i+0][1] == _pos[i+1][1])    --make sure the match index is the same
            assert(_neg[i+0][1] == _neg[i+1][1])    --make sure the match index is the same
            if 
                inBounds(_pos[i+0][4], bounds, padding) and
                inBounds(_pos[i+1][4], bounds, padding) and
                inBounds(_neg[i+0][4], bounds, padding) and
                inBounds(_neg[i+1][4], bounds, padding)
            then
                table.insert(poss, {sceneName, _pos[i], _pos[i+1]})
                table.insert(negs, {sceneName,_neg[i], _neg[i+1]})
                
                --[[
                print(_pos[i])
                print(_pos[i][4])
                print(_pos[i+1])
                print(_pos[i+1][4])
                print(_neg[i])
                print(_neg[i][4])
                print(_neg[i+1])
                print(_neg[i+1][4])
                io.read()--]]
            end
        end
    end
    return poss, negs
end

function markImage(img, locx, locy)
    local res = img:clone()
    res[1][locy][locx] = 1
    res[2][locy][locx] = 0
    res[3][locy][locx] = 0
    return res
end

-- Includes: matching tote image and product image, and non-matching tote image
function getTrainingExampleTriplet(path, kp_anc, kp_pos, kp_neg, patchSize)
    
    local str_anc = string.format("color-%02d-%06d.jpg", kp_anc[2], kp_anc[3])
    local str_pos = string.format("color-%02d-%06d.jpg", kp_pos[2], kp_pos[3])
    local str_neg = string.format("color-%02d-%06d.jpg", kp_neg[2], kp_neg[3])

    --local t =  torch.tic()
    local anchorImg = image.load(paths.concat(path,str_anc),3,'float')
    local matchImg = image.load(paths.concat(path,str_pos),3,'float')
    local nonMatchImg = image.load(paths.concat(path,str_neg),3,'float')
    --print('image load time:' , torch.toc(t) * 1000.0 .. ' ms')

    -- Pixel locations of patch centers (x,y)
    local anchorPixelLoc = torch.floor(kp_anc[4])
    local matchPixelLoc = torch.floor(kp_pos[4])
    local nonMatchPixelLoc = torch.floor(kp_neg[4])

    --[[--debugging
    image.save('orig_anchor.png', markImage(anchorImg, anchorPixelLoc[1], anchorPixelLoc[2]))
    image.save('orig_match.png', markImage(matchImg, matchPixelLoc[1], matchPixelLoc[2]))
    image.save('orig_nonmatch.png', markImage(nonMatchImg, nonMatchPixelLoc[1], nonMatchPixelLoc[2]))
    --debugging--]]
--[[
    print(str_anc .. '\t' .. anchorPixelLoc[1] .. ',' .. anchorPixelLoc[2])
    print(str_pos .. '\t' .. matchPixelLoc[1] .. ',' .. matchPixelLoc[2])
    print(str_neg .. '\t' .. nonMatchPixelLoc[1] .. ',' .. nonMatchPixelLoc[2])
--]]
    -- Extract 64x64 patches
    local matchPatch = image.crop(matchImg,matchPixelLoc[1]-patchSize/2,matchPixelLoc[2]-patchSize/2,matchPixelLoc[1]+patchSize/2,matchPixelLoc[2]+patchSize/2)
    local anchorPatch = image.crop(anchorImg,anchorPixelLoc[1]-patchSize/2,anchorPixelLoc[2]-patchSize/2,anchorPixelLoc[1]+patchSize/2,anchorPixelLoc[2]+patchSize/2)
    local nonMatchPatch = image.crop(nonMatchImg,nonMatchPixelLoc[1]-patchSize/2,nonMatchPixelLoc[2]-patchSize/2,nonMatchPixelLoc[1]+patchSize/2,nonMatchPixelLoc[2]+patchSize/2)

    --[[--debugging
    image.save('patch_anchor.png', markImage(anchorPatch, patchSize/2+1, patchSize/2+1))
    image.save('patch_match.png', markImage(matchPatch, patchSize/2+1, patchSize/2+1))
    image.save('patch_nonmatch.png', markImage(nonMatchPatch, patchSize/2+1, patchSize/2+1))
    --debugging --]]
    
    -- Preprocess image patches
    anchorPatch = preprocessImg(anchorPatch)
    matchPatch = preprocessImg(matchPatch)
    nonMatchPatch = preprocessImg(nonMatchPatch)

    return anchorPatch,matchPatch,nonMatchPatch
end

--remove trailing/leading whitespace from string
function trim(s)
    return (s:gsub("^%s*(.-)%s*$", "%1"))
end

-- read h5 filename list
-- specify base path to filter out scenes that don't exist
function getDataFiles(input_file, base_path)
    assert(paths.filep(input_file))
    local data_files = {}
    for line in io.lines(input_file) do
        local scene = trim(line)
        if base_path then
            if paths.filep(paths.concat(base_path, scene, 'matches.txt')) and paths.filep(paths.concat(base_path, scene, 'negatives.txt')) and paths.dirp(paths.concat(base_path, scene, 'images')) then
                data_files[#data_files+1] = scene
            else
                print('warning: skipping non-existent scene ' .. scene)
            end
        else 
            data_files[#data_files+1] = scene
        end
    end
    return data_files
end

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
    if img:size(1) ~= 224 then
        img = image.scale(img,224,224)
    end
    --local mean = {0.485,0.456,0.406}
    --local std = {0.229,0.224,0.225}
    local mean = {0.524978, 0.510486, 0.468835}
    local std = {0.2088, 0.2155, 0.2415}
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

function recursiveModelReset(model)
    for i = 1,model:size() do
        local tmpLayer = model:get(i)
        --print(i)
        --print(tmpLayer)
        if torch.type(tmpLayer):find('Convolution') or torch.type(tmpLayer):find('Linear') then
            --print('reset')
            tmpLayer:reset()
        end
        if torch.type(tmpLayer):find('Sequential') or torch.type(tmpLayer):find('ConcatTable') or torch.type(tmpLayer):find('ParallelTable') then
            recursiveModelReset(tmpLayer)
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



function serialize (o)
    if type(o) == "number" then
        io.write(o)
    elseif type(o) == "string" then
        io.write(string.format("%q", o))
    elseif type(o) == "boolean" then
        io.write(tostring(o))
    elseif type(o) == "table" then
        io.write("{\n")
        for k,v in pairs(o) do
            io.write("  ", k, " = ")
            serialize(v)
            io.write(",\n")
        end
        io.write("}\n")
    else
        error("cannot serialize a " .. type(o))
    end
end
