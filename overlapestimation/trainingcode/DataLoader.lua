require 'image'
local threads = require 'threads'
threads.Threads.serialization('threads.sharedserialize')

function file_exists(name)
   local f=io.open(name,"r")
   if f~=nil then io.close(f) return true else return false end
end

function lines_from(file)
  local lines = {}
  local numoflines = 0
  local maxfilename = 0
  if not file_exists(file) then 
    print('file not exists', file)
    return lines,numoflines,maxfilename
  end
  
  for line in io.lines(file) do 
    lines[#lines + 1] = line
    maxfilename = math.max(string.len(line),maxfilename)
    numoflines = numoflines + 1
  end
  maxfilename =maxfilename+1
  return lines,numoflines,maxfilename
end


local DataLoader = torch.class('DataLoader')

function  DataLoader:readOverLapFile(opt) 
    self.imagePairs = {}
    for index,sceneId in ipairs(self.sceneIds) do 
        -- local matchImg = image.load(self.dataRoot..'/'..sceneId..'/'..'/overlaps2/'..sceneId..'_iis_iou.png')
        local sceneId_c
        local sceneId_s = string.split(sceneId,'/')
        for i,s in ipairs(sceneId_s) do 
            sceneId_c = s
        end
         
        local matchImg = image.load(self.dataRoot..'/'..sceneId..'/'..'/overlaps/'..sceneId_c..'_iis_iou.png')
        local matchImg = matchImg:mul(65.536)
        self.imagePairs[index] = image.vflip(matchImg)
    end
end

function  DataLoader:readConfigFile(opt) 
    self.datasetSize = {}
    self.dataList = {{}}
    self.numofScene = 0;

    for index,sceneId in ipairs(self.sceneIds) do 
        local sceneId_c
        local sceneId_s = string.split(sceneId,'/')
        for i,s in ipairs(sceneId_s) do 
            sceneId_c = s
        end


        local lines,numoflines,maxPathLength = lines_from(self.dataRoot..'/'..sceneId..'/'..sceneId_c..'.conf')
        print(self.dataRoot..'/'..sceneId..'/'..sceneId_c..'.conf')
        if numoflines>0 then 
            local count = 0
            self.dataList[index] = {};
            for k=1,#lines do
                local linsplit = string.split(lines[k],' +')
                if linsplit[1] == 'scan' then 
                   count = count+1;
                   if opt.dataset =='sun3d'  then 
                      if count%10==1 then  -- overlap is compute every 10 frame
                         self.dataList[index][math.floor(count/10)+1] =  linsplit[3]
                      end
                   else
                      self.dataList[index][count] =  linsplit[3]
                   end

                end

                if linsplit[1] == 'image_directory' then 
                   self.image_directory = linsplit[2]
                end
                if linsplit[1] == 'color_directory' then 
                   self.image_directory = linsplit[2]
                end
                
            end
            if opt.dataset =='sun3d'  then 
                self.datasetSize[index] = math.floor(count/10)+1
            else
                self.datasetSize[index] = count
            end
            self.numofScene = self.numofScene+1
        end
    end    
    self.datasetTotalSize = torch.sum(torch.Tensor(self.datasetSize))
end

-- Initialize multi-thread data loader
function DataLoader:__init(opt)
    self.loadSize = opt.loadSize
    self.fineSize = opt.fineSize
    self.batchSize = opt.batchSize
    self.nc = opt.nc
    self.dataRoot = opt.dataRoot
    self.phase = opt.phase
    self.targetType  = opt.targetType
    self.sceneIds = {}
    if self.phase == 'train' then 
        -- List all paths to training images
        local lines,numoflines,maxPathLength = lines_from(opt.seqlist)
        for i =1,#lines do
            self.sceneIds[i] = lines[i]
        end
    else
        local lines,numoflines,maxPathLength = lines_from(opt.seqlist)
        for i =1,#lines do
            self.sceneIds[i] = lines[i]
        end
    end
    self.sceneCounter = 1
    self:readConfigFile(opt)
    print(self.datasetSize)
    self:readOverLapFile(opt)
end

function DataLoader:scanDir(directory,query)
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

function DataLoader:preprocessImg(img)
    img = image.scale(img,self.fineSize,self.fineSize)

    local mean = {0.485,0.456,0.406}
    local std = {0.229,0.224,0.225}
    for i=1,3 do
        img[i]:add(-mean[i])
        img[i]:div(std[i])
    end
    return img
end


function DataLoader:size()
    return self.datasetTotalSize
end
function DataLoader:SceneSize(insceneId)
    for i = 1,self.numofScene do
        if insceneId==self.sceneIds[i] then 
           return self.datasetSize[i]
        end
    end

end


-- Load training batch with multi-thread
function DataLoader:getBatch(inputtable)
    collectgarbage()
    local input_matchPatch, input_anchorPatch, input_nonMatchPatch, matchRatio


    if self.phase == 'train' then
        input_matchPatch = torch.zeros(self.batchSize,3,self.fineSize,self.fineSize)
        input_anchorPatch   = torch.zeros(self.batchSize,3,self.fineSize,self.fineSize)
        input_nonMatchPatch = torch.zeros(self.batchSize,3,self.fineSize,self.fineSize)
        matchRatio = torch.zeros(self.batchSize)
        
        for testIdx = 1,self.batchSize do
            -- find random pair 
            local rd_sIdx= math.floor(math.random()*self.numofScene)+1
            
            local pathtofolder = self.dataRoot..self.sceneIds[rd_sIdx]
            
            local rd_anchorIdx, rd_matchIdx,rd_nonmatchIdx,rd_idx,matchratio  = 0,0,0,0,0
            local matchcount = 0
            local matchArray


            while matchcount == 0 do 
                rd_anchorIdx =  math.floor(math.random()*self.imagePairs[rd_sIdx]:size(2))+1
                matchArray  = torch.zeros(self.imagePairs[rd_sIdx]:size(2))
                for i = 1,self.imagePairs[rd_sIdx]:size(2) do 
                    mean_matchratio = 0.5*(self.imagePairs[rd_sIdx][1][rd_anchorIdx][i]+self.imagePairs[rd_sIdx][1][i][rd_anchorIdx])
                    if mean_matchratio>0.1 and mean_matchratio < 0.85 then 
                       matchcount = matchcount+1
                       matchArray[matchcount] = i
                    end
                end
            end

            local pickid = math.floor(math.random()*matchcount)+1;
            local rd_matchIdx =  matchArray[pickid]
            
            matchratio = 0.5*(self.imagePairs[rd_sIdx][1][rd_anchorIdx][rd_matchIdx] + self.imagePairs[rd_sIdx][1][rd_matchIdx][rd_anchorIdx])
            matchRatio[testIdx] = matchratio

            while rd_nonmatchIdx == 0 do 
                rd_idx =  math.floor(math.random()*self.imagePairs[rd_sIdx]:size(2))+1
                matchratio = self.imagePairs[rd_sIdx][1][rd_anchorIdx][rd_idx]
                if matchratio<0.05 then
                   rd_nonmatchIdx = rd_idx
                   break
                end
            end

            -- print(pathtofolder)
            -- print(rd_sIdx..' '..rd_matchIdx ..' '.. rd_anchorIdx ..' '..rd_nonmatchIdx .. ' '.. matchRatio[testIdx] )
            -- print(self.dataList[rd_sIdx][rd_matchIdx])
            -- print(self.dataList[rd_sIdx][rd_anchorIdx])

            local matchImg = image.load(pathtofolder..'/'..self.image_directory..'/'.. self.dataList[rd_sIdx][rd_matchIdx],3)
            local anchorImg = image.load(pathtofolder..'/'..self.image_directory..'/'.. self.dataList[rd_sIdx][rd_anchorIdx],3)
            local nonMatchImg = image.load(pathtofolder..'/'..self.image_directory..'/'.. self.dataList[rd_sIdx][rd_nonmatchIdx],3)
            -- Preprocess image patches
            matchImg = self:preprocessImg(matchImg)
            anchorImg = self:preprocessImg(anchorImg)
            nonMatchImg = self:preprocessImg(nonMatchImg)

            
            -- Assign data to input
            input_matchPatch[testIdx] = matchImg
            input_anchorPatch[testIdx] = anchorImg
            input_nonMatchPatch[testIdx] = nonMatchImg
        end
    else
        local insceneId = inputtable[1]
        local rd_matchIdx = inputtable[2]
        local rd_anchorIdx = inputtable[3]
        local rd_sIdx = 0

        for i = 1,self.numofScene do
            if insceneId==self.sceneIds[i] then 
               rd_sIdx = i
               break
            end
        end

        local pathtofolder = self.dataRoot..self.sceneIds[rd_sIdx]
        local matchImg = image.load(pathtofolder..'/'..self.image_directory..'/'.. self.dataList[rd_sIdx][rd_matchIdx],3)
        local anchorImg = image.load(pathtofolder..'/'..self.image_directory..'/'.. self.dataList[rd_sIdx][rd_anchorIdx],3)    
        input_matchPatch = self:preprocessImg(matchImg)
        input_anchorPatch = self:preprocessImg(anchorImg)
        input_nonMatchPatch = nil
        matchRatio = nil 
    end
    collectgarbage()

    return input_matchPatch,input_anchorPatch,input_nonMatchPatch,matchRatio
end
