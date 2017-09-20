require 'image'
require 'cutorch'
require 'cunn'
require 'cudnn'
require 'optim'

-- Custom files
require 'model'
require 'util'
-- require 'qtwidget' -- for visualizing images

-- Set RNG seed
math.randomseed(os.time())

local basePath = '/mnt/raid/datasets/Matterport/Matching/17DRP5sb8fy'
local file_pos = paths.concat(basePath, 'matches.txt')
local file_neg = paths.concat(basePath, 'negatives.txt')
local imgPath = paths.concat(basePath, 'images/') 

local poss = loadMatchFile(file_pos)
local negs = loadMatchFile(file_neg)
--print(matches)


for i = 1, #poss, 2 do 
	print(i)
	print(poss[i])
	print(negs[i])
	--print(matches[i+1])

	local patch_anc, patch_pos, patch_neg = 
		getTrainingExampleTriplet(imgPath, poss[i], poss[i+1], negs[i+1])  
end

