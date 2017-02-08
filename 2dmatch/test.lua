require 'cutorch';
require 'cunn';
require 'cudnn'
require 'image'
require 'util'
-- require 'qtwidget'

-- Load training snapshot model
model = torch.load('snapshots/snapshot_100.net')

-- Extract one tower from the training model (3 cloned towers)
patchEncoder = model:get(1):get(1):clone()

-- Set to testing mode (disable dropout and batch normalization)
patchEncoder:evaluate()

-- Get example triplet patches
matchPatch,anchorPatch,nonMatchPatch = getTrainingExampleTriplet()

-- Compute features per patch
matchFeat = patchEncoder:forward(matchPatch:cuda()):clone()
anchorFeat = patchEncoder:forward(anchorPatch:cuda()):clone()
nonMatchFeat = patchEncoder:forward(nonMatchPatch:cuda()):clone()

-- Print descriptor distance between matching patches
print('Distance between matching patches: ' .. anchorFeat:dist(matchFeat))

-- Print descriptor distances between non-matching patches
print('Distance between non-matching patches 1: ' .. matchFeat:dist(nonMatchFeat))
print('Distance between non-matching patches 2: ' .. anchorFeat:dist(nonMatchFeat))
