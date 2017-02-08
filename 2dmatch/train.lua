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

-- Load model and criterion
model,criterion = getModel()

-- Construct window for visualizing training examples
-- visWindow = qtwidget.newwindow(1792,672)

for trainIter = 1,1000 do

	-- Create a mini-batch of training examples
	for i = 1,8 do 

		-- Get training triplet of anchor patch, matching patch, and non-matching patch
		local matchPatch,anchorPatch,nonMatchPath = getTrainingExampleTriplet()

		-- Reshape data for input to network
		matchPatch = matchPatch:reshape(1,3,224,224);
		anchorPatch = anchorPatch:reshape(1,3,224,224);
		nonMatchPath = nonMatchPath:reshape(1,3,224,224);
		if i == 1 then
			input = {matchPatch,anchorPatch,nonMatchPath}
			-- mosaic = anchorPatch:reshape(3,224,224):cat(matchPatch:reshape(3,224,224),2):cat(nonMatchPath:reshape(3,224,224),2)
		else
			input[1] = input[1]:cat(matchPatch,1)
			input[2] = input[2]:cat(anchorPatch,1)
			input[3] = input[3]:cat(nonMatchPath,1)
			-- mosaic = mosaic:cat(anchorPatch:reshape(3,224,224):cat(matchPatch:reshape(3,224,224),2):cat(nonMatchPath:reshape(3,224,224),2),3)
		end
	end

	-- Convert input to GPU memory
	input[1] = input[1]:cuda()
	input[2] = input[2]:cuda()
	input[3] = input[3]:cuda()

	-- Visualize training examples (console command: qlua -lenv train.lua)
	-- image.display{image=mosaic,offscreen=false,win=visWindow}

	params,gradParams = model:getParameters()
	local feval = function(x)

	    -- Update model parameters
	    if x ~= params then
	        params:copy(x) 
	    end

	    -- Reset gradients
	    gradParams:zero() 

	    -- Forward and backward pass
	    local output = model:forward(input)
	    local loss = criterion:forward(output)
	    print('Training iteration '..trainIter..': '..loss)
	    local dLoss = criterion:backward(output)
	    model:backward(input,dLoss)
	    return loss,gradParams
	end

	-- Update model parameters (SGD)
	config = {learningRate = 0.01,momentum = 0.9}
	optim.sgd(feval,params,config)

	-- Save training snapshot of model
	if trainIter%100 == 0 then
		local filename = paths.concat("snapshots","snapshot_"..trainIter..".net")
		os.execute('mkdir -p '..sys.dirname(filename))
		torch.save(filename, model:clearState()) 
	end
end
