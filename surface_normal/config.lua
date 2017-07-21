require "utils"
--- All parameters goes here
local config = config or {}

function config.parse(arg)
	local cmd = torch.CmdLine()
	cmd:text()
	cmd:text('Multi-Task Classification FCN')
	cmd:text()
	-- Parameters

    -- model configuration
    cmd:option('-model', 'model_deep.lua', 'model file')
    cmd:option('-input_channel', 3, '# of input channels')
    cmd:option('-output_channel', 3, '# of output channels')

    -- testing
    cmd:option('-test_model', '', 'model used for testing')
    cmd:option('-result_path', './result/', 'path to save result')
    cmd:option('-max_count', 1000000, 'max number of data to test')

    -- data loader
    cmd:option('-root_path', '~/matterport/v1/', 'path to the root of the dataset')
    cmd:option('-train_file', './data_list/train_list_noup.txt', 'train file, compulsory');
    cmd:option('-test_file',  './data_list/test_list_horizontal_small.txt', 'test file, compulsory');
	
    -- training and testing
    cmd:option('-gpuid', 1, 'gpu id')
    cmd:option('-optim_state', {rho=0.95, eps=1e-6, learningRate=1e-3, learningRateMin=1e-7, momentum=0.9}, 'optim state')
    cmd:option('-lr_decay', 150000, 'iterations between lr decreses')
    cmd:option('-lr_decay_t', 5, 'lr decay times')
    cmd:option('-nb_epoch', 20, 'number of epoches')
    cmd:option('-batch_size', 1, 'batch size')
    cmd:option('-pixel_means', {128, 128, 128}, 'Pixel mean values (RGB order)')

    -- resume
    cmd:option('-resume_training', false, 'whether resume training')
    cmd:option('-saved_model_weights', '', 'path to saved model weights')
    cmd:option('-saved_optim_state', '', 'path to saved model weights')

    -- finetune
    cmd:option('-finetune', false, '')
    cmd:option('-finetune_model', '', '')
    cmd:option('-finetune_init_lr', 1e-4, '')

    -- save/print/log
	cmd:option('-snapshot_iters', 10000, 'Iterations between snapshots (used for saving the network)')
	cmd:option('-print_iters', 20, 'Iterations between print')
	cmd:option('-log_iters', 20, 'Iterations between log')
	cmd:option('-log_path','./logs/','Path to be used for logging')
    cmd:option('-ps', '', 'prefix: path&name to model and snapshot')
    cmd:option('-verbose', false, 'show more message')

	-- Parsing the command line 
	config = cmd:parse(arg or {})
    config.colors = {{0, 0, 0}, -- black 
                     {1, 0, 0}, -- red
                     {0, 1, 0}, -- green
                     {0, 0, 1}, -- blue
                     {1, 1, 0}, -- yellow
                     {1, 0, 1}, -- magenta
                     {0, 1, 1}, -- cyan
                     {1, 1, 1}  -- white
                    }

    return config
end

return config
