% foldername = 'train_sun3ddiff_ft/';
% dataset = 'sun3d';
%train_sun3d %train_sun3ddiff  %train_sun3ddiff_ft
%train_sun3d_ftmpdiff_ft train_sun3d_ftmp_ft
function mpOverlapEval(foldername)
%foldername = 'train1diff/';
dataset = 'mp';
teslist = '../data/test';
dataRoot =  '/n/fs/rgbd/data/matterport/v1/';

fid = fopen(teslist);
line = fgetl(fid);
allsceneId ={};
while ischar(line)
    allsceneId{end+1}  = line;
    line = fgetl(fid);
end
    

nDCGScene = zeros(length(allsceneId),1);
for sid = 1:length(allsceneId)
    sceneId = allsceneId{sid};
    %sceneId = 'home_at/home_at_scan1_2013_jan_1';
    ind = find(sceneId=='/');
    if ~isempty(ind)
        sceneId_c = sceneId(ind(end)+1:end);
    else
        sceneId_c = sceneId;
    end
    overlapImgaeName = fullfile(dataRoot,sceneId,'overlaps',[sceneId_c '_iis_iou.png']);   
    IoUmap_flip = double(imread(overlapImgaeName))/1000;
    IoUmap = IoUmap_flip(end:-1:1,:);
    filename = sprintf('./checkpoints/%s/result/%sfeatAll.h5',foldername,strrep(sceneId,'/','--'));
    feat = hdf5read(filename,'featAll');
    feat = feat(:,1:size(IoUmap,1));
    feat = feat';
    distPredMap = pdist2(feat,feat);
    
    % read config file 
    if strcmp(dataset,'sun3d')
        configFileName = fullfile(dataRoot,sceneId,[sceneId_c '_ground_truth.conf']);   
        sample = 10;
    else
        configFileName = fullfile(dataRoot,sceneId,[sceneId_c '.conf']);   
        sample = 1;
    end
    data = readconfigFile(configFileName);
    % cacluate travel distance 
    centers = data.extrinsic(:,4,1:sample:end);
    centers = squeeze(centers);
    centers = centers(1:3,:);
    dis = centers(:,2:end) - centers(:,1:end-1);
    dis = sqrt(sum(dis.*dis));
    cdis = zeros(1,size(centers,2));
    for i=2:size(centers,2)
        cdis(i) = dis(i-1)+cdis(i-1);
    end
    
    [n,m] = ndgrid(1:size(centers,2),1:size(centers,2));
    disMatrix = abs(cdis(n)-cdis(m));
    
    dis_threhold = 0.5;
    
    %% for each qurey, caculate the Discounted cumulative gain
    a = 1:size(IoUmap,1);
    a = log(a)/log(2);
    a(1) = 1;
    nDCG = nan(1,size(IoUmap,1));
    for q = 1:size(IoUmap,1) 
        predict_dis = distPredMap(q,:);
        [~,predList] = sort(predict_dis);
        if strcmp(dataset,'sun3d')
           valid = find(disMatrix(q,:)>dis_threhold);
        else
           valid = find(IoUmap(q,:)<0.99); 
        end
        
        predList = predList(ismember(predList,valid));
        predict_dis = predict_dis(predList);

        predscores = IoUmap(q,predList);
        [~, goldList] = sort(IoUmap(q,:),'descend');
        goldList = goldList(ismember(goldList,valid));
        scores = IoUmap(q,goldList);
        DCGp = sum(predscores./a(1:length(scores)));
        DCGg = sum(scores./a(1:length(scores)));
        if DCGg>0
           nDCG(q) = DCGp/DCGg;
        end
    end
    nDCGScene(sid) = nanmean(nDCG);
end
nDCGScene
fprintf('mean nDCG:%f\n', mean(nDCGScene));
