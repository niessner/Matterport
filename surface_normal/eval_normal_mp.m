function [ total_error, angle_error ] = eval_normal_mp( folder, root_path )
%EVAL_NORMAL_MP Summary of this function goes here
%   Detailed explanation goes here
angle_error = [];
image_list = dir([folder '/*_normal_est.png']);
for a = 1:length(image_list)
    fprintf('Eval: %s\n', image_list(a).name);
    est_map = im2double(imread([folder image_list(a).name]));
    try
        p = strfind(image_list(a).name, '_');
        gnd_path_root = [root_path image_list(a).name(1:p(1)-1) '/undistorted_normal_images/' image_list(a).name(p(1)+1:p(4)-1)];
        nx = im2double(imread([gnd_path_root '_nx.png']));
        ny = im2double(imread([gnd_path_root '_ny.png']));
        nz = im2double(imread([gnd_path_root '_nz.png']));
        val_map = (nx.^2+ny.^2+nz.^2)>0.5;
        gnd_map = cat(3, nx, 1-nz, ny);
    
    catch
        continue;
    end
    
    
    val_map = val_map(1:2:end,1:2:end);
    gnd_map = gnd_map(1:2:end,1:2:end,:);
    
    est_norm = est_map * 2 - 1;
    gnd_norm = gnd_map * 2 - 1;
    
    est_norm = est_norm ./ repmat( sqrt(sum(est_norm.^2,3)), [1 1 3] );
    gnd_norm = gnd_norm ./ repmat( sqrt(sum(gnd_norm.^2,3)), [1 1 3] );
    
    error = acos(dot(reshape(est_norm,[],3), reshape(gnd_norm,[],3), 2));
    angle_error(a).error = error(val_map(:));
    
end
total_error = vertcat(angle_error.error);
fprintf('Mean: %4.2f\n', mean(total_error)*180/pi);

end

