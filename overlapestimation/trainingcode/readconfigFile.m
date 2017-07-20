function data = readconfigFile(configFileName)
fid =  fopen(configFileName,'r');
tline = fgetl(fid);
imagename={};
depthname={};
extrinsic =[];
while ischar(tline)
    str_s = strsplit(tline);
    if ~isempty(str_s)
        if strcmp(str_s{1},'scan')
            imagename{end+1} = str_s{3};
            depthname{end+1} = str_s{2};
            extrinsic = cat(3,extrinsic,reshape(str2double(str_s(4:end)),[4,4])');
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
data.imagename = imagename;
data.depthname = depthname;
data.extrinsic = extrinsic;


end