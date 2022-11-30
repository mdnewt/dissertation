function d = dir2cell(d,isDir)

if exist('isDir','var') && isDir
    d = d(3:end); d = cat(1,{d(cat(1,d(:).isdir)).name})';
else
    d = cat(1,{d(:).name})';
end

