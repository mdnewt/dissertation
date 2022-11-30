function mask = purify(mask)
% A MATLAB implementation of the imageJ 'purify' function.

mx = bwconncomp(mask);
mx = cellfun(@numel,mx.PixelIdxList);
mx = max(mx);

mask = bwareaopen(mask,mx-1);