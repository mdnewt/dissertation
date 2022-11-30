function [HD,mx,mn,b,d] = hausdorff(mask1,mask2,vox,dir)
% HAUSDORFF.M is a 3D Hausdorff distance mapping function for binary masks.
% This function computes the boundary of each mask and computes a map of
% Hausdorff distances for each voxel in the boundary layer, and also
% generates maximum and mean Hausdorff distances.
% 
% Inputs: mask1 - the 'fixed' mask, supplied as a logical
%         mask2 - the 'moving' mask, supplied as a logical
%           vox - optional input denoting voxel size
%           dir - 

% Copyright 2019 Michael Newton. E-mail: mnewton@ltu.edu.

%% Distance Map-Based Computation of Hausdorff Distance 

% tic

% Create structuring element for 1-pixel erosion
S = false(3,3,3) ;
S(2,2,1  ) = true;
S(2,1:3,2) = true;
S(1:3,2,2) = true;
S(2,2,3  ) = true;

% Use erosion to obtain 1-pixel mask boundary
b1 =  mask1 & ~imerode(mask1,S);
b2 =  mask2 & ~imerode(mask2,S);
b  = {b1,b2}                   ;

% Create BW distance maps for boundaries
d1 =  bwdist(b1);
d2 =  bwdist(b2);
d  = {d1,d2}    ;

% Check distances of each boundary at the other boundary
HD1 =  d2; HD1(~b1) = 0;
HD2 =  d1; HD2(~b2) = 0;
HD  = {HD1,HD2}        ;

% Compute mean and max distances
mn1 =  mean(HD1(b1));
mn2 =  mean(HD2(b2));
mx1 =  max( HD1(: ));
mx2 =  max( HD2(: ));
mn  = [mn1,mn2]     ;
mx  = [mx1,mx2]     ;

if exist('vox','var'), mn = mn .* vox; mx = mx .* vox; end

if exist('dir','var') && isnumeric(dir) && any(dir==[1,2])
    HD =     HD{dir}     ;
    mn =     mn(dir)     ;
    mx =     mx(dir)     ;
    b  =      b{dir}     ;
    d  =  d{rem(dir,2)+1};
end

% toc

%% Brute Force Computation of Hausdorff Distance (Validation)
% 
% tic
% 
% % List of all boundary voxels
% [x1,y1,z1] = ind2sub(size(b1),find(b1));
% [x2,y2,z2] = ind2sub(size(b2),find(b2));
% 
% % Loop 1: 
% fED = @(x1,x2,y1,y2,z1,z2) sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
% ed1 = zeros(length(x1),1,'single');
% parfor ii = 1:length(x1)
%     temp = zeros(length(x2),1,'single');
%     for jj = 1:length(x2)
%         temp(jj) = fED(x1(ii),x2(jj),y1(ii),y2(jj),z1(ii),z2(jj));
%     end
%     ed1(ii) = min(temp);
% end
% 
% % Loop 2: 
% ed2 = zeros(length(x2),1,'single');
% parfor ii = 1:length(x2)
%     temp = zeros(length(x1),1,'single');
%     for jj = 1:length(x1)
%         temp(jj) = fED(x2(ii),x1(jj),y2(ii),y1(jj),z2(ii),z1(jj));
%     end
%     ed2(ii) = min(temp);
% end
% 
% ed_mn = {mean(ed1),mean(ed2)};
% ed_mx = {max( ed1),max( ed2)};
% 
% toc