function [reg,regM,t,warp,regM_aff] = coregister2(mov,fix,vox,mov_roi)
% COREGISTER performs two-step coregistration of a moving volume onto a
% fixed volume. The first step is an affine registration to co-align the
% two volumes, and the second is a non-rigid registration to fine-tune the
% original alignment.
% 
% Input variables:
%     mov - preprocessed data for reference stack
%     fix - preprocessed data for unknown stack
%     vox - voxel sizes of the input images in the order {mov,fix}, where
%           mov and fix are 1x3 vectors
% mov_roi - optional cell of ROIs associated with mov. If supplied, each
%           ROI is registered onto fix.
% 
% This function is part of the iterative shape averaging (ISA) toolbox.
% Version history: V1 - 2017 May 02.

%% Create reference frames

% Check for proper vox input - currently only built for iso
if ~exist('vox','var') || isempty(vox), vox = {[1,1,1],[1,1,1]}; end

% Create reference frames
R_mov = imref3d(size(mov),vox{1}(1),vox{1}(2),vox{1}(3));
R_fix = imref3d(size(fix),vox{2}(1),vox{2}(2),vox{2}(3));

%% Affine registration on original inputs

t1 = tic;
[reg,TF] = affinereg(mov,fix,'affine',vox{1},vox{2});

if exist('mov_roi','var')
    L = length(mov_roi );
    regM_aff = cell(1,L);
    for ii = 1:L
        regM_aff{ii} = logical(imwarp(single(mov_roi{ii}),R_mov,TF,'nearest','OutputView',R_fix));
    end
end
t(1) = toc(t1);

%% NRR on top of rigid registration

[deform,reg] = imregdemons(reg,fix);
if exist('mov_roi','var')
	regM = cell(1,L);
    for ii = 1:L
        regM{ii} = logical(imwarp(single(regM_aff{ii}),deform,'nearest'));
    end
end
t(2) = toc(t1);
 
%% Create 'warp' output
 
warp = {TF,deform};