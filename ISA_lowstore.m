% ISA_2 script uses the deformations created in ISA_1 and applies them to a
% set of tissue masks to create the average atlas.
% 
% To create the average atlas successfully, the following steps must be
% adhered to:
% 
% 1) Create a main directory and save ISA_1.m to this directory. Then,
%    create the following sub-directories within this directory:
% 
%   '0 - Pre-training Images' - OPTIONAL - In cases in which the user
%                               desires to resample the raw training images
%                               to a lower resolution prior to drawing the
%                               training masks, the raw training images can
%                               be placed here. The code ISA_0.M can then
%                               be used to create training images from this
%                               list of pre-training images.
% 
%   '1 - Training Images'     - This will contain .MATs of all image stacks
%                               to be used to create the average atlas. If
%                               the tissue of interest has a laterality (ie
%                               right and left femora), these should be
%                               demarcated by a '_L' or '_R' at the end of
%                               the filename so that proper mirroring can
%                               be performed.
%   
%   '2 - Training Masks'      - This will contain SUBFOLDERS for each
%                               tissue included in the average atlas. Each
%                               subfolder will then contain .MATs of the
%                               corresponding tissue for each training
%                               image.
% 
%   '3 - Average Atlas'       - This fill contain information related to
%                               the ISA process, including intermediate
%                               steps and all deformations used in the
%                               creation of the final average atlas.
% 
%   '4 - Unknown Samples'     - This will contain .MATs of the unknown
%                               images upon which the average atlas will be
%                               co-registered.
% 
%   '5 - Registrations'       - Registrations of the average atlas onto
%                               unknown samples will be saved here.
% 
% 2) Follow the steps required to run ISA_1.m and run the function.
% 
% 3) Populate the subfolders of '2 - Training Masks' with the desired
%    tissue masks. The folder names will be converted into variable names,
%    so it is important that the folder names represent valid variable
%    names. It will be assumed that corresponding masks will be present for
%    each training image in every subfolder.
% 
% 4) Run ISA_2.M.
% 
% This function is part of the iterative shape averaging (ISA) toolbox.
% Version history: V1 - 2017 May 02.

%% Create directory lists

% Training images
p = '1 - Training\Stacks\';
d = dir([p_fN,p,'*.mat']); d = dir2cell(d);

% Training masks
p2 = '1 - Training\Masks\';
d2 = dir2cell_isdir(dir(p2));

% Average atlas iterations
p3 = '2 - Average Atlas\Iterations\';
d3 = dir([p_fN,p3,'*.mat'])   ; d3 = dir2cell(d3);

%% Load in voxel information and reference image

% Determine voxel information
vox_train = load('vox.mat','vox_train');
vox_train = vox_train.vox_train;

% Determine reference image
if ~exist('ref.mat','file'), I_ref = 1;
else I_ref = load('ref.mat'); I_ref = I_ref.I_ref;
end

%% Transform and average masks

% Compile training masks
avg = struct;
for ii = 1:length(d2)
    %%
    % Compile true masks
    masks = cell(length(d),1);
    for jj = 1:length(d)
        mask = load([p2,d2{ii},'\',d{jj}],'roivol');
        
        %LEFT TO RIGHT MIRRORING - LMJ
        if(strfind(d{jj}, 'L'))
            mask.roivol = flipdim(mask.roivol, 2);
        end
        
        masks{jj} = single(mask.roivol);
    end
    clear mask

    % Perform affine transformation
    deforms = load('2 - Average Atlas\Iterations\it0.mat','deforms','reg');
    reg  = deforms.reg; deforms = deforms.deforms;
    Rfix = imref3d(size(masks{I_ref}),vox_train,vox_train,vox_train);
    I = find(~cellfun(@isempty,deforms(2:end,2)));
    for jj = 1:length(I)
        Rmov = imref3d(size(masks{I(jj)}),vox_train,vox_train,vox_train);
        masks{I(jj)} = imwarp(masks{I(jj)},Rmov,deforms{I(jj)+1,2},...
            'OutputView',Rfix);
    end

    % Perform sequential NRR registration
    for it = 1:length(d3)-1
        deforms = load(['2 - Average Atlas\Iterations\it',...
            num2str(it),'.mat'],'deforms');
        deforms = deforms.deforms;
        
        for jj = 1:length(d)
            
            % Warp mask
            masks{jj} = imwarp(masks{jj},deforms{jj+1,end});
        end
        it
    end
    
    % Average masks into an average roi
    avgatlas = logical(round(mean(cat(4,masks{:}),4)));
    
    % Create variable
    avg.(d2{ii}) = avgatlas;
    
    % Save atlas
    if ~exist('2 - Average Atlas\avg_atlas.mat','file');
        save('2 - Average Atlas\avg_atlas.mat','avg');
    else
        save('2 - Average Atlas\avg_atlas.mat','avg','-append');
    end
    
%     mpr3D(avgatlas,add_avgatlas);
end