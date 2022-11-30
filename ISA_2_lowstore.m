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
p = 'Stacks\';
d = dir([p,'*.mat']); d = dir2cell(d);
d = d(whichiswhich{folderN,1});

% Training masks
p2 = 'Masks\';
d2 = dir2cell_isdir(dir(p2));

% Average atlas iterations
p3 = [p_fN,'1 - Average Atlas\Iterations\'];
d3 = dir([p3,'*.mat']); d3 = dir2cell(d3);

%% Load in voxel information and reference image

% Determine voxel information
vox_train = load('vox.mat','vox_train');
vox_train = vox_train.vox_train;

% Determine reference image
if   ~exist('ref.mat','file'), I_ref = 1;
else, I_ref = load('ref.mat'); I_ref = I_ref.I_ref;
end

%% Transform and average masks

% Compile training masks
avg = struct;
for ii = 1:length(d2)
    %%
    
    t = tic;
    
    % Compile true masks
    masks = cell(length(d),1);
    for jj = 1:length(d)
        
        mask = load([p2,d2{ii},'\',d{jj}],'roivol');
        
        % Take up to 24 um
        mask.roifull = mask.roivol;
        mask.roivol = stackinterp(mask.roivol,...
            round(size(mask.roivol)./2),'interp3','Nearest');
        
        %LEFT TO RIGHT MIRRORING - LMJ
        if(strfind(d{jj}, 'L'))
            mask.roivol = flip(mask.roivol, 2);
        end
        
        masks{jj} = single(mask.roivol);

    end
    clear mask

    % Perform affine transformation
    deforms = load([p3,'it0.mat'],'deforms','reg');
    reg  = deforms.reg; deforms = deforms.deforms;
    Rfix = imref3d(size(masks{I_ref}),vox_train,vox_train,vox_train);
    I = find(~cellfun(@isempty,deforms(2:end,2)));
    for jj = 1:length(I)
        Rmov = imref3d(size(masks{I(jj)}),vox_train,vox_train,vox_train);
        masks{I(jj)} = imwarp(masks{I(jj)},Rmov,deforms{I(jj)+1,2},...
            'OutputView',Rfix,'interp','nearest');
    end
%     looptrack(folderN,ii,t2,[d{ii},' affine']);

    % Perform sequential NRR registration
    for it = 1:length(d3)-1
        deforms = load([p3,'it',...
            num2str(it),'.mat'],'deforms');
        deforms = deforms.deforms;
        
        for jj = 1:length(d)
            
            % Warp mask
            masks{jj} = imwarp(masks{jj},deforms{jj+1,end},'interp','nearest');
        end
%         it
    end
%     looptrack(folderN,ii,t3,[d{ii},' NRR']);
    
    % Average masks into an average roi
    avg_atlas = round(mean(cat(4,masks{:}),4));
    
    % Back to 12 um
    avg_atlas = round(stackinterp(avg_atlas,size(avg_atlas) .* 2,...
        'interp3','Nearest'));
    
    % Create variable
    avg.(d2{ii}) = avg_atlas;
    
    % Save atlas
    if ~exist([p_fN,'1 - Average Atlas\avg_atlas.mat'],'file')
        save( [p_fN,'1 - Average Atlas\avg_atlas.mat'],'avg','-v7.3');
    else
        save( [p_fN,'1 - Average Atlas\avg_atlas.mat'],'avg','-append');
    end
    
    looptrack(ii,length(d2),t,[d2{ii},' total completion time']);
%     mpr3D(avgatlas,add_avgatlas);
end

avg_image = load([p_fN,'1 - Average Atlas\avg_image.mat'],'avg_image');
avg_image = avg_image.avg_image;
avg_image = round(stackinterp(avg_image,size(avg_image) .* 2,...
        'interp3','Linear'));
avg.Image = avg_image;
save([p_fN,'1 - Average Atlas\avg_atlas.mat'],'avg','-append');