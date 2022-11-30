% ISA_1 script creates an average image from a set of input images. ISA 
% refers to Iterative Shape Averaging, a process described by Rohlfing et
% al (2004, doi:10.1016/j.neuroimage.2003.11.010) - see also "Quo Vadus, 
% Atlas-Based Segmentation?" by Rohlfing et al from The Handbook of Medical
% Image Analysis: Segmentation and Registration Methods.
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
% 2) Create a .mat save file called 'vox.mat' with 2 variables saved to it,
%    and place this file in the main directory. These are the 2 variables:
% 
%       vox_train - The voxel size of the training image stacks and masks.
%                   The code assumes that all tissue masks are stored at 
%                   the same resolution. The training data should be stored
%                   at the same resolution that the average atlas will be
%                   created and used at.
% 
%       vox_unk   - The voxel size of the unknown samples. If
%                   co-registration is performed at a different resolution,
%                   the resultant atlas will be resampled to match the
%                   original resolution of the unknown image.
% 
% 3) Populate folder 1. This can be done directly by saving the raw image 
%    stacks to folder 1. Alternately, the raw images can be saved to 
%    folder 0, and ISA_0.m can be used to resample these and save them to 
%    folder 1.
% 
% 4) Designate a reference image by saving the variable I_ref to the file
%    'ref.mat', where I_ref refers to the index of the training image in
%    folder 1 that should be used as reference image. The reference image
%    is the image that the others are registered on top of during the
%    initial affine transformation to start the ISA process. Ideally this
%    is a sample of excellent alignment and representative geometry. If no
%    I_ref is specified, the first image in the folder will be used.
% 
% 5) Run ISA_1.M.
% 
% This function is part of the iterative shape averaging (ISA) toolbox.
% Version history: V1 - 2017 May 02.

%% Create directory list of training images

%% '1-Training Images' changed to '1 - Training\Stacks\'

% Generate directory
p = 'Stacks\';
d = dir([p,'*.mat']); d = cat(1,{d(:).name})';
d = d(whichiswhich{folderN,1});
f = @(x) x(1:end-4);
name = cellfun(f,d,'UniformOutput',false); clear f

%% Create cell to hold all deformations for later application to atlas

deforms = cell(length(d)+1,2);
deforms{1,1} = 'Samples'; deforms(2:end,1) = d;
deforms{1,2} = 'Affine';

%% Load in voxel information and reference image

% **MDN - commented out since you define vox_train in master.m
% % Determine voxel information
% vox_train = load('vox.mat','vox_train');
% vox_train = vox_train.vox_train;

% Determine reference image
% if ~exist('ref.mat','file'), I_ref = 1;
% else I_ref = load('ref.mat'); I_ref = I_ref.I_ref;
% end

% Assign firstR as reference
firstR = strfind(d,'R');
firstR = cellfun(@(x)isequal(x, 6), firstR);
if(find(firstR))
    [row] = find(firstR);
    I_ref = row(1);
else
    I_ref = 1;
end

%% Step 1 - Perform affine registration to roughly co-align

% Load reference image
fix   = load([p,d{I_ref}]);

% **MDN change size temporarily to 24 um
fix.stackfull = fix.stack;
fix.stack = stackinterp(fix.stack,round(size(fix.stack)./2));

% LEFT TO RIGHT MIRRORING - LMJ
if(strfind(d{I_ref}, 'L'))
    fix.stack = flip(fix.stack, 2);
end
fix.R = imref3d(size(fix.stack),vox_train,vox_train,vox_train);

% Create cell for registrations
reg = cell(length(d),1);
reg{I_ref} = fix.stack;

% Register the other stacks onto the reference
I = find(1:length(d) ~= I_ref);
for ii = 1:length(I)
    %%
    
    t = tic;
    
    % Load moving sample in
    mov   = load([p,d{I(ii)}],'stack');
    
    % **MDN change size temporarily to 24 um
    mov.stackfull = mov.stack;
    mov.stack = stackinterp(mov.stack,round(size(mov.stack)./2));

    % LEFT TO RIGHT MIRRORING - LMJ
    if(strfind(d{I(ii)}, 'L'))
        mov.stack = flip(mov.stack, 2);
    end
   
    mov.R = imref3d(size(mov.stack),vox_train,vox_train,vox_train);
    
    % Configure and perform registration
   [optimizer,metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = 500;
    geomtform = imregtform(mov.stack,mov.R,fix.stack,fix.R,...
        'affine',optimizer,metric);
    reg{I(ii)} = imwarp(mov.stack,mov.R,geomtform,'linear','OutputView',...
        fix.R);
    
    deforms{I(ii)+1,2} = geomtform;
    
    looptrack(ii,length(I),t,'Rigid');
end

% Create 'iteration 0' average atlas
avg_image = {mean(cat(4,reg{:}),4)};

% Save iteration 0
%changed '3' to '2'
save([p_fN,'1 - Average Atlas\Iterations\it0.mat'],'reg','deforms',...
    'avg_image', '-v7.3');

fprintf('Iteration 0 complete.\n');

% Optional visualization
% for ii = 1:(length(d) - 1)
%     figure,imshowpair(reg{I_ref}(:,:,100),reg{I(ii)}(:,:,100));
%     title(d{I(ii)});
%     w = waitforbuttonpress;
%     close(gcf);
% end

%% Step 2 - Iterative rounds of NRR

% Run loop for different iterations - the loop will run a minimum of 2
% iterations of NRR, stopping when the average deformation reaches below
% 1/4 voxel, up to a maximum of 10 iterations
it = 1;
%%
while (it < 3 || diff(avgdef(end:-1:end-1)) > .1) && it < 10
    %%
    
    % Registrations from previous step
    reg_prev = reg;
    
    % Setup new row of deforms
    deforms{1,2} = ['NRR It. ',num2str(it)];
    deforms(2:end,2) = {[]};
    
    % Reference image
    fix.stack = avg_image{it};
    
    % Register all stacks onto the reference
    for ii = 1:length(d)
        
        t = tic;
        
        % Load moving sample in from last iteration
        mov.stack = reg_prev{ii};
        
        % Perform NRR
       [deform,reg{ii}] = imregdemons(mov.stack,fix.stack);
        deforms{ii+1,2} = deform;
        
        looptrack(ii,length(d),t,['It. ',num2str(it)]);

    end
    
    % Create average image
    avg_image{it+1} = mean(cat(4,reg{:}),4);
    
    % Compute average deformation
    ad = sqrt(sum((mean(cat(5,deforms{2:end,2}),5)).^2,4));
    avgdef(it) = mean(ad(:));
    
    % Save iteration
    save([p_fN,'1 - Average Atlas\Iterations\it',num2str(it),'.mat'],...
        'reg','deforms','avg_image','avgdef','-v7.3');
    
%     % Optional visualization
%     for jj = 1:(length(d_stack) - 1)
%         figure,imshowpair(avgatlas{it}(:,:,100),reg{jj}(:,:,100));
%         title(d_stack{I(jj)});
%         w = waitforbuttonpress;
%         close(gcf);
%     end

    fprintf(['NRR iteration ',num2str(it),' complete.\n']);
    it = it+1;
    
end

% Save final average image
avg_image = avg_image{end};
save([p_fN,'1 - Average Atlas\avg_image.mat'],'avg_image','-v7.3');