% Setup
p_A       = 'Analysis\';

folderN = figgen{ii,1}              ;
p_fN    = [p_A,num2str(folderN),'\'];

% Load stuff
avg   = load([p_fN,'1 - Average Atlas\avg_atlas.mat']); avg = avg.avg;
stack = load(['Stacks\'    ,name,'.mat']); stack = stack.stack;
tru   = load(['Masks\Cart\',name,'.mat']);
epi   = load(['Masks\Epi\' ,name,'.mat']);
regM  = load([p_fN,'2 - Registrations\',name,'.mat'],'regM'); regM = regM.regM;

% Left to right mirroring
if(strfind(ukEpi(ukN).name, 'L'))
    stack      = flip(stack     ,2); epi.roivol = flip(epi.roivol,2);
    epi.bone   = flip(epi.bone  ,2); tru.roivol = flip(tru.roivol,2);
    tru.mask   = flip(tru.mask  ,2);
end

% Perform segmentation and compute Hausdorff distance map
if ii <= 10
    % Bone Segmentation
    if(strfind(ukEpi(ukN).name, '_')), bone = regM{2} & stack > ACLR2W_TH;
    else, bone = regM{2} & stack > stan_TH;
    end

    % Compute Hausdorff distance map
   [HD,mn,mx,b,d] = hausdorff(epi.bone,bone);
else
    % AC segmentation
    if(strfind(figgen{ii,2},'_')), bone2 = stack > ACLR2W_TH;
    else, bone2 = stack > stan_TH;
    end
    bone2     = imj_purify(bone2);
    bone_fill = imclose(bone2,strel('disk',22));
    bone_fill = imfill(bone_fill,'holes');
    vol       = logical(regM{1}) & ~bone_fill & stack >= 100;
    vol       = purify(vol);

    % Compute Hausdorff distance map
   [HD,mn,mx,b,d] = hausdorff(tru.mask,vol);
end

% Create patch
FV   = isosurface(vol);
cdat = interp3(d{1},FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));
figure,P = patch(FV,'EdgeColor','none','FaceColor','interp','CData',cdat);
axis equal

% Save
savefig(['S:\OrthopaedicsProjects\Newton\Automated Segmentation\',...
         'CE-uCT Knee - Atlas Registration\Validation\Manuscript',...
         '\Figure\',fg_names{ii},'.fig']);