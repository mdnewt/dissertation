
%   This is the final script which traverses the folders created by setup.m
%   expects that folders in the base directory are set up in the following way
%   only folders present are 'Stacks', 'Masks', and folders numbered (1 to n)
%   where n is total number of atlases generated
%   otherwise vox.mat is the only file in the directory

% Initialize vox_train (all testing done on same res images)
% vox_train = 0.024;
vox_train = 0.012;
stan_TH   = 3800; %stanphd & aclr2wk threshold
ACLR2W_TH = 5000;
n         = 40;
p_A       = 'Analysis\';
load('rando.mat');

%% Loop 1: ISA 1
for folderN = 1:50
    %%
%     %navigate to correct folder
%     cd([p,num2str(folderN)]);
    
    t_ISA1 = tic;

    % MDN - above code errors bc changing the current directory navigates
    % away from the ISA_1LJ function. Revised so that the analysis
    % directory 'folderN' is supplied to the function and the cd remains
    % the same.
    p_fN = [p_A,num2str(folderN),'\'];
    
    if folderN > 10
        % Run ISA_1
        ISA_1_lowstore;
    else
        % Training images
        p = 'Stacks\';
        d = dir([p,'*.mat']); d = dir2cell(d);
        d = d(whichiswhich{folderN,1});

        % Training masks
        p2 = 'Masks\';
        d2 = dir2cell_isdir(dir(p2));
        
        avg = struct;
        avg.Image = load([p,d{1}]);
        avg.Image = avg.Image.stack;
        for ii = 1:length(d2)
            avg.(d2{ii}) = load([p2,d2{ii},'\',d{1}],'roivol');
            avg.(d2{ii}) = avg.(d2{ii}).roivol;
        end
        save([p_fN,'1 - Average Atlas\avg_atlas.mat'],'avg','-v7.3');
    end
    
    looptrack(folderN,n,t_ISA1,'ISA 1 complete.');
    
end

%% Loop 2: ISA 2
for folderN = 1:30
    %%
    
    if folderN > 10
        % MDN - See above
        p_fN = [p_A,num2str(folderN),'\'];

        % Run ISA_2
        t_ISA2 = tic;
        ISA_2_lowstore;
        looptrack(folderN,n,t_ISA2,'ISA 2 complete.');
    end
        
end
    
%% Loop to delete iterations
for folderN = 1:n
    %%
    %delete iterations after ISA_2 to save space
    delete '2 - Average Atlas\Iterations\*.mat'
end

%% Loop 3: Unknown registration and data computation
%mkdir 'Excel Output';
for folderN = 1:50
    %%
    t = tic;
    
    p_fN = [p_A,num2str(folderN),'\'];
    
    %Initializes directories
    ts = whichiswhich{folderN,1};
    uk = find(~ismember(1:50,whichiswhich{folderN,1}))';
    Stacks   = dir(fullfile('Stacks\','*.mat'));
    ukStacks = Stacks(uk);
    Stacks   = Stacks(ts);
    ukCart   = dir(fullfile('Masks\Cart\','*.mat')); ukCart = ukCart(uk);
    ukEpi    = dir(fullfile('Masks\Epi\' ,'*.mat')); ukEpi  = ukEpi( uk);
    
    %Initialize XLS and Records on first column the list of samples used 
    %to construct the atlas (each atlas will have its own xls)
    % MDN - Put all excel documents in new 'Excel Output' subdirectory
    xlswrite(['Excel Output\',num2str(folderN) '.xlsx'], {Stacks.name});
    
    %write output measure headings to xls
    table = {'Unknowns', 'sens_E_TV', 'spec_E_TV', 'DSC_E_TV', ...
        'sens_E_BVTV', 'spec_E_BVTV', 'DSC_E_BVTV', ...
        'bvtv', 'bvtv_R', 'tmd', 'tmd_R', 'bmd', 'bmd_R', ...
        'tbth', 'tbth_R', 'tbthstd', 'tbthstd_R', 'tbn', ...
        'tbn_R', 'tbsp', 'tbsp_R', 'SMI', 'SMI_R', 'sens_C', ...
        'spec_C', 'DSC_C', 'cartVol', 'cartVol_R', 'cartTh', ...
        'cartTh_R', 'cartThSt', 'cartThSt_R', 'bvtv_PE', 'tmd_PE', ...
        'bmd_PE', 'tbth_PE', 'tbthstd_PE', 'tbn_PE', 'tbsp_PE', ...
        'SMI_PE', 'cartVol_PE', 'cartTh_PE', 'cartThSt_PE'};
    xlswrite(['Excel Output\',num2str(folderN) '.xlsx'], table,  1, 'A2');
    
    %load average masks and stack
    avg = load([p_fN,'1 - Average Atlas\avg_atlas.mat']); avg = avg.avg;
    % Because sometimes avg image didn't save correctly in avg_atlas.mat
    if ~isfield(avg,'Image')
        avg.Image = load([p_fN,'1 - Average Atlas\avg_image.mat']);
        avg.Image = avg.Image.avg_image;
        avg.Image = stackinterp(avg.Image,size(avg.Image) .* 2,'interp3','Linear');
        save([p_fN,'1 - Average Atlas\avg_atlas.mat'],'avg','-v7.3');
    end
    
%     load([p_fN,'2 - Average Atlas\avg_image.mat']);
    
    %loops through unknowns
    numUnknowns = length(ukStacks);
    
    % MDN - added in case not all unknowns are computed
    nUr = randperm(numUnknowns,10);
    
    % Loop of each individual unknown sample
    %%
    for ukN = 1:numUnknowns
        %%
        t2 = tic;
        % Write current unknown ID to XLS 
        %LMJ - erased nUr  to read from all unknowns
        xlswrite(['Excel Output\',num2str(folderN) '.xlsx'], {ukStacks(ukN).name}, 1, ['A', num2str(ukN+2)]);
        
        % Load unknown stack - variable is 'stack'
        % Load random unknown
        stack = load(fullfile(ukStacks(ukN).folder, ukStacks(ukN).name));
%         stack = load(fullfile(ukStacks(ukN).folder, ukStacks(ukN).name));
        stack = stack.stack;
        
        %LEFT TO RIGHT MIRRORING - LMJ
        if(strfind(ukStacks(ukN).name, 'L'))
            stack = flip(stack, 2);
        end
        
        % Registration 
        % Register avg_image onto unknown stack and apply x-form on avg_masks 
        % coregister script modified into coregister2.m to apply x-form onto
        % both avg masks (epi and cart)
        
        % MDN - Bring down to 24 um *had to update stack interp*
        stack_24       = stackinterp(stack,round(size(stack)./2));
          avg_24       = avg;
          avg_24.Image = stackinterp(avg.Image,round(size(avg.Image)./2),'interp3','linear' );
          avg_24.Cart  = stackinterp(avg.Cart ,round(size(avg.Cart )./2),'interp3','nearest');
          avg_24.Epi   = stackinterp(avg.Epi  ,round(size(avg.Epi  )./2),'interp3','nearest');
        
       [reg_24,regM_24,~,warp,regM_aff] = coregister2(avg_24.Image,...
           stack_24,[],{avg_24.Cart,avg_24.Epi});
        
        % MDN - Resize back up to 12 um *had to update stack interp*
        reg     = stackinterp(reg_24    ,size(stack),'interp3','linear' );
        regM{1} = stackinterp(regM_24{1},size(stack),'interp3','nearest');
        regM{2} = stackinterp(regM_24{2},size(stack),'interp3','nearest');
       
        % Save intermediate files
        % regM{1} is the morphed AC mask
        % regM{2} is the morphed epi mask
        save(fullfile(p_fN,'2 - Registrations',ukStacks(ukN).name),'reg',...
            'regM','t','warp','regM_aff','-v7.3');
      
        % Processes epi mask
        % variables: bone is a trabecular bone volume mask, roivol is a
        % total trabecular volume mask
        
        epi = load(fullfile(ukEpi(ukN).folder,ukEpi(ukN).name));
        
        %LEFT TO RIGHT MIRRORING - LMJ
        if(strfind(ukEpi(ukN).name, 'L'))
            epi.roivol = flip(epi.roivol,2);
            epi.bone   = flip(epi.bone  ,2);
        end
        
       [sens_E_TV,spec_E_TV,DSC_E_TV] = senSpec(regM{2},epi.roivol);
       
        %TH MATCH - LMJ
        if(strfind(ukEpi(ukN).name, '_'))
            regbone = regM{2} & stack > ACLR2W_TH;
        else
            regbone = regM{2} & stack > stan_TH;
        end
        
       [sens_E_BVTV,spec_E_BVTV,DSC_E_BVTV] = senSpec(regbone,epi.bone);  
        
        % MDN - added this code inline in the above TH MATCH section
%         regbone(~regM{2}) = 0; %remove bone that's not in morphed Epi mask
        
        %morphed epi mask stats
        bvtv_R=numel(find(regbone==1))/numel(find(regM{2}==1)); % Whole BVTV
        tmd_R = mean(stack(regbone==1)); %TMD
        bmd_R = mean(stack(regM{2}==1)); %BMD
       [tbth_R,tbthstd_R] = boneJ(regbone,'1',vox_train); %Tb.Th.
        tbn_R = bvtv_R/tbth_R; %Tb.N.
        T = uint8(regbone);
        T(~regM{2}==1)=2; % This makes a uint8 matrix where background can be removed using a "0/0" threshold in Miji
       [tbsp_R,tbspstd_R] = boneJSp(T,'1',vox_train); %Tb.Sp.
        SMI_R = boneJSMI(regbone,1); %SMI
        
        %true mask stats
        % MDN - should rewrite so this isn't redone every time
        bvtv=numel(find(epi.bone==1))/numel(find(epi.roivol==1)); % Whole BVTV
        tmd = mean(stack(epi.bone==1)); %TMD
        bmd = mean(stack(epi.roivol==1)); %BMD
        [tbth,tbthstd] = boneJ(epi.bone,'1',vox_train); %Tb.Th.
        tbn = bvtv/tbth; %Tb.N.
        T = uint8(epi.bone);
        T(~epi.roivol==1)=2; % This makes a uint8 matrix where background can be removed using a "0/0" threshold in Miji
       [tbsp,tbspstd] = boneJSp(T,'1',vox_train); %Tb.Sp.
        SMI = boneJSMI(epi.bone,1); %SMI
        
        % Processes cartilage mask
        
        % load true trimmed cartilage mask - variable is 'mask' 
        % MDN - loaded in as struct called 'tru'
        tru = load(fullfile(ukCart(ukN).folder, ukCart(ukN).name));
        
        % LEFT TO RIGHT MIRRORING - LMJ
        if(strfind(ukCart(ukN).name, 'L'))
            tru.roivol = flip(tru.roivol, 2);
            tru.mask   = flip(tru.mask  , 2);
        end
        
        % TH MATCH - LMJ
        % MDN - this is a bone threshold so naming the mask 'cart' was
        % confusing - changed to bone2 for clarity
        if(strfind(ukCart(ukN).name, '_'))
            bone2 = stack > ACLR2W_TH;
        else
            bone2 = stack > stan_TH;
        end
        
        % trim cart roivol regM.Cart into true_cart_vol
        bone2     = imj_purify(bone2);
        bone_fill = imclose(bone2,strel('disk',22)); %22 size of strel used to generate true masks 
        bone_fill = imfill(bone_fill,'holes');
        vol       = logical(regM{1}) & ~bone_fill & stack >= 100;
        
%         vol = roivol;
%         vol(bone_fill > 0) = 0;
%         vol(stack < 100) = 0;
        
        % senSpec trimmed cartilage volume onto true trimmed cartilage volume.
       [sens_C,spec_C,DSC_C] = senSpec(vol,tru.mask);
        cartVol_R = sum(vol(:)>0)      ;
        cartVol   = sum(tru.mask(:)>0);
       [cartTh  ,cartThSt  ] = boneJ(tru.mask,'0',vox_train);
       [cartTh_R,cartThSt_R] = boneJ(vol      ,'0',vox_train);
        
        % write results to xls on next row
        table = {sens_E_TV, spec_E_TV, DSC_E_TV, ...
           sens_E_BVTV, spec_E_BVTV, DSC_E_BVTV ...
           bvtv, bvtv_R, tmd, tmd_R, bmd, bmd_R, tbth, tbth_R, ...
           tbthstd, tbthstd_R, tbn, tbn_R, tbsp, tbsp_R, SMI, ...
           SMI_R, sens_C, spec_C, DSC_C, cartVol, cartVol_R, cartTh, ...
           cartTh_R, cartThSt, cartThSt_R, PE(bvtv_R, bvtv), PE(tmd_R, tmd), ...
           PE(bmd_R, bmd), PE(tbth_R, tbth), PE(tbthstd_R, tbthstd), ...
           PE(tbn_R, tbn), PE(tbsp_R, tbsp), PE(SMI_R, SMI), ...
           PE(cartVol_R, cartVol), PE(cartTh_R, cartTh), ... 
           PE(cartThSt_R, cartThSt)};
        
        xlswrite(['Excel Output\',num2str(folderN) '.xlsx'], table, 1, ['B', num2str(ukN+2)]);
        
        looptrack(ukN,numUnknowns,t2,'Unknown Registration');
        
    end
    looptrack(folderN,40,t,['folderN = ',num2str(folderN),' complete.']);
end