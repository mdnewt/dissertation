%% Setup
pA        = 'Analysis\';
pM        = 'Masks\'   ;
vox_train = 0.012      ;
stan_TH   = 3800       ;
ACLR2W_TH = 5000       ;

%% Loop
for ii = 31:50
    %%
    
    % Check for the presence of registrations
    pFN = [pA,num2str(ii),'\'];
    d   = dir2cell(dir([pFN,'2 - Registrations\*.mat']));
    L   = length(d);
    
    % Prepare CSV header
    output = {'Unknown','AC Atten (HU)','Mean AC HD (µm)',...
        'Max AC HD (µm)','Mean Bone HD (µm)','Max Bone HD (µm)',...
        'Mean Fill HD (µm)','Max Fill HD (µm)'};
    
    %%
    if L==0, continue
    else
        %% If present, run loop
        
        avg = load([pFN,'\1 - Average Atlas\avg_atlas.mat']);
        avg = avg.avg;
        
        % Create directory
        if ~exist([pFN,'3 - Extras'],'dir'), mkdir([pFN,'3 - Extras']); end
        
        %% Inner Loop
        for jj = 1:L
            %%
            
            t = tic;
            
            % Load stuff
            stack = load(['Stacks\',d{jj}]); stack = stack.stack;
            tru   = load(['Masks\Cart\',d{jj}]);
            epi   = load(['Masks\Epi\' ,d{jj}]);
            regM  = load([pFN,'2 - Registrations\',d{jj}]);
            reg   = regM.reg; regM = regM.regM;

            % Left to right mirroring
            if(strfind(d{jj}, 'L'))
                stack      = flip(stack     ,2); epi.roivol = flip(epi.roivol,2);
                epi.bone   = flip(epi.bone  ,2); tru.roivol = flip(tru.roivol,2);
                tru.mask   = flip(tru.mask  ,2);
            end
            
            % Bone Segmentation
            if(strfind(d{jj}, '_')), bone = regM{2} & stack > ACLR2W_TH;
            else, bone = regM{2} & stack > stan_TH;
            end

            % Compute Hausdorff distance map
           [HD,mx,mn,b,dist] = hausdorff(epi.bone,bone,vox_train*1000,1);
            HD_bone = struct('HD',HD,'mn',mn,'mx',mx,'b',b,'dist',dist);

            % AC segmentation
            if(strfind(d{jj},'_')), bone2 = stack > ACLR2W_TH;
            else, bone2 = stack > stan_TH;
            end
            bone2     = purify(bone2);
            bone_fill = imclose(bone2,strel('disk',22));
            bone_fill = imfill(bone_fill,'holes');
            vol       = logical(regM{1}) & ~bone_fill & stack >= 100;
            vol       = purify(vol);
            
            % Compute Hausdorff distance map
           [HD,mx,mn,b,dist] = hausdorff(tru.mask,vol,vox_train*1000,1);
            HD_AC = struct('HD',HD,'mn',mn,'mx',mx,'b',b,'dist',dist);
            
            % Reg Bone Fill
            if(strfind(d{jj},'_')), bone2 = reg > ACLR2W_TH;
            else, bone2 = reg > stan_TH;
            end
            bone2    = imj_purify(bone2);
            reg_fill = imclose(bone2,strel('disk',22));
            reg_fill = imfill(reg_fill,'holes');
            
            % Compute Hausdorff distance map
           [HD,mx,mn,b,dist] = hausdorff(bone_fill,reg_fill,vox_train*1000,1);
            HD_fill = struct('HD',HD,'mn',mn,'mx',mx,'b',b,'dist',dist);
            
            % Compute AC attenuation
            atten = mean(stack(vol));
            
            % Place row of output
            output = [output;{d{jj},atten,HD_AC.mn,HD_AC.mx,HD_bone.mn,...
                      HD_bone.mx,HD_fill.mn,HD_fill.mx}];
            
            % Organize
            cart      = struct('mov',vol     ,'fix',tru      );
            bone      = struct('mov',bone    ,'fix',epi.bone );
            bone_fill = struct('mov',reg_fill,'fix',bone_fill);
            
            % Save extras
            save([pFN,'3 - Extras\',d{jj}],'HD_AC','HD_bone','HD_fill',...
                'atten','cart','bone','bone_fill','-v7.3');
            
            looptrack(jj,L,t,['It. ',num2str(ii),': ',d{jj}]);
            
        end
    end
    
    % Save output
    if ~exist('Excel Output\Extras','dir'), mkdir('Excel Output\Extras'); end
    xlswrite('Excel Output\Extras\extras.xls',output,num2str(ii));
    
end