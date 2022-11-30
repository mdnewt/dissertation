 %% 0. Define directories, file locations, and in-line functions
 
 % Define directories and file locations
 xlsx    = [pwd,'..\..\..\Scan Log.xlsx'];
 p       =  struct  ; p.nir    = '..\NIR\'; p.vct = '..\VivaCT Alignment\';
 p.ceuct = '..\CE-µCT\'; p.output = '1 - Images\';
 
 % Load and preprocess crop and scale files
 crop  = load('crop.mat'); crop = crop.crop;
 c     = round([crop.holder(2),crop.holder(2)+crop.holder(4),...
                crop.holder(1),crop.holder(1)+crop.holder(3)]);
%  scale = load('scale.mat','scalingfactor','rN','rV');
 projectionTF = load([pwd,'\..\..\..\Aim 2\Projection Phantom\projectionTF.mat']);
 
 % Define In-Line Functions
 se  = strel('disk',7); se2 = strel('disk',21);
 f     = @(x)     x(1:end-4);
 f2    = @(x)     x(1:2    );
 f3    = @(x)     x(4:7    );
 fCROP = @(x)     x(c(1):c(2),c(3):c(4))             ;
 fPROJ = @(x)     sqrt(sum(x.^2,3))                  ;
 fOC   = @(x)     imopen(x,se).*imclose(x,se)        ;
 fER   = @(x)     imreconstruct(imerode(x,se2),x)    ;
 fNORM = @(x)    (x-min(x(:)))./(max(x(:))-min(x(:)));
 fMM   = @(x,xp) [min([x;xp]),max([x;xp])]           ;
[opts,met] = imregconfig('multimodal');
 
 %% 1. Parse XLS file to extract image identifiers
[n1,t1] = xlsread(xlsx,'Endpoint Log'  ,'B1:J7' );
[n2,t2] = xlsread(xlsx,'Endpoint Log'  ,'O1:R26');
[n3,t3] = xlsread(xlsx,'CE-µCT Log'    ,'A1:E25');
 id = struct;
 id.nir        = num2cell(n1(~any(isnan(n1),2),[1:2,4:5,7:8]))    ;
 id.ep         = num2cell(n2(~any(isnan(n2),2), :           ))    ;
 id.vct        = num2cell(n1(~any(isnan(n1),2),[3,6,9]      ))    ;
 id.ceuct      = num2cell(n3(~any(isnan(n3),2), :           ))    ;
 id.nir        = cellfun(@num2str,id.nir  ,'UniformOutput',false) ;
 id.ep         = cellfun(@num2str,id.ep   ,'UniformOutput',false) ;
 id.vct        = cellfun(@num2str,id.vct  ,'UniformOutput',false) ;
 id.ceuct      = cellfun(@num2str,id.ceuct,'UniformOutput',false) ;
 id.name.nir   = t1(:,[1:2,4:5,7:8])                              ;
 id.name.ep    = t2                                               ;
 id.name.vct   = t1(1,[1,4,7]            )                        ;
 id.name.ceuct = t3                                               ;
 id.list.vct   = id.vct(  :    ); id.list.vct   = id.list.vct(  :);
 id.list.ceuct = id.ceuct(:,2:5); id.list.ceuct = id.list.ceuct(:);
 clear n1 n2 n3 t1 t2 t3
 
 %% 2. Generate a full catalog of available scans
 d           = struct                   ;
 d.nir       = dir2cell(dir(p.nir   ),1);
 d.vct       = dir2cell(dir(p.vct   ),1);
 d.ceuct     = dir2cell(dir(p.ceuct ),1);
 d.out.nir   = dir2cell(dir([p.output,p.nir  ,'*.mat']),1);
 d.out.vct   = dir2cell(dir([p.output,p.vct  ,'*.mat']),1);
 d.out.ceuct = dir2cell(dir([p.output,p.ceuct,'*.mat']),1);
 d.out.nir   = unique(cellfun(f,d.out.nir  ,'UniformOutput',false));
 d.out.vct   = unique(cellfun(f,d.out.vct  ,'UniformOutput',false));
 d.out.ceuct = unique(cellfun(f,d.out.ceuct,'UniformOutput',false));
 
 %% 3. Extract information about already-completed scans
 sl = load('subjectlist.mat','scanlist');
 scanlist = sl.scanlist; clear sl
 L = numel(scanlist); status = false(L,3);
 for ii = 1:L
     %%
     fn = [p.output,scanlist{ii},'.mat'];
     if exist(fn,'file')
         w = whos('-file',fn); w = cat(1,{w(:).name})'   ;
         status(ii,:) = ismember({'nir','vct','ceuct'},w);
     end
 end
 
 %% 4. Process NIR and VivaCT images by endpoint
 
 % Pull a list of missing NIR and VivaCT for which scan IDs are available
 tbd = scanlist(any(~status(:,1:2),2));
 tbd = unique(cellfun(f2,tbd,'UniformOutput',false));
 tbd = tbd(ismember(tbd,id.ep(:,1)));
 
 % Which endpoints do these scans come from?
 eps = ismember(id.ep(:,1),tbd); eps = cellfun(@str2num,unique(id.ep(eps,2)))';
 L  = length(eps); numproc = 1 ;
 eps = 1:5; L = 5; numproc = 1;
%  eps = 1; L = 1; numproc = 1; ii = 1; % Manual override for first endpoint
 
 %%
 for ii = eps
     %%
     if isempty(ii), continue; end % Loop isn't smart enough not to run
     t = tic;
     
     % Only proceed if full scan data is available for this endpoint
     idcheck = cellfun(f3,d.nir,'UniformOutput',false);
     if any(~contains(id.nir(ii,:),idcheck)), continue; end
     if any(~contains(id.vct(ii,:),d.vct  )), continue; end
%%
     %% Load in NIR scan data
     NIR = struct;
     for jj = 1:length(id.nir)
         %%         
         % Determine field name
         if ~isempty(id.name.nir{1,jj})
             fname1 = genvarname(id.name.nir{1,jj});
         end
         name2 = genvarname(id.name.nir{2,jj}); name2 = name2(4:end);
         
         % Import NIR scans in nested structure
         sn = contains(d.nir,id.nir(ii,jj));
        [NIR.im.(fname1).(name2).x700,...
         NIR.im.(fname1).(name2).x800,...
         NIR.im.(fname1).(name2).W]   = NIR_batch(p.nir,d.nir{sn})        ;
     
         NIR.im.(fname1).(name2).x700 = fCROP(NIR.im.(fname1).(name2).x700);
         NIR.im.(fname1).(name2).x800 = fCROP(NIR.im.(fname1).(name2).x800);
         NIR.im.(fname1).(name2).W    = fCROP(NIR.im.(fname1).(name2).W   );
         
         % Apply median filter to white light image
         NIR.im.(fname1).(name2).W = medfilt2(NIR.im.(fname1).(name2).W,[5,5]);
     end
     NIR.info = imfinfo([p.nir,d.nir{sn},'\',d.nir{sn},'_White.TIF']);
     NIR.pix  = 10/NIR.info(1).XResolution; pix = NIR.pix; % Convert to mm/pix
     NIR.r2d  = imref2d(size(NIR.im.(fname1).(name2).W),pix,pix);

     %% Load in VivaCT scan data
     VivaCT = struct;
     for jj = 1:size(id.vct,2)
         %%
         % Determine field name; import VivaCT scans in nested structure
         fname1 = genvarname(id.name.vct{jj});
        [stack,vox,dinfo] = importdicom([p.vct,id.vct{ii,jj}],'gauss',1.5);
         stack = permute(stack,[2,3,1]);
         VivaCT.stack.(fname1) = stack; 
         VivaCT.vox.  (fname1) = vox  ; 
         VivaCT.info. (fname1) = dinfo;
         
%          % Create reference frames for registration
%          VivaCT.r2d.(name1) = imref2d(size(VivaCT.proj. (name1)),vox,vox);
%          VivaCT.r3d.(name1) = imref3d(size(VivaCT.stack.(name1)),vox,vox,vox);
     end
     
     % Manually correct angulation and specify sample holder floor
     for jj = 1:size(id.vct,2)
        %%
        fname1     = genvarname(id.name.vct{jj});
        stack     = VivaCT.stack.(fname1);
        vox       = VivaCT.vox.  (fname1);
        RSS_imshow(stack); g = ginput(4); close(gcf);
        ang1      = atan2d(g(2,2)-g(1,2),g(2,1)-g(1,1));
        ang2      = atan2d(g(4,2)-g(3,2),g(4,1)-g(3,1));
        TF_adj    = affineTF_Maker([0,0,0],[1,1,1],[0,ang2,ang1]);
        stack_adj = imwarp(stack,TF_adj,'linear');
        RSS_imshow(stack_adj); g2 = ginput(1); close(gcf);
        VivaCT.floor.(fname1) = round(g2(2));
        VivaCT.adj.  (fname1) = {TF_adj,ang1,ang2};
        VivaCT.stack.(fname1) = stack_adj;
        VivaCT.proj. (fname1) = fPROJ(stack_adj(:,:,1:VivaCT.floor.(fname1)));
        VivaCT.r2d.  (fname1) = imref2d(size(VivaCT.proj. (fname1)),vox,vox);
        VivaCT.r3d.  (fname1) = imref3d(size(VivaCT.stack.(fname1)),vox,vox,vox);
     end
     
%      save([p.output,'Endpoints\',num2str(ii),'.mat'],'VivaCT','NIR',...
%           'id','vox','pix','-v7.3');
%  end
%  %%
%  for ii = eps
%      %% Correct projection distortion, register, and crop
%      t = tic; clear VivaCT NIR vox pix
%      load([p.output,'Endpoints\',num2str(ii),'.mat'],'VivaCT','NIR',...
%           'id');
      
     nir_proj = struct; height = struct;
     %%
     for jj = 1:size(id.vct,2)
         %%
         fname1 = genvarname(id.name.vct{jj});
         
         % 1. Obtain initial rough coalignment
         
         % 1a. Preprocess images to aid registration; create ref frames
         stack  = VivaCT.stack.(fname1)(:,:,1:VivaCT.floor.(fname1));
         if isstruct(VivaCT.vox), vox = VivaCT.vox.(fname1); else, vox = VivaCT.vox; end
         r3d_in = imref3d(size(stack),vox,vox,vox);
         proj   = stack; proj(proj < 0) = 0; proj = fPROJ(proj);
         mov    = sqrt(proj);
         r_m    = VivaCT.r2d. (fname1);
         fix    = fNORM(sqrt(NIR.im.(fname1).Post.x800));
         pix    = NIR.pix;
         r_f    = imref2d(size(fix),pix,pix);
         r      = pix/vox;
        
         % 1b. Prealign image centers via center of mass
        [tx,ty] = align_centers(mov,r_m,fix,r_f);
         iTF  = affine2d();  iTF.T(3,1:2) = [tx,ty];
         iTFi = iTF; iTFi.T = round(inv(iTFi.T),4);

%          % 2. Determine heights of each bone sample and use that to create
%          %    a projective correction of that sample
%          if   vox==.0624, se = ones(7,7,7); thresh = 2800;
%          else             se = ones(9,9,9); thresh = 3500;
%          end
%          bone  = imfill(imclose(bwareaopen(stack > thresh,...
%                         500),se),'holes');
%          cropLR = crop.(fname1);
%          % Note: Brute force fix for inconsistent orientation between
%          %       first endpoint and subsequent endpoints for FemCond
%          cropLR = round([cropLR(:,2),cropLR(:,2)+cropLR(:,4),...
%                          cropLR(:,1),cropLR(:,1)+cropLR(:,3)]...
%                        - repmat(c([1,1,3,3]),[size(cropLR,1),1]));
%          if ii~=1 && jj==1, cropLR(:,3:4) = cropLR(:,3:4)+60; end
%          fC = @(x,c,kk) x(c(kk,1):c(kk,2),c(kk,3):c(kk,4),:);
%          I  = cellfun(@str2num,id.ep(:,2))==ii;
%          epnames = id.ep(I,:);
%          epid    = cellfun(@str2num,id.ep(I,3:4));
%          nir     = NIR.im.(fname1);
%          fname2  = fieldnames(NIR.im.(fname1));
%          fname3  = fieldnames(NIR.im.(fname1).(fname2{1}));
%          
%          % 2a. Loop of each sample to determine 1) height, and 2) xy center
         ctrV = zeros(max(epid(:)),3); ctrN = ctrV(:,1:2);
         for kk = 1:max(epid(:))
             %%
             % Isolate individual sample and determine its height
             cLR   = cropLR(kk,:) + round([tx,tx,ty,ty]);             
             cmask = false(size(fix));
             cmask(cLR(1):cLR(2),cLR(3):cLR(4),:) = true;
             cmask2 = imwarp(cmask,r_f,iTFi,'Nearest','OutputView',r_m);
             bone2  = purify(imreconstruct(repmat(cmask2,...
                     [1,1,size(stack,3)]),bone));
            [x,y,z] = ind2sub(size(bone2),find(bone2));
             height.(fname1)(kk,1) = (max(z)-min(z)+1) .* vox;
             hHL    = [floor(height.(fname1)(kk)),ceil(height.(fname1)(kk))];
             hprop  =  height.(fname1)(kk)/hHL(end);
             
             % Generate a projective transform based on sample height -
             % this transform is individual to each sample and replaces the
             % original sample on the transformed image piece-wise.
             %
             % Note: Currently also adds TF0, which projects the sample all 
             %       the way back to the LiCOR imaging bed floor rather  
             %       than the sample holder floor - this theoretically 
             %       deals with all projection-related issues without an
             %       additional scaling factor being necessary.
             pTF1  = projectionTF.TF(hHL(1)); pTF2 = projectionTF.TF(hHL(2));
             pTF   = pTF1; pTF.T = pTF1.T .* (1-hprop) + pTF2.T .* hprop;
             pTF.T = pTF.T * projectionTF.iTF0.T;
             pTF.T(3,1:2) = pTF.T(3,1:2) .* pix;
             clear pTF1 pTF2
             
             % Write projected sample back into the image
             temp = imwarp(NIR.im.(fname1).Post.x800,...
                           r_f,pTF,'OutputView',r_f);
            [x,y] = ind2sub(size(temp),find(purify(...
                            cmask & (otsu(fNORM(temp)) | otsu(fix)))));
             x = [min(x)-30,max(x)+30]; y = [min(y)-30,max(y)+30];
             x(x < 1) = 1; y(y < 1) = 1;
             cropLR(kk,:) = [x,y];
             for mm = 1:length(fname2)
                 for nn = 1:length(fname3) % ISSUE HERE
                     nir.(fname2{mm}).(fname3{nn})(x(1):x(2),y(1):y(2))...
                         = temp(x(1):x(2),y(1):y(2));
                 end
             end
             nir_proj.(fname1) = nir;
             
             % Obtain centroids for each sample
             N = purify(imreconstruct(cmask,otsu(nir_proj.(fname1).Post.x800)));
             ctrV(kk,:) = regionprops(bone2,'Centroid').Centroid;
             ctrN(kk,:) = regionprops(N    ,'Centroid').Centroid;
         end
%          clear cmask hprop hHL bone2 temp
         % Note: the resulting variable 'nir_proj' now replaces the orignal
         %       NIR image for future analysis steps, as the piece-wise
         %       mapping makes it challenging to apply the inverse
         %       transform later on in the analysis.
         
         % 3. Perform registrations

         % 3a. Register VivaCT 2D projection onto NIR image
         clear ctrV2 ctrN2
        [ctrV2(:,1),ctrV2(:,2)] = intrinsicToWorld(r_m,ctrV(:,1),ctrV(:,2));
        [ctrN2(:,1),ctrN2(:,2)] = intrinsicToWorld(r_f,ctrN(:,1),ctrN(:,2));
         iTF  = fitgeotrans(ctrV2,ctrN2,'similarity');
         fix  = fNORM(sqrt(nir_proj.(fname1).Post.x800));
         TF   = imregtform(mov,r_m,fix,r_f,'similarity',opts,met,...
                          'InitialTransformation',iTF);
         prev = imwarp(mov,r_m,TF,'Linear','OutputView',r_f);
         nir_proj.(fname1).TF  = TF ;
         nir_proj.(fname1).r2d = r_f;
         
         % 3a. Print a preview of the registration
         figure,imshowpair(prev,fix,'ColorChannels','red-cyan');
         title(fname1);
%          print([p.output,'Previews\Endpoint\Endpoint',num2str(ii),'_',...
%                 fname1,'.tif'],'-dtiff','-r200');
%          close(gcf);

        % 3a. Register VivaCT 2D projection onto NIR image
         clear ctrV2 ctrN2
        [ctrV2(:,1),ctrV2(:,2)] = intrinsicToWorld(r_m,ctrV(:,1),ctrV(:,2));
        [ctrN2(:,1),ctrN2(:,2)] = intrinsicToWorld(r_f,ctrN(:,1),ctrN(:,2));
         iTF  = fitgeotrans(ctrV2,ctrN2,'similarity');
         fix  = fNORM(sqrt(NIR.im.(fname1).Post.x800));
         TF   = imregtform(mov,r_m,fix,r_f,'similarity',opts,met,...
                          'InitialTransformation',iTF);
         prev = imwarp(mov,r_m,TF,'Linear','OutputView',r_f);
         nir_proj.(fname1).TF  = TF ;
         nir_proj.(fname1).r2d = r_f;
         
         % 3a. Print a preview of the registration
         figure,imshowpair(prev,fix,'ColorChannels','red-cyan');
         title([fname1,' - Uncorrected']);
     end
     %%
     for jj = 1:1
         %%
         % 3b. Create reference frames for 3D VivaCT registration at both 
         %     VivaCT and NIR resolutions
         r = pix/vox;
         r3d_in   = imref3d( size(stack),vox,vox,vox);
         r3d_low  = imref3d([size(prev), round(size(stack,3)/r)],pix,pix,pix);
         r3d_high = imref3d([round(size(prev).*r),size(stack,3)],vox,vox,vox);
         r2d_low  = imref2d(r3d_low. ImageSize(1:2),pix,pix);
         r2d_high = imref2d(r3d_high.ImageSize(1:2),vox,vox);
         
         % 3c. Perform 3D VivaCT registrations
         TF3d = affine3d();
         TF3d.T(1:2,1:2) = TF.T(1:2,1:2);
         TF3d.T(end,1:2) = TF.T(end,1:2);
         reg_low  = imwarp(stack,r3d_in,TF3d,'Linear','OutputView',r3d_low );
         reg_high = imwarp(stack,r3d_in,TF3d,'Linear','OutputView',r3d_high);
         reg_proj = imwarp(proj,VivaCT.r2d.(fname1),...
                                        TF  ,'Linear','OutputView',r2d_low );
         
         % Crop and Save Images
         
         % Preprocess crop structure
         cropHR = round(cropLR .* r);
         %%
         % Crop images via pre-determined crop boxes
         for kk = 1:max(epid(:))
             %%
             % Which sample is this?
            [x,y] =  find(epid==kk);
             sampname = [epnames{x,1},id.name.ep{2,y+2},fname1(1)];
             
             % Does a partial image file already exist? If so, load it
             savefile = [p.output,sampname,'.mat'];
             if exist(savefile,'file')
                 if    isempty(whos('-file',savefile,'nir')), nir = struct;
                 else, nir = load(savefile,'nir'); nir = nir.nir;
                 end
                 if    isempty(whos('-file',savefile,'vct')), vct = struct;
                 else, vct = load(savefile,'vct'); vct = vct.vct;
                 end
             else, nir = struct; vct = struct;
             end
             
             % Perform cropping on all images
             nir.(fname1).pre. x700 = fC(nir_proj.(fname1).Pre. x700,cropLR,kk);
             nir.(fname1).pre. x800 = fC(nir_proj.(fname1).Pre. x800,cropLR,kk);
             nir.(fname1).pre. W    = fC(nir_proj.(fname1).Pre. W   ,cropLR,kk);
             nir.(fname1).post.x700 = fC(nir_proj.(fname1).Post.x700,cropLR,kk);
             nir.(fname1).post.x800 = fC(nir_proj.(fname1).Post.x800,cropLR,kk);
             nir.(fname1).post.W    = fC(nir_proj.(fname1).Post.W   ,cropLR,kk);
             nir.(fname1).TFproj    = TF ;
             nir.(fname1).pix       = pix;
             vct.(fname1).vox       = vox;
             vct.(fname1).LR        = fC(reg_low ,cropLR,kk);
             vct.(fname1).HR        = fC(reg_high,cropHR,kk);
             vct.(fname1).proj      = fC(reg_proj,cropLR,kk);
             %%
%              % Save sample
%              if   exist([p.output,sampname,'.mat'],'file'               )
%                    save([p.output,sampname,'.mat'],'nir','vct','-append');
%              else, save([p.output,sampname,'.mat'],'nir','vct','-v7.3'  );
%              end
         end
         
     end
%      save([p.output,'Endpoints\',num2str(ii),'.mat'],'nir_proj','height','-append');
     looptrack(numproc,L,t,['Endpoint ',num2str(ii),' (NIR/VCT)']);
     numproc = numproc + 1;
 end
 fprintf(['Successfully loaded ',num2str(numproc-1),' Endpoint.\n']);
 
 %% 5. Load in CE-µCT Scans
 
%  % Pull a list of unloaded CE-µCT scans
%  tbd = scanlist(status(:,3));
%  L = length(tbd); numproc = 1;
%  
%  % Loop through each and process if possible
%  for ii = 1:L
%      %%
%      % Search identifier to find the correct scan number
%      t = tic;
%      x = find(ismember(id.ceuct(:,1),tbd{ii}(1:2)));
%      y = find(ismember(id.name.ceuct,tbd{ii}(3:4)));
%      
%      % Check whether the scan file is listed on the identifier
%      if ~any([isempty(y),isempty(x)])
%          %%
%          % Check if the DICOM files are present; if so, load and save
%          sn = num2str(id.ceuct{x,y});
%          if ~isempty(dir([p.ceuct,sn,'\*.dcm']))
%             [ceuct.stack,ceuct.vox,ceuct.info] = importdicom([p.ceuct,sn]);
%              if   exist([p.output,tbd{ii},'.mat'],'file'           )
%                    save([p.output,tbd{ii},'.mat'],'ceuct','-append');
%              else, save([p.output,tbd{ii},'.mat'],'ceuct','-v7.3'  );
%              end
%              looptrack(numproc,nnz(~status(:,3)),t,[tbd{ii},' (CE-µCT)']);
%              numproc = numproc + 1;
%          end
%      end
%  end
%  fprintf(['Successfully loaded ',num2str(numproc-1),' CE-µCT scans.\n']);
% 
%  % Apply any necessary corrections
%  correctionator;