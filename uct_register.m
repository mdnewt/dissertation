 %% 
 % Setup
 redo = true; % Should already-processed files be redone?
 p = '1 - Images\'              ;
 f = @(x,y) cellfun(@(x)isempty(whos('-file',[p,x],y)),x);
 d =  dir2cell(dir([p,'*.mat'])); d = d(~f(d,'vct'));
 if  ~redo, d = d(f(d,'reg')); end
 L =  length(d)                       ;
[opts,met] = imregconfig('monomodal'); 
load('ep1.mat'); % Brute force fix
 
 % Define In-Line Functions
 fCENT = @(x) regionprops(single(x>3000),'Centroid').Centroid;
 fPROJ = @(x) sqrt(sum(x(:,:,1:round(size(x,3)/2)).^2,3));
 
 % Define rotation required to orient CE-µCT and VivaCT scans consistently
 imrot = {'FemCond','FemTroch','Tibia';90,-90,-90};
  
 %%
 for ii = 9:12
     %%
     t = tic;
     
     % Load image file, perform initial queries, create reg structure
     im   = load([p,d{ii}]); ceuct = im.ceuct; vct = im.vct; clear im
     voxF = ceuct.vox;
     name = fieldnames(vct);
     reg = struct;
     
     %%
     % Loop of all image subsets
     for jj = 1:length(name)
         %%
         % Create variables necessary for registration
         voxM = vct.(name{jj}).vox;
         r    = voxF/voxM         ;
         mov  = vct.(name{jj}).HR ;
         if ~any(contains(ep1,d{ii}(1:2))) && strcmp(name{jj},'FemCond')
             rot = -imrot{2,cellfun(@(x)strcmp(x,name{jj}),imrot(1,:))};
         else
             rot =  imrot{2,cellfun(@(x)strcmp(x,name{jj}),imrot(1,:))};
         end
         fROT = @(x) imrotate(x,rot);
         fix  = fROT(stackinterp(ceuct.stack,round(size(ceuct.stack).*r))); % Adjust ceuct to vct res
         rM   = imref3d(size(mov),voxM,voxM,voxM);
         rF   = imref3d(size(fix),voxM,voxM,voxM);
         
         % Compute initial transforma by aligning centers of bone threshold
         T    = (fCENT(fix) - fCENT(mov)) .*voxM ;
         iTF  =  affine3d(); iTF.T(end,1:3) = T  ;
%         [tx,ty,tz] = align_centers(mov,rM,fix,rF,'geo');
%          iTF  =  affine3d(); iTF.T(end,1:3) = [tx,ty,tz];
         
%          test = imwarp(mov,rM,iTF ,'Linear','OutputView',rF);
%          mpr3D(test,fix > 2800);
        
        %%
         % Compute and apply forward and inverse registration
         % Note: Study of the registrations done using a rigid transform
         %       demonstrated that the CE-uCT scan was slightly smaller 
         %       than the vct scan, though the reason these do not align 
         %       precisely is not known. This has been fixed by using a
         %       similarity transform instead of rigid.
         TF  = imregtform(mov,rM,fix,rF,'similarity',opts,met,...
              'InitialTransformation',iTF);
%          TF  = imregtform(mov,rM,fix,rF,'rigid',opts,met,...
%               'InitialTransformation',iTF);
         TFi = TF; TFi.T = round(inv(TFi.T),4);


        test2 = imwarp(mov,rM,TF ,'Linear','OutputView',rF);
        mpr3D(test2,fix > 2800);
%%
         reg.(name{jj}).forward = imwarp(mov,rM,TF ,'Linear','OutputView',rF);
         reg.(name{jj}).inverse = imwarp(fix,rF,TFi,'Linear','OutputView',rM);
         reg.(name{jj}).TF. forward = TF ;
         reg.(name{jj}).TF. inverse = TFi;
         reg.(name{jj}).r2d.forward = rF ;
         reg.(name{jj}).r2d.inverse = rM ;
         temp = reg.(name{jj}).inverse; temp(temp==0) = -1000; sl = 96;
         figure,imshowpair(temp(:,:,sl),mov(:,:,sl),'Montage');
         figure,imshowpair(temp(:,:,sl),mov(:,:,sl),'ColorChannel','red-cyan');
%          figure,imshowpair(fPROJ(reg.(name{jj}).forward),fPROJ(fix),'ColorChannel','red-cyan');
%          figure,imshowpair(fPROJ(fNORM(reg.(name{jj}).inverse)),fPROJ(fNORM(mov)),'ColorChannel','red-cyan');
%          figure,imshowpair(sum(reg.(name{jj}).inverse,3),sum(mov,3),'Montage');
         title([d{ii}(1:end-4),' - ',name{jj}]);
         print([p,'Previews\Registration\reg_',d{ii}(1:end-4),'_',...
             name{jj},'.tif'],'-dtiff','-r150');
         close(gcf);
     end
     save([p,d{ii}],'reg','-append');
     looptrack(ii,L,t,d{ii}(1:end-4));
 end