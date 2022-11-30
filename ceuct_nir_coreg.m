 %% Define paths, directories, in-line functions, etc.
  
 % Paths and directories
 pI = [pwd,'\..\1 - Images\'                   ]       ;
 pM = [pwd,'\..\3 - Masks\'                    ]       ;
 pP = [pwd,'\..\4 - ACT Analysis\3 - Params\'  ]       ;
 pA = [pwd,'\..\4 - ACT Analysis\5 - ACT Maps\']       ;
 p1 = '1 - Registration\'                              ;
 p2 = '2 - Mesh Mapping\'                              ;
 fn = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false);
 dI =  dir2cell(dir([pI,'*.mat']))                     ;
 d5 =  dir2cell(dir([pA,'*.mat']))                     ;
 d  =  dI(contains(fn(dI,[1,2,4]),fn(d5,1:3)))         ;
 L  =  length(d); clear dI d5 fn
 
 % Define In-Line Functions
 fCENT = @(x  )  regionprops(single(x>3000),'Centroid').Centroid;
 fPROJ = @(x  )  sqrt(sum(x.^2,3));
 fCTR  = @(x,v)  repmat(x([2,1,3])./2,[length(v),1]);
 fNORM = @(x  ) (x-min(x(:)))./(max(x(:))-min(x(:)));
 
 opt = cell(1,2); met = opt;
[opt{1},met{1}] = imregconfig('multimodal');
[opt{2},met{2}] = imregconfig('monomodal' );

 projectionTF = load([pwd,'\..\..\..\..\Aim 2\Projection Phantom\projectionTF.mat']);

 % Other definitions
 imrot  = {'FemCond','FemTroch','Tibia';90,-90,-90};
 pairs  = {'Med',{'FemCond','Tibia'};'Lat',{'FemCond','Tibia'};'Tro',{'FemTroch'}};
 nir_tp = {'pre','post'};
 par    =  parula; par(1,:) = [0,0,0];
 % Note: 'imrot' defines the z-axis rotation necessary to align the CE-µCT
 %        scan to the VivaCT.
 %       'pairs' specifies which image orientation (FemCond/FemTroch/Tibia)
 %        is used for each compartment of the mesh.

 %% Loop of All Specimens
 for ii = 1:L
    %% Load and parse inputs
    t  = tic;
    im = load([pI,d{ii}]); ceuct = im.ceuct; reg = im.reg;
                           vct   = im.vct  ; nir = im.nir;
    fnames = fieldnames(vct)       ; clear im
    
    % Create output structures
    TF3    = struct; cect_proj = TF3; pre = TF3; post = TF3;
    
    % For each scan orientation, perform registrations of: 
    % 1) CE-µCT image onto the post-wash NIR image (fine-tuned from vct)
    % 2) Pre-wash NIR image onto the post-wash NIR image
    for jj = 1:length(fnames)
        %%
        
        fname = fnames{jj};
        
        % 1. CE-µCT onto post-wash NIR image
        mov  = reg.(fname).inverse; mov(mov<0) = 0;
        mov  = fNORM(sqrt(sqrt(fPROJ(mov))));
        vox  = reg.(fname).r2d.inverse.PixelExtentInWorldX;
        rM   = imref2d(size(mov),vox,vox);
        x800 = nir.(fname).post.x800; post.(fname) = x800;
        fix  = fNORM(sqrt(sqrt(x800)));
        pix  = nir.(fname).pix;
        rF   = imref2d(size(fix),pix,pix);
%         figure,imshowpair(mov,rM,fix,rF);
%         W    = imopen(nir.(fname).post.W,ones(5,5,5));
%         W2   = imregister(W,x800,'rigid',opt{1},met{1});
%         fix  = fNORM(W2); fix(~otsu(x800)) = 0;
        %
        TF = imregtform(single(otsu(mov)),rM,single(otsu(fix)),rF,...
                       'rigid',opt{2},met{2});
        cect_proj.(fname) = imwarp(reg.(fname).inverse,rM,TF,'Linear','OutputView',rF);
        check             = imwarp(mov                ,rM,TF,'Linear','OutputView',rF);
        figure,imshowpair(check,fix,'ColorChannels','red-cyan');
        mov2  = reg.(fname).inverse; mov2(mov2<0) = 0;
        mov2  = fNORM(fPROJ(mov2));
        check             = imwarp(mov2                ,rM,TF,'Linear','OutputView',rF);
        figure,imshowpair(check,x800,'Montage');
%         figure,subplot(1,2,1); imshowpair(check ,fix,'ColorChannels','red-cyan');
        title('CEuCT on NIR x800');
        
        % Convert 2d affine transformation to 3d
%         TF3.(fname) = affine3d([TF.T(1:2,[1:3,3]);0,0,1,0;TF.T(3,1:2),0,1]);
        
%         % 2. Pre-wash onto post-wash NIR image
%         mov         = imgaussfilt(fNORM(sqrt(sqrt(nir.(fname).pre .x800))),3);
%         fix         = imgaussfilt(fNORM(sqrt(sqrt(x800                 ))),3);
%         rM          = imref2d(size(mov),pix,pix);
%         TF          = imregtform(mov,rM,fix,rF,'rigid',opt{2},met{2});
%         pre.(fname) = imwarp(nir.(fname).pre.x800,rM,TF,'Linear','OutputView',rF);
%         check       = imwarp(mov                 ,rM,TF,'Linear','OutputView',rF);
% %         subplot(1,2,2),imshowpair(check,fix,'ColorChannels','red-cyan');
% %         title('Pre on Post NIR x800 (w/ Gauss s=3)');
% %         print([p1,'Previews\',d{ii}(1:end-4),'_',fname,'.tif'],'-dtiff','-r200');
%         close(gcf);
%         clear mov fix pix TF rM tF
    end
%     save([p1,d{ii}],'cect_proj','pre','post','TF3');
    looptrack(ii,L,t,d{ii}(1:end-4));
 end