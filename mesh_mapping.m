%% Define paths, directories, in-line functions, etc.
 
% Paths and directories
pI  = [pwd,'\..\1 - Images\'                   ]       ;
pM  = [pwd,'\..\3 - Masks\'                    ]       ;
pP  = [pwd,'\..\4 - ACT Analysis\3 - Params\'  ]       ;
pU  = [pwd,'\..\4 - ACT Analysis\3 - UVH\'     ]       ;
pA  = [pwd,'\..\4 - ACT Analysis\5 - ACT Maps\']       ;
p1  = '1 - Registration\'                              ;
p2  = '2 - Mesh Mapping\'                              ;
fn  = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false);
dI  =  dir2cell(dir([pI,'*.mat']))                     ;
d5  =  dir2cell(dir([pA,'*.mat']))                     ;
d   =  dI(contains(fn(dI,[1,2,4]),fn(d5,1:3)))         ;
Lii =  length(d); clear dI d5 fn

% Define In-Line Functions
fCENT = @(x  )  regionprops(single(x>3000),'Centroid').Centroid;
fPROJ = @(x  )  sqrt(sum(x(:,:,1:round(size(x,3)/2)).^2,3));
fCTR  = @(x,v)  repmat(x([2,1,3])./2,[length(v),1]);
fNORM = @(x  ) (x-min(x(:)))./(max(x(:))-min(x(:)));
 
opt = cell(1,2); met = opt;
[opt{1},met{1}] = imregconfig('multimodal');
[opt{2},met{2}] = imregconfig('monomodal' );

projectionTF = load([pwd,'\..\..\..\..\Aim 2\Projection Phantom\projectionTF.mat']);

id = load('S:\OrthopaedicsProjects\Data\OrthopaedicsLab\Projects\Newton\Disseration Research\Aim 3\Imaging Data\The MATLABs\1 - Images\Endpoints\5.mat','id'); id.id;
id = id.id;
% id  = load([pI,'Endpoints\5.mat'],'id'); id = id.id;
ep1 = id.ep(contains(id.ep(:,2),'1'),1);

% Other definitions
imrot  = {'FemCond','FemTroch','Tibia';90,-90,-90};
pairs  = {'Med',{'FemCond','Tibia'};'Lat',{'FemCond','Tibia'};'Tro',{'FemTroch'}};
nir_tp = {'pre','post'};
par    =  parula; par(1,:) = [0,0,0];
% Note: 'imrot' defines the z-axis rotation necessary to align the CE-µCT
%        scan to the VivaCT.
%       'pairs' specifies which image orientation (FemCond/FemTroch/Tibia)
%        is used for each compartment of the mesh.

%% Outer Loop of All Specimens
%17 issue
for ii = 1:Lii
    %% Load and parse inputs
    t = tic;
    im     = load([pI,d{ii}]); ceuct = im.ceuct; reg = im.reg;
                               vct   = im.vct  ; nir = im.nir;
    masks  = load([pM,d{ii}([1:2,4:end])]);                       
    ACT    = load([pA,d{ii}([1,2,4:end])],'ACT'); ACT = ACT.ACT;
    Meshes = load([pP,d{ii}]); Meshes = Meshes.Meshes;
    step1  = load([p1,d{ii}]); TF3 = step1.TF3; cect_proj = step1.cect_proj;
                               pre = step1.pre; post      = step1.post     ;
    Ljj = length(Meshes)-1; LR = d{ii}(3);
    fnames = fieldnames(vct) ;
    out = cell(1,length(fnames));
    Mrot = Meshes(1:Ljj); ang = cell(Ljj,1); Ivis = ang; orient = ang; nir_uv = ang;
    clear im step1
    
    % Transfer mesh from CE-µCT Image onto NIR Image
    %
    F = struct; fA = @(x) findall(x,'type','axes'); clear Mvis b
    for jj = 1:Ljj
        %%
        % Define variables from original CEµCT mesh
        v = Meshes(jj).v; f = Meshes(jj).f; n = Meshes(jj).n;
        Mrot(jj).u = Mrot(jj).uv(:,1:2);
        
        % Determine which image orientation will be used
        fname      = pairs{contains(pairs(:,1),Meshes(jj).name),2:end};
        fname      = fname{ismember(fname,fnames)};
        orient{jj} = fname;
        TFcorr     = TF3.(fname);
        x800       = nir.(fname).post.x800;
        if ~isfield(F,fname)
            F.(fname)      = figure; imshow(x800,[]); hold on;
            F.(fname).Name = fname;
        end

        % Define transformation matrices, reference frames, and sizes
        if ~any(contains(ep1,d{ii}(1:2))) && strcmp(fname,'FemCond')
            rot = -imrot{2,cellfun(@(x)strcmp(x,fname),imrot(1,:))};
        else
            rot =  imrot{2,cellfun(@(x)strcmp(x,fname),imrot(1,:))};
        end
        fROT = @(x) imrotate(x,rot);
        vox1 = ceuct.vox                                       ;
        vox2 = vct.(fname).vox                                 ;
        vox3 = nir.(fname).pix                                 ;
        TF1  = affineTF_Maker([0,0,0],[1,1,1],[rot,0,0])       ;
        TF2  = reg.(fname).TF. inverse; TF2.T = TF2.T*TFcorr.T ;
        R1   = imref3d(size(ceuct.stack),vox1,vox1,vox1)       ;
        R2   = reg.(fname).r2d.forward                         ;
        R3   = reg.(fname).r2d.inverse                         ;
        sz1  = R1.ImageSize .* vox1; ctr1 = fCTR(sz1,v)        ;
        sz2  = R2.ImageSize .* vox2; ctr2 = fCTR(sz2,v)        ;
%         sz3  = R3.ImageSize .* vox2; ctr3 = fCTR(sz3,v)        ;

        % Step 1: Rotate vertices from CE-µCT to VivaCT orientation
       [x,y,z] = intrinsicToWorld(R1,v(:,1),v(:,2),v(:,3)); v2 = [x,y,z];
        v2     = transformPointsForward(TF1,v2-ctr1)+ctr2;
%        [x,y,z] = worldToIntrinsic(R2,v2(:,1),v2(:,2),v2(:,3));
%         RSS_imshow(reg.(fname).forward,3); hold on; scatter(x,y);
%
        % Step 2: Perform inverse registration on the vertices
        v3     = transformPointsForward(TF2,v2);
       [x,y,z] = worldToIntrinsic(R3,v3(:,1),v3(:,2),v3(:,3));
%         RSS_imshow(reg.(fname).inverse,3); hold on; scatter(x(b),y(b));
%         RSS_imshow(reg.(fname).inverse,1); hold on; scatter(z(b),x(b));
        %
        % Step 3: Change the resolution to match the NIR image
        v4 = [x,y,z] .* (vox2/vox3);
%         figure,imshow(nir.(fname).post.x800,[]); hold on; scatter(v4(b,1),v4(b,2));
%         figure,imshow(vct.(orient{jj}).proj,[]); hold on; scatter(v4(:,1),v4(:,2));

        % Write fully translated vertices to 'v' reference frame
        Mrot(jj).v = v4;

        % Rotate normal vectors and calculate angle in sagittal plane
        n = n * (TF1.T(1:3,1:3) * TF2.T(1:3,1:3));
        if ii==90 && jj==1, n = -n; end
        
        ang{jj} = atan2d(n(:,3),n(:,1));
        Mrot(jj).n = n;

        % Eliminate any vertices with positive normal angles (these
        % correspond to surfaces 'underneath' the imaging surface)
        Iv = ang{jj} < -18; If = find(any(ismember(f,find(Iv)),2));
       [Mvis(jj),Ivis{jj}] = subMesh(Mrot(jj),If);
        b{jj} = compute_boundary(Mvis(jj).f);
        plot(fA(F.(fname)),Mvis(jj).v(b{jj},1),Mvis(jj).v(b{jj},2),...
            'LineWidth',0.75);
        
%         figure,quiver(Mvis(jj).v(:,3),Mvis(jj).v(:,1),Mvis(jj).n(:,3),Mvis(jj).n(:,1)); axis equal
%         figure,patch('Vertices',Mvis(jj).v,'Faces',Mvis(jj).f); axis equal

    end
    
    % Save progress and print previews
    save([p2,d{ii}],'Mvis','Mrot','b','orient','Ivis','ang');
    Ffnames = fieldnames(F);
    for jj = 1:length(Ffnames)
        print(F.(Ffnames{jj}),[p2,'Previews\Boundaries\',d{ii}(1:end-4),...
             '_',Ffnames{jj},'.tif'],'-dtiff','-r200');
        close(F.(Ffnames{jj}));
    end
    % Write boundaries to Mvis
    for jj = 1:Ljj, Mvis(jj).b = b{jj}; end
    looptrack(ii,Lii,t);
end
%%
for ii = 1:Lii
    %%
    % Loop 3: Map NIR Intensity Values onto Parameterized Image
    t=tic;
    ACT = load([pA,d{ii}([1,2,4:end])],'ACT'); ACT = ACT.ACT;
    nir = load([pI,d{ii}],'nir'); nir = nir.nir;
    LR  = d{ii}(3);
    M = load(  [p2,d{ii}],'Mvis','Mrot','b','Ivis','orient');
    B = M.b; Mrot = M.Mrot; Mvis = M.Mvis; Ivis = M.Ivis; orient = M.orient; clear M
    Ljj = length(orient);
    nir_interp = cell(1,Ljj);
    Fsp = figure;
    F = struct; fA = @(x) findall(x,'type','axes');
    for jj = 1:Ljj
        %%
        
        % Display stuff
        fname = orient{jj};
        x800  = nir.(fname).post.x800;
        if ~isfield(F,fname)
            F.(fname)      = figure;
            subplot(1,2,1); imshow(x800,[]); hold on;
            subplot(1,2,2); imshow(x800,[]); hold on;
            F.(fname).Name = fname;
        end
        
        % Note: NIR projections measure at the cartilage surface, not the
        %       BCI - therefore, the vertices, which are located on the
        %       BCI, need to be moved outwards along the surface normals as
        %       a function of AC.Th.
        uvACT = Mrot(jj).u; % Unchanged from original mesh
        mx    = max(uvACT,[],1) ;
        matSize = size(ACT.(LR){jj},[2,1]);
        uvACT2       = zeros(size(uvACT))             ;
        uvACT2(:,1)  = uvACT(:,1).*(matSize(1)/mx(1)) ;
        uvACT2(:,2)  = uvACT(:,2).*(matSize(2)/mx(2)) ;
        uvACT2(:,1:2)= uvACT2(:,1:2) + 2              ;
        ACTh{jj} = interp2(ACT.(LR){jj}./vox3,uvACT2(:,1),uvACT2(:,2));
        ACTh{jj}(isnan(ACTh{jj})) = 0;
%         figure,scatter(uvACT(:,1),uvACT(:,2),[],ACTh);
%         clear uvACT uvACT2 matSize mx

        % Define mesh parameters 
        v = Mvis(jj).v; f = Mvis(jj).f; uv = Mvis(jj).u;
        n = Mvis(jj).n; b = B{   jj};
        v = v + (n .* repmat(ACTh{jj}(Ivis{jj}(:,2)),[1,3]));
        uvROT = Mrot(jj).uv(:,1:2);
        Mvis(jj).v2 = v;
        
        % Resize uv
        uv2       = zeros(size(uv))            ;
        uv2(:,1)  = uv(:,1).*(matSize(1)/mx(1));
        uv2(:,2)  = uv(:,2).*(matSize(2)/mx(2));
%         uv2(:,1:2)= uv2(:,1:2)                 ;
        uv2(:,1:2)= uv2(:,1:2) + 2             ;

        % Discretize
%        [x,y] = meshgrid(1:matSize(1)+3,1:matSize(2)+3);
       [x,y] = meshgrid(1:matSize(1),1:matSize(2));
        %
%         nir_interp = cell(1,2);
        for kk = 1:2 % Broken, only writes nir_interp for post
            %%
            % Define NIR image for this iteration
            nir_im = nir.(orient{jj}).(nir_tp{kk}).x800;
            
            % Interpolate NIR values at each vertex and discretize
            nir_interp{jj} = double(interp2(nir_im,v(:,1),v(:,2),'linear'));
            sInt = scatteredInterpolant(uv2(:,1),uv2(:,2),nir_interp{jj},...
                                       'linear','none')                ;
            nir_uv{jj} = sInt(x,y)                                     ;
            
            % Remove values outside of the boundary
            sz                  = size(ACT.(LR){jj})                      ;
            bwMask              = poly2mask(uv2(b,1),uv2(b,2),sz(1),sz(2));
            perim               = bwperim(bwMask)                         ;
            nir_uv{jj}(~bwMask) = NaN                                     ;
            
            % Write temporary variables into data structures
            nir.(orient{jj}).(nir_tp{kk}).uv = nir_uv{jj}    ;
            Mvis(jj).nir.    (nir_tp{kk})    = nir_interp{jj};
        end
    end
    for jj = 1:Ljj
        %%
        % Figure Stuff
        fname = orient{jj};
        v = Mvis(jj).v2;
        fMM = @(x,mn,mx) (x-mn)/(mx-mn);
        figure(Fsp),subplot(1,Ljj,jj),imshow([fNORM(ACT.(LR){jj}),...
                    fNORM(nir_uv{jj})]); colormap(parula);
        clear RGB
       [RGB(:,1),RGB(:,2),RGB(:,3)] = ind2rgb(uint8(...
           fMM(ACTh{jj}(Ivis{jj}(:,2))*vox3,0,.2)*255),par); % 
        A = fA(F.(fname));
        patch(A(1),'Vertices',v(:,1:2),'Faces',Mvis(jj).f,'EdgeColor',...
            'none','FaceVertexCData',RGB,'FaceColor','interp');
        clear RGB
        if jj==3, mn = min(nir_interp{jj});
                  mx = max(nir_interp{jj});
        else,     mn = min([nir_interp{1};nir_interp{2}]);
                  mx = max([nir_interp{1};nir_interp{2}]);
        end
       [RGB(:,1),RGB(:,2),RGB(:,3)] = ind2rgb(uint8(...
           fMM(nir_interp{jj},mn,mx)*255),plasma);
        patch(A(2),'Vertices',v(:,1:2),'Faces',Mvis(jj).f,'EdgeColor',...
            'none','FaceVertexCData',RGB,'FaceColor','interp');
        clear RGB

    end
    save([p2,d{ii}],'nir_uv','Mvis','ACTh','nir_interp','-append')    ;
    Ffnames = fieldnames(F);
    for jj = 1:length(Ffnames)
        print(F.(Ffnames{jj}),[p2,'Previews\Overlays\',d{ii}(1:end-4),...
             '_',Ffnames{jj},'.tif'],'-dtiff','-r200');
        close(F.(Ffnames{jj}));
    end
    print(Fsp,[p2,'Previews\Maps\',d{ii}(1:end-4),'.tif'],'-dtiff','-r200');
    close(Fsp);
    
    looptrack(ii,Lii,t,d{ii}(1:end-4));
end