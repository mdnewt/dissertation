%% Startup

onlynew = false; % If true omits redos, if false analyzes all samples

% Paths
pI = [pwd,'\..\1 - Images\'];
pM = [pwd,'\..\3 - Masks\' ];
p1 = '1 - STLs\'; p2 = '2 - PLYs\'    ; p3 = '3 - Params\';
p4 = '4 - UVH\' ; p5 = '5 - ACT Maps\';
% skip = load('skip.mat','skip'); skip = skip.skip;

%% Loop 1 - BCI Isolation and Meshing

% Determine which samples should be processed
dI = dir2cell(dir([pI,'*.mat']));
dM = dir2cell(dir([pM,'*.mat']));
d1 = dir2cell(dir([p1,'*.stl']));
fn = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false);
if onlynew, d = dM(contains(fn(dM,1:3),fn(dI,[1,2,4])) & ...
                  ~contains(fn(dM,1:3),fn(d1,[1,2,4])));
else,       d = dM(contains(fn(dM,1:3),fn(dI,[1,2,4])));
% d = d(contains(cellfun(fn,d   ,'UniformOutput',false),... % Remove skips
%         unique(cellfun(fn,skip,'UniformOutput',false))));
end; clear dI dM d1 fn
LL   =  length(d)  ;
side = {'L','R'}   ;
se   =  ones(5,5,5);

%%
for ii = 1:LL
    %%
%     try
        t  = tic;
        FT = contains(d{ii},'F'); % 0 = Tibia, 1 = Femur

        stack   = struct;
        stack.R = load([pI,d{ii}(1:2),'R',d{ii}(3:end)],'ceuct'); 
        stack.L = load([pI,d{ii}(1:2),'L',d{ii}(3:end)],'ceuct');
        stack.R = stack.R.ceuct.stack; stack.L = stack.L.ceuct.stack;
        masks   = load([pM,d{ii}]);
%
        % Isolate BCI
        for jj = 1:2
            %%
            mask = imopen(imclose(masks.cartroi.(side{jj}),se),se);
           [f,v] = computeBCI(stack.(side{jj}),mask,FT);
            figure,patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceColor','b');
            axis equal;
            print([p1,'Previews\',d{ii}(1:3),side{jj},'.jpg'],'-djpeg','-r150');
            close(gcf);
            stlwrite([p1,d{ii}(1:2),side{jj},d{ii}(3),'.stl'],f,v);
        end

        % Record iteration
        looptrack(ii,LL,t,d{ii}(1:3));
        
%     catch, warning([d{ii}(1:3),' no bueno']);
%     end
end
if LL > 0
    fprintf(['1. BCI Isolation step completed on ',num2str(LL),'sets of samples.\n']);
end

%% Loop 2 - Mesh Processing in Meshlab
% Note: Meshlab v2020.07 required to guarantee compatibility with MLX files

% Determine which samples should be processed
d1 = dir2cell(dir([p1,'*.stl']));
d2 = dir2cell(dir([p2,'*.ply']));
fn = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false);
if onlynew, d = d1(~contains(fn(d1,1:4),fn(d2,1:4)));
else,       d = d1;
end; clear d1 d2 fn
% d  = d(~contains(cellfun(fn,d   ,'UniformOutput',false),... % Remove skips
%           unique(cellfun(fn,skip,'UniformOutput',false))));
LL = length(d);
M  = '"C:\Program Files\VCG\MeshLab\meshlabserver"';

%%
for ii = 1:LL
    %%
    t = tic;
    
    % Filename and info
    name = strcrop(d{ii},'.',1);
    LR   = contains(name,'R'); % 0 = Left , 1 = Right
    FT   = contains(name,'F'); % 0 = Tibia, 1 = Femur
    
    % Process in Meshlab and write PLY
    I = ['"',pwd,'\',p1,name,'.stl"'];
    O = ['"',pwd,'\',p2,name,'.ply"'];
    if FT, S = ['"',pwd,'\meshclean_1_fem.mlx"'];
    else , S = ['"',pwd,'\meshclean_1_tib.mlx"'];
    end
   [~,log1]  = system([M,' -i ',I,' -o ',O,' -m vn -s ',S]);
    S = ['"',pwd,'\meshclean_2.mlx"'];
   [~,log2]  = system([M,' -i ',O,' -o ',O,' -m vn -s ',S]);
   
    % Record iteration
    looptrack(ii,LL,t,name);
end
if LL > 0
    fprintf(['2. Meshlab step completed on ',num2str(LL),' samples.\n']);
end
clear log1 log2 I O M S

%% Loop 3 - Parameterization of BCI

% Determine which samples should be processed
d2 = dir2cell(dir([p2,'*.ply']));
d3 = dir2cell(dir([p3,'*.mat']));
fn = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false);
if onlynew, d = d2(~contains(fn(d2,1:4),fn(d3,1:4)));
else,       d = d2; 
end; clear d2 d3 fn
% d  = d(~contains(cellfun(fn,d   ,'UniformOutput',false),... % Remove skips
%           unique(cellfun(fn,skip,'UniformOutput',false))));
LL   = length(d);

%%
for ii = 1:LL
    %%
%     try
    t = tic;
    
    % Filename and info
    name = strcrop(d{ii},'.',1);
    side = name(end-1)         ;
    LR   = contains(name,'R'    ); % 0 = Left , 1 = Right
    FT   = contains(name,'Femur'); % 0 = Tibia, 1 = Femur
    
    % Load in and process PLY
   [v,f,n] = read_ply([p2,d{ii}]);
    Mesh = makeMesh(v,f,n);
    Mesh = clipEars(Mesh) ; Mesh.name = name;
    v = Mesh.v; f = Mesh.f; n = Mesh.n;
    
    % Load in mask structure (for compartment masks)
    masks =  load([pM,d{ii}([1,2,4]),'.mat'])        ;
    comps = {'med','lat','tro'}                      ;
    comps =  comps(contains(comps,fieldnames(masks)));
    ML = struct; NC = length(comps);
    for jj = 1:NC, ML.(comps{jj}) = masks.(comps{jj}).(side); end
	
    % Run conformal parameterization
    Meshes = computeSCP(Mesh,ML);
    
    F = figure;
    if contains(Mesh.name,'L'), JJ = [2,1,3]; else, JJ = 1:3; end
    for jj = 1:NC
        subplot(1,NC,jj); M = Meshes(JJ(jj));
        patch('Faces',M.f,'Vertices',M.uv,'EdgeColor',[0,0,.5],...
              'FaceColor','b'); axis equal;
        xlim([0,max(M.uv(:,1))]); ylim([0,max(M.uv(:,2))]);
        title(M.name);
    end
    F.Units = 'centimeters'; F.PaperUnits = 'centimeters';
    F.Position(3:4) = [10*NC,12]; F.PaperPosition = F.Position;
    print([p3,'Previews\',name,'.jpg'],'-r150','-djpeg');
    %
    close(gcf);
    
    save([p3,d{ii}(1:end-3),'mat'],'Meshes');
    
    looptrack(ii,LL,t,name);
%     catch, warning([name,' no bueno']);
%     end
end
if LL > 0
    fprintf(['3. Parameterization step completed on ',num2str(ii),' samples.\n']);
end

%% Loop 4 - UVH Parameterized Image Stacks

% Determine which samples should be processed
d  = dir2cell(dir([p3,'*.mat']));
if onlynew, d2 = dir2cell(dir([p4,'*.mat'])); d = d(~contains(d,d2)); end
% d = d(~contains(cellfun(fn,d   ,'UniformOutput',false),... % Remove skips
%          unique(cellfun(fn,skip,'UniformOutput',false))));
clear d2
LL = length(d)        ;

%%
for ii = 1:LL
    %%
    t = tic;
    
    % Load data
    LR     = d{ii}(3);
    ceuct  = load([pI,d{ii}],'ceuct'  ); ceuct  = ceuct.ceuct  ;
    Meshes = load([p3,d{ii}],'Meshes' ); Meshes = Meshes.Meshes;
    mask   = load([pM,d{ii}([1:2,4:end])],'cartroi')           ;
    mask   = single(mask.cartroi.(LR))                         ;
    
    % Note: xyz_to_uvh includes application of threshold for AC isolation.
    %       For this study (conducted at 70 kVP on a Scanco µCT-40),
    %       200 and 2400 were empirically determined to be suitable.
    uv = struct;
   [uv.stack,uv.cart,uv.bwm] = xyz_to_uvh(ceuct.stack,mask,Meshes,...
                                          ceuct.vox,[200,2400]);
    uv.vox = ceuct.vox;
                                      
    save([p4,d{ii}],'uv','-v7.3');
    
    looptrack(ii,LL,t,d{ii});
end
if LL > 0
    fprintf(['4. UVH step completed on ',num2str(LL),' sets of samples.\n']);
end
gongsound;

%% Loop 5 - ACT Maps

% Determine which samples should be processed
% Note: for comparison and color scaling, lefts and rights will be
%       processed and saved together
fn = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false);
dL = dir2cell(dir([p4,'*L*.mat'])); dR = dir2cell(dir([p4,'*R*.mat']));
d  = fn(dL(contains(fn(dL,[1,2,4]),fn(dR,[1,2,4]))),[1,2,4:8]);
LL = length(d);
par = parula; par(1,:) = [0,0,0] ; % Colormap with black for 0 thickness
lat = {'L','R'};

%%
for ii = 1:LL
    %%
    
    t = tic; 
    
    % Filename
    name = strcrop(d{ii},'.',1);
    
    % Loop of left and right limbs
    uv = struct; submask = uv; ACT = uv; atten = uv;
    for jj = 1:2
        %%
        
        % Load files
        LR = lat{jj};
        uv.(LR)     = load([p4,d{ii}(1:2),LR,d{ii}(3:end)]);
        if isfield(uv.(LR),'submask'), submask.(LR) = uv.(LR).submask; end
        uv.(LR)     = uv.(LR).uv;
%         uv.(LR).vox = .012;

        % Apply submasks, if present
        if exist('submask','var') &&  isfield(submask,LR)
            for kk = 1:length(submask.(LR))
                if ~isempty(submask.(LR){kk})
                    uv.(LR).cart{kk} = uv.(LR).cart{kk} & ~submask.(LR){kk};
                end
            end
        end

        % Compute cartilage thickness map
        L2  = length(    uv.(LR).stack     );
        mid = round(size(uv.(LR).stack,3)/2);
        ACT.(LR) = cell(1,L2);
        for kk = 1:L2 
            ACT.(LR){kk} = sum(uv.(LR).cart{ kk},3) .* uv.(LR).vox;
    %             ACT.(lat{jj}){kk}(isnan(uv.(lat{jj}).stack{kk}(:,:,mid))) = NaN      ;
        end

        % Compute cartilage attenuation map
        atten.(lat{jj}) = cell(1,L2);
        for kk = 1:L2
            temp = uv.(lat{jj}).stack{kk}       ;
            temp( ~uv.(lat{jj}).cart{ kk}) = NaN;
            atten.(lat{jj}){kk} = mean(temp,3,'omitnan')       ;
            atten.(lat{jj}){kk}(isnan(ACT.(lat{jj}){kk})) = NaN;
        end
    end
    
    % Define scaling factor based on F/T
    % Note: Scaling factor was empirically determined for rat AC at the
    %       voxel size used in this study - it should be adjusted for use
    %       in other situations.
    if strcmp(d{ii}(3),'F'), sf = 5; else sf = 3.5; end
    
    % Generate Preview Figures
    figure;
    for jj = 1:L2-1
        subplot(L2-1,1,jj);
        imshowpair(ACT.L{jj}.*sf,ACT.R{jj}.*sf,'Scaling','none','Method','Montage');
        colormap(par); caxis('auto'); title('L <-> R');
    end
    set(gcf,'PaperUnits','inches','PaperSize',[2.5,8],'PaperPositionMode',...
        'manual','PaperPosition',[.1,.1,2.3,7.8]);
    print([p5,'Previews\',name,'.tif'],'-dtiff','-r150');
    close(gcf);
    
    save([p5,d{ii}],'ACT','atten');
    
    looptrack(ii,LL,t,d{ii});
end
if LL > 0
    fprintf(['5. ACT Map step completed on ',num2str(LL),' sets of samples.\n']);
end
gongsound;

% Note: at this point, run 'touchup_loop.mat' and re-compute ACT maps

%% Loop 6 - Compute Results

% Determine which samples should be processed
f   = @(x,y) [x(1:2),y,x(3:end)];
d   =   dir2cell(dir([p5,'*.mat'])); L1 = length(d);
d2  =  [cellfun(f,d,repmat({'L'},[L1,1]),'UniformOutput',false),...
        cellfun(f,d,repmat({'R'},[L1,1]),'UniformOutput',false)];
L2  =   numel(d2);
d2  =   reshape(d2',[L2,1]);
lat = {'L','R'};

% Prepare results cell
comps = cell(L2+1,5,5);
comps(2:end,1,:) = repmat(d2,[1,1,5]);
comps(1,1,:) = {'Mean Thickness','Sa','normSa','Cart Atten','BCI SA'};
comps(1,2:end,:) = repmat({'Med','Lat','Tro','Whole'},[1,1,5]);
names = cell(L2+1,3); names(1,:) = {'An.#','Lat','F/T'};

%%
%% 
for ii = 1:L1
    %%
    
    t = tic;
    
    % Get filename
    file = strsplit(d{ii},'.'); file = file{1};
    name = strsplit(file,'_');
    
    % Load in sample
    H = load([p5,d{ii}]);
    
    for kk = 1:2
        %%
        ACT   = H.ACT.(  lat{kk});
        atten = H.atten.(lat{kk});
        Meshes = load([p3,d{ii}(1:2),lat{kk},d{ii}(3:end)],'Meshes'); Meshes = Meshes.Meshes;

        N = length(ACT)-1;

        % Prepare mask and voxlist
        masks = cell(2,N); voxList = cell(2,N);
        for jj = 1:N
            masks{  1,jj} = ~isnan(    ACT{    jj});
            masks{  2,jj} = ~isnan(    atten{  jj});
            voxList{1,jj} =  ACT{  jj}(masks{1,jj});
            voxList{2,jj} =  atten{jj}(masks{2,jj});
        end

        % Compute compartments
        area = mesh_SA(Meshes,0.012); % Should load voxel size
        for jj = 1:N
            JJ         = find(contains(squeeze(comps(1,:,1)),Meshes(jj).name));
            cart_thick = mean(voxList{1,jj});
            Sa         = mean(abs(voxList{1,jj}-cart_thick));
            normSa     = Sa/cart_thick      ;
            cart_atten = mean(voxList{2,jj});
            comps((2*ii)+kk-1,JJ,:) = {cart_thick,Sa,normSa,cart_atten,area(jj)};
        end

        % Compute whole joint
        allVox = cell(2,1);
        for jj = 1:2 % For height and atten
           allVox{jj} = cat(1,voxList{jj,:});
        end
        cart_thick = mean(allVox{1});
        Sa         = mean(abs(allVox{1}-cart_thick));
        normSa     = Sa/cart_thick      ;
        cart_atten = mean(allVox{2});
    %     SCBthick   = mean(allVox{3});
    %     SCBdens    = mean(allVox{4});
        comps((2*ii)+kk-1,5,:) = {cart_thick,Sa,normSa,cart_atten,sum(area)};
%         names((2*ii)+kk-1,:  ) = {[name{1},lat{kk}]};
    end
    looptrack(ii,L1,t,d{ii});
end
%
% Restructure
toprow = {'','AC.Th. (mm)'   ,'','','','Sa (mm)'   ,'','','',...
             'normSa'        ,'','','','Atten (HU)','','','',...
             'BCI Area (mm2)','','',''};
% comps = cat(1,toprow,cat(2,comps(:,1,1),names,...
%               comps(:,2:5,1),comps(:,2:5,2),comps(:,2:5,3),...
%               comps(:,2:5,4),comps(:,2:5,5)));
comps = cat(1,toprow,cat(2,comps(:,1,1),...
              comps(:,2:5,1),comps(:,2:5,2),comps(:,2:5,3),...
              comps(:,2:5,4),comps(:,2:5,5)));
comps{2,1} = 'Sample';
%%
save([p2,'femur_paramdata.mat'],'comps');