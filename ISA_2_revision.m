%% Startup

% Load in pre-requisite variables
load('revision.mat','whichiswhich','d','L');
d{1} = d{1}(:,1); d{2} = d{2}(:,1); d_save = d;

% Generate directory and name list
p    = 'Masks\'              ;
f    = @(x) x(1:end-4)       ;
d2   = dir2cell_isdir(dir(p));
L2   = length(d2)            ;
name = {cellfun(f,d{1},'UniformOutput',false),...
        cellfun(f,d{2},'UniformOutput',false)}; clear f

% Determine voxel information
vox = load('vox.mat','vox_train');
vox = vox.vox_train;

% Training masks
% p = 'Masks\';

% Average atlas iterations
p2 = 'Analysis - Revision\1 - Average Atlas\Atlas\'          ;
p3 = 'Analysis - Revision\1 - Average Atlas\Registrations\'  ;
d3 = dir([p3,'*.mat']); d3 = dir2cell(d3);

CI = 'CI';

%% Transform and average masks

for cc = 1:2
for at = 1:25
%%
    
% Generate filename
fn = [CI(cc),num2str(at,'%02.f')];

% Apply randomizer for this iteration
d = d_save                       ;
d = d{cc}(whichiswhich{cc}{at,1});
L = length(d)                    ;

% Load average image
avg = load([p2,fn,'.mat'],'avg'); avg = avg.avg; sz = size(avg.Image);
avg.Image = round(stackinterp(avg.Image,sz.*2,'interp3','Linear'))   ;

% For 1-sample atlases
if L==1
    
    % Load masks
    for ii = 1:length(d2)
        mask = load([p,d2{ii},'\',d{1}],'roivol');
        if contains(d{1},'L'),mask.roivol = flip(mask.roivol, 2); end
        avg.(d2{ii}) = mask.roivol;
    end
    clear mask
    
    % Save and continue
    save([p2,fn,'.mat'],'avg','d','-append');
    fprintf(['Atlas iteration ',fn,' complete\n']);
    continue
    
end

% Compile training masks
for ii = 1:length(d2)
    %%
    
    t = tic;
    
    % Compile true masks
    masks = cell(length(d),1);
    for jj = 1:length(d)
        
        mask = load([p,d2{ii},'\',d{jj}],'roivol');
        
        % Change res to 24 um
        mask.roifull = mask.roivol;
        mask.roivol = stackinterp(mask.roivol,...
            round(size(mask.roivol)./2),'interp3','Nearest');
        if contains(d{jj},'L'), mask.roivol = flip(mask.roivol, 2); end % L?R flip
        
        masks{jj} = single(mask.roivol);

    end
    clear mask

    % Perform affine transformation
    deforms = load([p3,fn,'.mat'],'deforms'); deforms = deforms.deforms;
    Rfix = imref3d(size(masks{1}),vox,vox,vox)                         ;
    I = find(~cellfun(@isempty,deforms(2:end,2)))                      ;
    for jj = 1:length(I)
        Rmov = imref3d(size(masks{I(jj)}),vox,vox,vox);
        masks{I(jj)} = imwarp(masks{I(jj)},Rmov,deforms{I(jj)+1,2},...
            'OutputView',Rfix,'interp','nearest');
    end

    % Perform sequential NRR registration
    for it = 1:size(deforms,2)-2
        for jj = 1:length(d)
            masks{jj} = imwarp(masks{jj},deforms{jj+1,it+2},'interp','nearest');
        end
    end
    
    % Average masks into an average VOI, return to 12 µm, and add to struct
    avg_atlas    = round(mean(cat(4,masks{:}),4));
    avg_atlas    = round(stackinterp(avg_atlas,size(avg_atlas) .* 2,...
                         'interp3','Nearest'));
    avg.(d2{ii}) = avg_atlas;
    
    looptrack(ii,L2,t,[d2{ii},' total completion time']);
%     mpr3D(avgatlas,add_avgatlas);
end

% Save atlas
save([p2,fn,'.mat'],'avg','d','-append');

fprintf(['Atlas iteration ',fn,' complete\n']);

end
end