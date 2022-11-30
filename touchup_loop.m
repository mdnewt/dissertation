%% Startup

onlynew = true;

% Paths
pI = [pwd,'\..\1 - Images\'];
pM = [pwd,'\..\3 - Masks\' ];
p1 = '1 - STLs\'; p2 = '2 - PLYs\'    ; p3 = '3 - Params\';
p4 = '4 - UVH\' ; p5 = '5 - ACT Maps\';

% Directories
f  = @(x  ) cellfun(@(x) isempty(whos('-file',[p4,x],'submask')),x);
fn = @(x,y) cellfun(@(x) x(y),x,'UniformOutput',false             );
d4 = dir2cell(dir([p4,'*.mat']));
d  = dir2cell(dir([p5,'*.mat']));
if onlynew, d4 = d4(f(d4)); d = d(contains(fn(d,1:3),fn(d4,[1,2,4]))); end 
LL  = length(d); %clear d4 f fn

% Colormap and anonymous fxns
par = parula; par(1,:) = [0,0,0]            ;
aim = @(x,y) any(ismember(x,y))             ;
tbx = @(ax ) set(ax.Toolbar,'Visible','off');

ldl_string = {'1. L med','2. R med','3. L lat',...
              '4. R lat','5. L tro','6. R tro'};

%% Loop
for ii = 1:LL
%%   
    t = tic;
    
    % Load thickness maps
    ACT = load([p5,d{ii}],'ACT'); ACT = ACT.ACT; LC = length(ACT.L)-1;
    
    % Define scaling factor based on F/T
    % Note: Scaling factor was empirically determined for rat AC at the
    %       voxel size used in this study - it should be adjusted for use
    %       in other situations.
    if strcmp(d{ii}(3),'F'), sf = 5; else sf = 3.5; end
    
    % Display ACT map
    figure;
    for jj = 1:LC
        subplot(LC,1,jj); n = [jj*2-1,jj*2];
        imshowpair(ACT.L{jj}.*sf,ACT.R{jj}.*sf,'Scaling','none','Method','Montage'); % PROBLEM HERE
        tbx(gca); colormap(par); caxis('auto');
        title([num2str(n(1)),' <-> ',num2str(n(2))]);
    end
    qdl = questdlg('Is touch up required?','Touchup?','Yes','No','Rescan','No');
    
    % Check for 'yes' or 'no' responses and proceed accordingly
    switch qdl
    case 'No'    , submask = []; close(gcf);
    case 'Rescan', submask = []; close(gcf);
    case 'Yes'
        %%
        % Load UV stack for this sample
        uv = struct;
        uv.R = load([p4,d{ii}(1:2),'R',d{ii}(3:end)],'uv');
        uv.L = load([p4,d{ii}(1:2),'L',d{ii}(3:end)],'uv');
        uv.R = uv.R.uv; vox = uv.R.vox; uv.L = uv.L.uv; vox = uv.L.vox;
        
        % Outer while loop allows entire sample to be retried or abandoned
        submask =  struct;
        qdl     = 'Retry';
        while strcmp(qdl,'Retry')
            %%
            % Select which regions require touch-ups
            
            ldl = listdlg('ListString',ldl_string(1:LC*2),...
                          'PromptString','Which ones?');
            close(gcf);

            % Loop through each region and perform corrections
            for jj = 1:length(ldl)
                %%

                % Determine which region is up next and display it
                if aim([1,3,5],ldl(jj))
                      I  = (ldl(jj)-1)/2+1; act_temp = ACT.L{I}; LR = 'L';
                      im = uv.L.stack{I}  ; ac = uv.L.cart{  I};
                else, I =  ldl(jj)/2      ; act_temp = ACT.R{I}; LR = 'R';
                      im = uv.R.stack{I}  ; ac = uv.R.cart{  I};
                end
                submask.(LR){I} = false(size(im));
                figure,imshow(act_temp.*sf); colormap(par); tbx(gca); % REFERENCE

                % While loop gives the user a chance to approve each correction
                w = false;
                while ~w
                    %%

                    % Determine where to touch up and how many slices to draw
                    title('Draw a rectangle around the area you want to touch up');
                    h   = drawrectangle    ;
                    pos = round(h.Position);
                    imc = im(pos(2):pos(2)+pos(4)-1,pos(1):pos(1)+pos(3)-1,:);
                    sz  = size(imc); mask = false(sz);
                    idl = str2double(inputdlg(['How many slices do you want ',...
                                    'to outline?'],'# Slices',1,{'5'}));
                    sli   = round(1:(sz(1)-1)/(idl-1):sz(1));
                    imd   = imc(sli,:,:); imd(isnan(imd))= 0;
                    maskd = false(size(imd));
                    close(gcf);

                    % Loop each slice and manually draw exclusion zone
                    fig = figure;
                    clear possave
                    for kk = 1:idl
                        %%
                        w2 = true;
                        imshow(imrotate(squeeze(imd(kk,:,:)),-90),[]); tbx(gca);
                        if exist('possave','var'), fig.Position = possave; end
                        while w2
                            title('Outline the exclusion zone');
                            h = drawpolygon;
                            title('Mouse to keep, key to redo');
                            w2 = waitforbuttonpress;
                            if w2, delete(h); end
                        end
                        maskd(kk,:,:) = imrotate(h.createMask,90);
                        if kk > 1
                            mask(sli(kk-1):sli(kk),:,:) = permute(morph_binary(...
                               squeeze(maskd(kk-1,:,:)),squeeze(maskd(kk,:,:)),...
                               sli(kk)-sli(kk-1)-1),[3,1,2]);
                        end
                        possave = fig.Position;
                    end
                    close(gcf);

                    % Confirm successful edits and if editing is complete
                    submask_temp = submask.(LR){I};
                    submask_temp(pos(2):pos(2)+pos(4)-1,pos(1):pos(1)+pos(3)-1,:) = mask;
                    figure,imshowpair(act_temp .* sf,sum(ac & ~submask_temp,3) .* vox .* sf,...
                   'Scaling','none','Method','Montage'); colormap(par);
                    tbx(gca); title('Mouse to keep changes, key to discard');
                    w2 = waitforbuttonpress;
                    if ~w2, submask.(LR){I} = submask_temp;
                        act_temp = sum(ac & ~submask_temp,3) .* vox;
                    end
                    imshow(act_temp.*sf); colormap(par); tbx(gca);
                    title('Mouse to continue, key to end');
                    w = waitforbuttonpress;

                end
                close(gcf);            
            end

            % Ask the user if they are happy with their touch-ups, if they
            % want to retry, or if they have decided a rescan is necessary
            qdl = questdlg('Was touch-up successful?','Success!?',...
                'All Good','Retry','Rescan','Retry');
        end
    end
    
    % Check for 'rescan' response and act accordingly (there are two
    % opportunities to select rescan - this is why this code is not just
    % part of the above SWITCH)
    if strcmp(qdl,'Rescan')
        if exist('rescan.mat','file')
              rescan = load('rescan.mat'); rescan = rescan.rescan;
              rescan = unique([rescan;{d{ii}(1:end-4)}]);
              save('rescan.mat','rescan');
        else, rescan = {d{ii}(1:end-4)}; save('rescan.mat','rescan');
        end
    end
    
    % Save submask in UVH folder
    temp = submask;
    if ~isempty(temp)
        if isfield(temp,'L') && any(~cellfun(@isempty,temp.L))
            submask = temp.L;
            save([p4,d{ii}(1:2),'L',d{ii}(3:end)],'submask','-append');
        end
        if isfield(temp,'R') && any(~cellfun(@isempty,temp.R))
            submask = temp.R;
            save([p4,d{ii}(1:2),'R',d{ii}(3:end)],'submask','-append');
        end
    else
        save([p4,d{ii}(1:2),'L',d{ii}(3:end)],'submask','-append');
        save([p4,d{ii}(1:2),'R',d{ii}(3:end)],'submask','-append');
    end
    
end