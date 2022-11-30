% New samples only? Or all?
onlynew = true;

% Create directory of reoriented scans
p1   =  '1 - Images\'             ;
p2   =  '2 - Atlas Registrations\';
p3   =  '3 - Masks\'              ;
bone = {'Femur\','Tibia\'}        ;

% In-line function
f1 = @(x,y) cellfun(@(x)  x(y),x,'UniformOutput',false             );
f2 = @(x  ) cellfun(@(x) ~isempty(who('-file',[p1,x],'ceuct_re')),x);

% Load atlas
avg =  load('avgatlas.mat');
aN  = {'fem','tib'}        ;
F   = {'cartroi','epi'}    ;
C   = {'med','lat'}        ;

LR = 'LR';
id = affineTF_Maker([0,0,0],[1,1,1],[0,0,0]);

%%
    
d2 = dir2cell(dir([p2,'*.mat']));
d3 = dir2cell(dir([p3,'*.mat']));
if onlynew, d = d2(~contains(d2,d3));
else,       d = d2;
end
LL = length(d); clear d2 d3
    
%% Step 1: Create the reoriented masks
for ii = 1:LL
    %%

    t = tic;
    
    % Femur or tibia?
    A = find(contains({'F','T'},d{ii}(3)));
    
    % Load the things
    R   = load([p1,d{ii}(1:2),'R',d{ii}(3:end)],'ceuct','ceuct_re','TF');
    L   = load([p1,d{ii}(1:2),'L',d{ii}(3:end)],'ceuct','ceuct_re','TF');
    reg = load([p2,d{ii}]);
    nrr = reg.nrr; sm = reg.sm; TF = reg.TF; df = reg.df; reg = reg.reg;

    % Generate full resolution masks that match to stack_re
    mask = struct; mask_re = mask;
    R.iTF = R.TF; R.iTF.T = inv(R.iTF.T); 
    L.iTF = L.TF; L.iTF.T = inv(L.iTF.T);
    
    %%
    for jj = 1:2
        %%

        % Compute the correct masks in lower resolution
        mask_re.(F{jj}).L = imwarp(avg.(aN{A}).(F{jj}),avg.(aN{A}).R,...
                            TF.L,'Linear','OutputView',sm.L.R);
        mask_re.(F{jj}).L = imwarp(mask_re.(F{jj}).L,df.L,'Nearest');
        mask_re.(C{jj}).L = imwarp(avg.(aN{A}).comps.(C{jj}),avg.(aN{A}).R,...
                            TF.L,'Linear','OutputView',sm.L.R);
        mask_re.(F{jj}).R = imwarp(avg.(aN{A}).(F{jj}),avg.(aN{A}).R,...
                            TF.R,'Linear','OutputView',sm.R.R);
        mask_re.(F{jj}).R = imwarp(mask_re.(F{jj}).R,df.R,'Nearest');
        mask_re.(C{jj}).R = imwarp(avg.(aN{A}).comps.(C{jj}),avg.(aN{A}).R,...
                            TF.R,'Linear','OutputView',sm.R.R);
%%
        % Upscale masks to full resolution
        mask_re.(F{jj}).L = imwarp(mask_re.(F{jj}).L,sm.L.R,id,'Linear',...
                           'OutputView',L.ceuct_re.R_re);
        mask_re.(C{jj}).L = imwarp(mask_re.(C{jj}).L,sm.L.R,id,'Linear',...
                           'OutputView',L.ceuct_re.R_re);
        mask_re.(F{jj}).L = flip(mask_re.(F{jj}).L,2);
        mask_re.(C{jj}).L = flip(mask_re.(C{jj}).L,2);
        mask_re.(F{jj}).R = imwarp(mask_re.(F{jj}).R,sm.R.R,id,'Linear',...
                           'OutputView',R.ceuct_re.R_re);
        mask_re.(C{jj}).R = imwarp(mask_re.(C{jj}).R,sm.R.R,id,'Linear',...
                           'OutputView',R.ceuct_re.R_re);

        % Derotate
        mask.(F{jj}).L = imwarp(mask_re.(F{jj}).L,L.ceuct_re.R_re,...
                         L.iTF ,'Linear','OutputView',L.ceuct_re.R) > 0.5;
        mask.(C{jj}).L = imwarp(mask_re.(C{jj}).L,L.ceuct_re.R_re,...
                         L.iTF ,'Linear','OutputView',L.ceuct_re.R) > 0.5;
        mask.(F{jj}).R = imwarp(mask_re.(F{jj}).R,R.ceuct_re.R_re,...
                         R.iTF ,'Linear','OutputView',R.ceuct_re.R) > 0.5;
        mask.(C{jj}).R = imwarp(mask_re.(C{jj}).R,R.ceuct_re.R_re,...
                         R.iTF ,'Linear','OutputView',R.ceuct_re.R) > 0.5;

    end
    %
    %%
    % Create extended troch mask using plane fitting (femur only)
    if contains(d{ii},'F')
        for jj = 1:2
            %%
            ml     =  mask.med.(LR(jj)) | mask.lat.(LR(jj));
            sz     =  size(ml)                             ;
            sli    =  round(sz(3) * [.25,.25,.5,.75,.75])  ;
            crpS   =  round(sz(2) * [.4,.6]     )          ;
            crpC   =  round(sz(1) * [.3,.7]     )          ;
            ml2    =  bwperim(ml(crpC(1):crpC(2),...
                                 crpS(1):crpS(2),sli),4)   ;
            ml2    =  ml2(2:end-1,2:end-1,2:end-1)         ;
            ml2    =  padarray(ml2,[1,1,0],'both')         ;
           [x,y,z] =  ind2sub(size(ml2),find(ml2(:,:,:)))  ;

            pts    = {x(z==1),y(z==1),x(z==2),y(z==2),x(z==3),y(z==3)};
            pts    = [pts{1}(10 )+crpC(1)-1,pts{2}(10 )+crpS(1)-1,sli(2);...
                      pts{3}(1  )+crpC(1)-1,pts{4}(1  )+crpS(1)-1,sli(3);...
                      pts{5}(end)+crpC(1)-1,pts{6}(end)+crpS(1)-1,sli(4)];
            v      = [pts(1,:)-pts(2,:);pts(1,:)-pts(3,:)];
            pln    =  cross(v(1,:),v(2,:)); pln = pln/norm(pln);
            k      =  sum(pln .* pts(1,:));
            f3     =  @(x,y,z) pln(1)*x + pln(2)*y + pln(3)*z - k > 0;
           [y,x,z] =  meshgrid(1:sz(2),1:sz(1),1:sz(3));
            tro   =   f3(x,y,z);
            if nnz(tro & ml) > nnz(~tro & ml), tro = ~tro; end
            mask.tro.(LR(jj)) = tro                      ;
        end
    end
    %%
    % Extend medial and lateral masks using plane fitting
    for jj = 1:2
        %%
        ml     =  bwperim(imdilate(mask.med.(LR(jj)),ones(3,3,3)) & ...
                          imdilate(mask.lat.(LR(jj)),ones(3,3,3)),4);
        sz     =  size(ml);
       [x,y,z] =  ind2sub(size(ml),find(ml));
        sli    =  round((max(z)-min(z)).*[.25,.5,.75] + min(z));
        pts    = {x(z==sli(1)),y(z==sli(1)),...
                  x(z==sli(2)),y(z==sli(2)),...
                  x(z==sli(3)),y(z==sli(3))};
        pts    = [pts{1}(10 ),pts{2}(10 ),sli(1);...
                  pts{3}(1  ),pts{4}(1  ),sli(2);...
                  pts{5}(end),pts{6}(end),sli(3)];
        v      = [pts(1,:)-pts(2,:);pts(1,:)-pts(3,:)];
        pln    =  cross(v(1,:),v(2,:)); pln = pln/norm(pln);
        k      =  sum(pln .* pts(1,:));
        f3     =  @(x,y,z) pln(1)*x + pln(2)*y + pln(3)*z - k > 0;
       [y,x,z] =  meshgrid(1:sz(2),1:sz(1),1:sz(3));
        med    =  imclose( f3(x,y,z),ones(5,5,5));
        lat    =  imclose(~f3(x,y,z),ones(5,5,5));
        if isfield(mask,'tro')
            med = med & ~mask.tro.(LR(jj));
            lat = lat & ~mask.tro.(LR(jj));
        end
        if nnz(med & mask.med.(LR(jj))) < nnz(~med & mask.med.(LR(jj)))
            temp = med; med = lat; lat = temp; 
        end
        mask.med.(LR(jj)) = med; mask.lat.(LR(jj)) = lat;
    end
    
%%
    % Create preview
    R.I = quick_prev(R.ceuct.stack);
    L.I = quick_prev(L.ceuct.stack);
    im  = struct('R',repmat(R.ceuct.stack(:,:,R.I),[1,2]),...
                 'L',L.ceuct.stack(:,:,L.I));
    rgb = struct;
    rgb.R = cat(3,mask.epi.R(:,:,R.I),mask.cartroi.R(:,:,R.I),...
                  mask.med.R(:,:,R.I));
    rgb.L = cat(3,mask.epi.L(:,:,L.I),mask.cartroi.L(:,:,L.I),...
                  mask.med.L(:,:,L.I));

    figure;
    i1 = imshow(L.ceuct.stack(:,:,L.I),[]); hold on; i2 = imshow(single(rgb.L));
    a  = (rgb.L(:,:,1) | rgb.L(:,:,2) | rgb.L(:,:,3)) .* 0.4; set(i2,'AlphaData',a);
    title('L');
    print([p3,'Previews\',d{ii}(1:2),'L',d{ii}(3),'.tif'],'-dtiff','-r150');
    hold off
    i1 = imshow(R.ceuct.stack(:,:,L.I),[]); hold on; i2 = imshow(single(rgb.R));
    a  = (rgb.R(:,:,1) | rgb.R(:,:,2) | rgb.R(:,:,3)) .* 0.4; set(i2,'AlphaData',a);
    title('R');
    print([p3,'Previews\',d{ii}(1:2),'R',d{ii}(3),'.tif'],'-dtiff','-r150');
    close(gcf);
%
    % Save
    save([p3,d{ii}],'-struct' , 'mask'  );
    save([p2,d{ii}], 'mask_re','-append');

    looptrack(ii,LL,t,d{ii}(1:3));

end
fprintf(['Processed ',num2str(LL),' scans.\n']);