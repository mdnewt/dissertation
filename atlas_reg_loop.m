 %% Atlas Registration loop
 
 % Onlynew
 onlynew = true;
 
 % Create directory of reoriented scans
 p1   =  '1 - Images\'             ;
 p2   =  '2 - Atlas Registrations\';

 % In-line function
 f1 = @(x,y) cellfun(@(x)  x(y),x,'UniformOutput',false             );
 f2 = @(x  ) cellfun(@(x) ~isempty(who('-file',[p1,x],'ceuct_re')),x);

 % Load in atlas
 avg =  load('avgatlas.mat');
 aN  = {'fem','tib'}        ;

 id = affineTF_Maker([0,0,0],[1,1,1],[0,0,0]);
[opts,metric] = imregconfig('monomodal')     ;
 opts.MaximumStepLength = .03125             ;
 opts.MaximumIterations =  500               ;
 
 % Which samples need to be processed?
 dR  = dir2cell(dir([p1,'*R*.mat'])); dR = dR(f2(dR));
 dL  = dir2cell(dir([p1,'*L*.mat'])); dL = dL(f2(dL));
 dC  = dir2cell(dir([p2,  '*.mat']));
 dLR = dR(contains(f1(dR,[1,2,4]),f1(dL,[1,2,4]))); 
 if onlynew, d = dLR(~contains(f1(dLR,[1,2,4]),f1(dC,1:3)));
 else      , d = dLR;
 end
 clear dR dL dLR dC
 LL =  length(d);
    
 %%
 for ii = 1:LL % Loop of samples
    %%
    t = tic;

    % Femur or tibia?
    A = find(contains({'F','T'},d{ii}(4)));
    
    % Load in reoriented stacks
    Ln = d{ii}; Ln(3) = 'L';
    R  = load([p1,d{ii}],'ceuct_re');
    L  = load([p1,Ln   ],'ceuct_re');
    if ~isfield(R.ceuct_re,'vox') % Temp fix
        ceuct_re = R.ceuct_re; ceuct_re.vox = .012; R.ceuct_re.vox = .012;
        save([p1,d{ii}],'ceuct_re','-append'); clear ceuct_re
    end
    if ~isfield(L.ceuct_re,'vox') % Temp fix
        ceuct_re = L.ceuct_re; ceuct_re.vox = .012; L.ceuct_re.vox = .012;
        save([p1,Ln],'ceuct_re','-append'); clear ceuct_re
    end
    
    % Restructure into the format used by this code
    R = struct('stack_re',R.ceuct_re.stack_re,'R_re',R.ceuct_re.R_re,...
       'vox',R.ceuct_re.vox);
    L = struct('stack_re',L.ceuct_re.stack_re,'R_re',L.ceuct_re.R_re,...
       'vox',L.ceuct_re.vox);

    %%
    % Flip the left (the atlas is right-sided)
    L.stack_re = flip(L.stack_re,2);

    % Adjust stack to resolution of atlas
    R.sz_sm = round(size(R.stack_re) .* (R.vox/avg.vox));
    L.sz_sm = round(size(L.stack_re) .* (L.vox/avg.vox));
    R.R_sm = R.R_re; R.R_sm.ImageSize = R.sz_sm;
    L.R_sm = L.R_re; L.R_sm.ImageSize = L.sz_sm;
    R.stack_sm = imwarp(R.stack_re,R.R_re,id,'Linear','OutputView',R.R_sm);
    L.stack_sm = imwarp(L.stack_re,L.R_re,id,'Linear','OutputView',L.R_sm);

    % Crop stack down based on tissue thresholding
    fill = struct('R',imopen(imclose(R.stack_sm > 500,ones(9,9,9)),ones(9,9,9)),...
                  'L',imopen(imclose(L.stack_sm > 500,ones(9,9,9)),ones(9,9,9)));
   [~,~,tempR] = bwtrim(fill.R,15,R.R_sm);
   [~,~,tempL] = bwtrim(fill.L,15,L.R_sm);
    R.stack_sm = imwarp(R.stack_sm,R.R_sm,...
            affineTF_Maker([0,0,0],[1,1,1],[0,0,0]),'OutputView',tempR);
    L.stack_sm = imwarp(L.stack_sm,L.R_sm,...
            affineTF_Maker([0,0,0],[1,1,1],[0,0,0]),'OutputView',tempL);
    R.R_sm = tempR; L.R_sm = tempL; clear tempR tempL fill
    
    % Define moving and fixed images and respective reference frames
    mov = avg.(aN{A}).stack; rM   = avg.(aN{A}).R;
    fix = struct('R',R.stack_sm,'L',L.stack_sm);
    rF  = struct('R',R.R_sm    ,'L',L.R_sm    );
    %%
    % Affine transformation (using center alignment for initial guess)
    clear TF reg
   [tx,ty,tz] = align_centers(mov,rM,fix.L,rF.L);
    iTF       = affine3d(); iTF.T(end,1:3) = [tx,ty,tz];
     TF.R     = imregtform(mov,rM,fix.R,rF.R,'affine',opts,metric,...
               'InitialTransformation',iTF);
    reg.R     = imwarp(mov,rM,TF.R,'OutputView',rF.R);
   [tx,ty,tz] = align_centers(mov,rM,fix.L,rF.L);
    iTF       = affine3d(); iTF.T(end,1:3) = [tx,ty,tz];
     TF.L     = imregtform(mov,rM,fix.L,rF.L,'affine',opts,metric,...
               'InitialTransformation',iTF);
    reg.L     = imwarp(mov,rM,TF.L,'OutputView',rF.L);
    %%
    % NRR registrations
    clear df nrr
   [df.L,nrr.L] = imregdemons(reg.L,L.stack_sm);
   [df.R,nrr.R] = imregdemons(reg.R,R.stack_sm);

    R.I = quick_prev(R.stack_sm);
    L.I = quick_prev(L.stack_sm);

    % Export preview
    figure;
    imshowpair(R.stack_sm(:,:,R.I),nrr.R(:,:,R.I));
    title([d{ii}(1:2),'R',d{ii}(4)]);
    print([p2,'Previews\',d{ii}(1:2),'R',d{ii}(4),'.tif'],'-dtiff','-r150');
    imshowpair(L.stack_sm(:,:,L.I),nrr.L(:,:,L.I));
    title([d{ii}(1:2),'L',d{ii}(4)]);
    print([p2,'Previews\',d{ii}(1:2),'L',d{ii}(4),'.tif'],'-dtiff','-r150');
    close(gcf);

    % Save
    sm = struct; sm.R.stack = R.stack_sm; sm.R.R = R.R_sm;
                 sm.L.stack = L.stack_sm; sm.L.R = L.R_sm;
    save([p2,d{ii}([1:2,4:end])],'reg','TF','df','nrr','sm','-v7.3');
    looptrack(ii,LL,t,d{ii});

 end
 fprintf(['Processed ',num2str(LL),' scans.\n']);