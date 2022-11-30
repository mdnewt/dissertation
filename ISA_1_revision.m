%% Startup

% Load in pre-requisite variables
load('revision.mat','whichiswhich','d','L');
d{1} = d{1}(:,1); d{2} = d{2}(:,1); d_save = d;

% Generate directory and name list
p    = 'Stacks\'      ;
f    = @(x) x(1:end-4);
name = {cellfun(f,d{1},'UniformOutput',false),...
        cellfun(f,d{2},'UniformOutput',false)}; clear f

% Determine voxel information
vox = load('vox.mat','vox_train');
vox = vox.vox_train;

% Configure options for affine registration
[opts,met]              = imregconfig('monomodal');
 opts.MaximumIterations = 500                     ;

CI = 'CI';

%%
for cc = 1  % Loop 1 - Control or Injury
for at = 17:25 % Loop 2 - Atlas Iteration #
%%

% STEP 0: SET UP NEW ATLAS ITERATION

% Generate filename
fn = [CI(cc),num2str(at,'%02.f')];

% Apply randomizer for this iteration
d = d_save                       ;
d = d{cc}(whichiswhich{cc}{at,1});
L = length(d)                    ;

% Load and prepare reference image (first on the list)
fix = load([p,d{1}]);
if contains(d{1},'L'), fix.stack = flip(fix.stack,2); end % L/R flip
fix.stack_full = fix.stack                                       ;
fix.stack      = stackinterp(fix.stack,round(size(fix.stack)./2));
fix.R          = imref3d(size(fix.stack),vox,vox,vox)            ;

% For 1-sample atlaes
if L==1
    avg = struct('Image',fix.stack);
    save(['Analysis - Revision\1 - Average Atlas\Atlas\',fn,'.mat'],'avg','-v7.3');
    fprintf(['Atlas iteration ',fn,' complete\n']);
    continue
end

% Create deformation and registration placeholder cells
deforms      = cell(L+1,2);
deforms{1,1} = 'Samples'; deforms(2:end,1) = d;
deforms{1,2} = 'Affine' ;
reg          = deforms  ;
reg{2,2}     = fix.stack;

% Register the other stacks onto the reference
for ii = 2:L % Loop 3.1 - AFFINE REGISTRATION
    %%
    
    t = tic;
    
    % Load and prepare moving sample
    mov = load([p,d{ii}],'stack');
    if contains(d{ii},'L'), mov.stack = flip(mov.stack, 2); end % L/R flip
    mov.stack_full = mov.stack                                       ;
    mov.stack      = stackinterp(mov.stack,round(size(mov.stack)./2));
    mov.R          = imref3d(size(mov.stack),vox,vox,vox)            ;

    % Perform the registration
    TF              = imregtform(mov.stack,mov.R,fix.stack,fix.R,'affine',opts,met);
    reg{ii+1,2}     = imwarp(    mov.stack,mov.R,TF,'linear','OutputView',fix.R   );
    deforms{ii+1,2} = TF;
    
    looptrack(ii-1,L-1,t,'Rigid');
end

% Create 'iteration 0' average atlas
avg_image = {mean(cat(4,reg{2:end,2}),4)};

% Save iteration 0
save(['Analysis - Revision\1 - Average Atlas\Registrations\',fn,'.mat'],...
    'reg','deforms','avg_image', '-v7.3');

fprintf('Iteration 0 complete.\n');

it = 1; avgdef = [];
while (it < 3 || diff(avgdef(end:-1:end-1)) > .1) && it < 10 % Loop 3.2 - NRR 'outer' iterations
    %%
    
    % Setup new row of deforms
    deforms{1    ,it+2} = ['NRR It. ',num2str(it)];
    deforms(2:end,it+2) = {[]}                    ;
    reg(:,it+2) = deforms(:,it+2);
    
    % Reference image
    fix = avg_image{it};
    
    for ii = 1:L % Loop 4.2 - Discrete NRR registrations
        %%
        t = tic;

        % Perform the registration
        mov                                = reg{ii+1,it+1}      ;
       [deforms{ii+1,it+2},reg{ii+1,it+2}] = imregdemons(mov,fix);

        looptrack(ii,length(d),t,['It. ',num2str(it)]);

    end
    
    % Create average image and average deformation
    avg_image{it+1} = mean(cat(4,reg{2:end,it+2}),4)                      ;
    ad              = sqrt(sum((mean(cat(5,deforms{2:end,it+2}),5)).^2,4));
    avgdef(it)      = mean(ad(:))                                         ;
    
    % Save iteration
    save(['Analysis - Revision\1 - Average Atlas\Registrations\',fn,'.mat'],...
          'reg','deforms','avg_image','avgdef','-append');

    fprintf(['NRR iteration ',num2str(it),' complete.\n']);
    it = it+1;
    
end

% Save final average image
avg = struct('Image',avg_image{end});
save(['Analysis - Revision\1 - Average Atlas\Atlas\',fn,'.mat'],'avg','-v7.3');

fprintf(['Atlas iteration ',fn,' complete\n']);

end
end