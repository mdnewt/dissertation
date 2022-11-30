% Batch reorientation of .mat image stacks

% New samples only? Or all?
onlynew = true;

p = '1 - Images\'              ;
d =  dir2cell(dir([p,'*.mat']));

% Cull samples that have already been processed
if onlynew
    f1 = @(x) isempty(who('-file',[p,x],'TF'));
    I1 = cellfun(f1,d);
    d  = d(I1); L = length(d); clear f1 I1
else,           L = length(d);
end

%% Loop of reorientations
for ii = 1:L
    %% Perform the manual step
    t = tic;
    m = matfile([p,d{ii}],'Writable',true);
    ceuct = m.ceuct;
   [~,TF] = imstack_reorient(ceuct.stack,[],[],[],true);
    m.TF  = TF;
%     save([p,d{jj}],'TF','-append');
    looptrack(ii,L,t,d{ii}(1:end-4));
%     if ii ~= L
%         q = questdlg([num2str(ii),' / ',num2str(L),' ',...
%             bone{ii}(1:3),'s completed. ','Proceed?'],'Yes','No');
%         switch q, case 'No', break; end
%     end
end
fprintf(['Manual step completed on ',num2str(ii),' scans.\n']);

%% Second half of loop    
% Reassess which scans require manual reorientation
d  =  dir2cell(dir([p,'*.mat']));
f1 = @(x) isempty(who('-file',[p,x],'TF'));
I1 = cellfun(f1,d);
if onlynew
      f2 = @(x) isempty(who('-file',[p,x],'ceuct_re'));
      I2 = cellfun(f2,d);
      d  = d(~I1 & I2); L = length(d);
else, d  = d(~I1     ); L = length(d);
end
%%
for ii = 1:L
    %%
    t = tic;
    m = matfile([p,d{ii}],'Writable',true);
    ceuct = m.ceuct;
   [stack_re,~,R_re,R] = imstack_reorient(ceuct.stack,[],[],m.TF,false,ceuct.vox);
    ceuct_re   = struct('stack_re',stack_re,'R_re',R_re,'R',R,'vox',ceuct.vox);
    m.ceuct_re = ceuct_re;
%     save([p,d{ii}],'stack_re','R_re','R','-append');
    looptrack(ii,L,t,d{ii}(1:end-4));
end
fprintf(['Successfully reoriented ',num2str(L),' scans.\n']);