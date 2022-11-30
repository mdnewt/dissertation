%% Create 'figgen.mat' containing worst DSCs
% p_E = ['S:\OrthopaedicsProjects\Newton\Automated Segmentation\CE-uCT ',...
%        'Knee - Atlas Registration\Validation\Manuscript\ISA Data ',...
%        'Compilation 1-15-19.xlsx'];
% groups = {1:10,11:20,21:30,31:40,41:50};
% figgen          = cell(6,13);
% figgen(1,:    ) = {'# Inputs','Worst Epi it','Worst Epi S','Worst Epi DSC',...
%                               'Best Epi it' ,'Best Epi S' ,'Best Epi DSC' ,...
%                               'Worst AC it' ,'Worst AC S' ,'Worst AC DSC' ,...
%                               'Best AC it'  ,'Best AC S'  ,'Best AC DSC'};
% figgen(2:end,1) = {1;3;5;10;20};
% for ii = 2:5
%     nums = [];
%     txts = cell(0,2);
%     for jj = 1:10
%         t = tic;
%        [temp1,temp2] = xlsread(p_E,[num2str(groups{ii}(jj)),'.xlsx'],'A:Z');
%         nums = [nums;temp1(:,[6,25])];
%         txts = [txts;repmat({num2str(groups{ii}(jj))},[size(temp1,1),1]),...
%                 temp2(3:end,1)];
%         looptrack(groups{ii}(jj),50,t);
%     end
%    [mn,mnI] = min(nums,[],1);
%    [mx,mxI] = max(nums,[],1);
%     figgen(ii+1,2:end) = [txts(mnI(1),:),{mn(1)},...
%                           txts(mxI(1),:),{mx(1)},...
%                           txts(mnI(2),:),{mn(2)},...
%                           txts(mxI(2),:),{mx(2)}];
% end
% save('figgen.mat','figgen');
% clear temp1 temp2 nums txts mn mnI groups p_E t

%% Generate HD figures

% % Setup
% load('figgen.mat')     ; load('vox.mat');
% vox_train = 0.012      ;
% stan_TH   = 3800       ;
% ACLR2W_TH = 5000       ;
% p_A       = 'Analysis\';
% 
% for jj = 1:4
% for ii = 1:size(figgen,1)-1
%     
%     p_fN  = [p_A,figgen{ii+1,3*jj-1},'\'];
%     d_Reg =  dir([p_fN,'2 - Registrations\*.mat']);
%     if ~isempty(d_Reg)
%         
%         % Bone
%         if any(jj==[1,2])
%             
%             % Load stuff
%             HD = load([p_fN,'3 - Extras\',figgen{ii+1,3*jj}],...
%                 'bone','HD_bone','bone_fill','HD_fill');
%             
%             % Left to right mirroring
%             if contains(figgen{ii+1,3*jj},'L')
%                 HD.bone.fix      = flip(HD.bone.fix     ,2);
%                 HD.HD_bone.dist  = flip(HD.HD_bone.dist ,2);
%                 HD.bone_fill.fix = flip(HD.bone_fill.fix,2);
%                 HD.HD_fill.dist  = flip(HD.HD_fill.dist ,2);
%             end
%             
%             FV   = isosurface(bwareaopen(HD.bone.fix,100));
%             cdat = interp3(HD.HD_bone.dist .* vox_unk,...
%                    FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));
%             figure,P = patch(FV,'EdgeColor','none','FaceColor','interp','CData',cdat);
%             axis equal; view([-155,40]); colorbar;
% %             view([0,-90]);
% 
%             FV2   = isosurface(bwareaopen(HD.bone_fill.fix,100));
%             cdat2 = interp3(HD.HD_fill.dist .* vox_unk,...
%                    FV2.vertices(:,1),FV2.vertices(:,2),FV2.vertices(:,3));
%             figure,P2 = patch(FV2,'EdgeColor','none','FaceColor','interp','CData',cdat2);
%             axis equal; view([0,-90]); colorbar;
% 
% %             savefig(['S:\OrthopaedicsProjects\Newton\Automated Segmentation\',...
% %                      'CE-uCT Knee - Atlas Registration\Validation\Manuscript',...
% %                      '\Figure\HD\',figgen{1,3*jj}(1:end-2),'_',...
% %                      num2str(figgen{ii+1,1}),'_fill.fig']);
%             savefig([figgen{1,3*jj}(1:end-2),'_',...
%                      num2str(figgen{ii+1,1}),'_fill.fig']);
%             close(gcf);
%             
%         % AC
%         elseif any(jj==[3,4])
%             
%             % Load stuff
%             HD = load([p_fN,'3 - Extras\',figgen{ii+1,3*jj}],'cart','HD_AC');
%             
%             % Left to right mirroring
%             if contains(figgen{ii+1,3*jj},'L')
%                 HD.cart.fix.mask = flip(HD.cart.fix.mask,2);
%                 HD.HD_AC.dist    = flip(HD.HD_AC.dist   ,2);
%             end
%             
%             % Create patches and save figures
%             FV   = isosurface(HD.cart.fix.mask);
%             cdat = interp3(HD.HD_AC.dist .* vox_unk,...
%                    FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));
%             figure,P = patch(FV,'EdgeColor','none','FaceColor','interp','CData',cdat);
%             axis equal; view([-135,20]); colorbar;
% %             view([0,-90]);
%             
%         else, error('Sum Ting Wong.');
%         end
%         
% %         savefig(['S:\OrthopaedicsProjects\Newton\Automated Segmentation\',...
% %                  'CE-uCT Knee - Atlas Registration\Validation\Manuscript',...
% %                  '\Figure\HD\',figgen{1,3*jj}(1:end-2),'_',...
% %                  num2str(figgen{ii+1,1}),'.fig']);
%         savefig([figgen{1,3*jj}(1:end-2),'_',...
%                  num2str(figgen{ii+1,1}),'.fig']);
%         close(gcf)
%     end
%     looptrack(ii,jj);
% end
% end

% %% Save AC images
% p1 = ['C:\Users\127907\Downloads\hausdorff output for newton-',...
%       '20190816T123539Z-001\hausdorff output for newton\'];
% p2 = ['S:\OrthopaedicsProjects\Newton\Automated Segmentation\',...
%       'CE-uCT Knee - Atlas Registration\Validation\Manuscript\Figure\HD\'];
% p3 = 'C:\Users\127907\Desktop\HD Figs\';
% 
% d1 = dir2cell(dir([p1,'*AC*.fig']));
% d2 = dir2cell(dir([p2,'*AC*.fig']));
% 
% for ii = 1:length(d1)
%     t = tic;
%     name = strsplit(d1{ii},'.'); name = name{1};
%     openfig([p1,d1{ii}]);
%     caxis([0,0.7]);
%     print([p3,name,'_1.eps'],'-depsc2');
%     view([0,-90]);
%     print([p3,name,'_2.eps'],'-depsc2');
%     looptrack(ii,length(d1),t,name);
%     close(gcf);
% end
% for ii = 1:length(d2)
%     t = tic;
%     name = strsplit(d2{ii},'.'); name = name{1};
%     openfig([p2,d2{ii}]);
%     caxis([0,0.7]);
%     print([p3,name,'_1.eps'],'-depsc2');
%     view([0,-90]);
%     print([p3,name,'_2.eps'],'-depsc2');
%     looptrack(ii,length(d2),t,name);
%     close(gcf);
% end

%% Save Epi images
p1 = ['C:\Users\127907\Downloads\hausdorff output for newton-',...
      '20190816T123539Z-001\hausdorff output for newton\'];
p2 = ['S:\OrthopaedicsProjects\Newton\Automated Segmentation\',...
      'CE-uCT Knee - Atlas Registration\Validation\Manuscript\Figure\HD\'];
p3 = 'C:\Users\127907\Desktop\HD Figs\';

d1 = dir2cell(dir([p1,'*Epi*.fig'])); d1 = d1(1:2:end);
d2 = dir2cell(dir([p2,'*Epi*.fig'])); d2 = d2(1:2:end);

% for ii = 1:length(d1)
%     t = tic;
%     name = strsplit(d1{ii},'.'); name = name{1};
%     openfig([p1,d1{ii}]);
%     caxis([0,0.7]);
%     print([p3,name,'_1.eps'],'-depsc2');
%     view([0,-60]);
%     print([p3,name,'_2.eps'],'-depsc2');
%     looptrack(ii,length(d1),t,name);
%     close(gcf);
% end
for ii = 1:length(d2)
    t = tic;
    name = strsplit(d2{ii},'.'); name = name{1};
    openfig([p2,d2{ii}]);
    caxis([0,0.7]);
%     print([p3,name,'_1.eps'],'-depsc2');
    view([0,-60]);
    print([p3,name,'_2.eps'],'-depsc2');
    looptrack(ii,length(d2),t,name);
    close(gcf);
end