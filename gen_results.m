%%
pI = [pwd,'\..\1 - Images\'  ]  ;
pA = [pwd,'\..\4 - ACT Analysis\5 - ACT Maps\']       ;
p3 = '3 - Cleanups\'            ;
d  = dir2cell(dir([p3,'*.mat']));
L = length(d);
f = @(x) cellfun(@(x)x(1:4),x,'UniformOutput',false);
results = cell(L+1,9);
results(1,:) = {'ID','AC.Th MFC' ,'AC.Th LFC' ,'AC.Sa MFC'  ,'AC.Sa LFC',...
                     'NIR MFC Av','NIR LFC Av','NIR MFC Dev','NIR LFC Dev'};
results(2:end,1) = f(d);
cond = {'MFC','LFC'};
%%
for ii = 1:L
    %%
    switch d{ii}(4)
        case 'F', fname = 'FemCond'; %mn = 0; mx = 3;
        case 'T', fname = 'Tibia'  ; %mn = 0; mx = 5;
    end
    t = tic;
    C   = load([p3,d{ii}],'Mkeep','ACTh','nir_interp');
    nir = load([pI,d{ii}],'nir'); nir = nir.nir.(fname).post.x800;
    ACT = load([pA,d{ii}(1:2),d{ii}(4)],'ACT'); ACT = ACT.ACT.(d{ii}(3));
    figure;
    res_row = cell(1,8);
    nir_interp = cell(1,2); ACTh = nir_interp;
    %
    for jj = 1:2
        %%
%         % Regenerate data
%         nir_interp{jj} = interp2(nir    ,Mkeep(jj).v2(:,1),Mkeep(jj).v2(:,2),'Linear');
%         ACTh{      jj} = interp2(ACT{jj},Mkeep(jj).v2(:,1),Mkeep(jj).v2(:,2),'Linear');
        I = C.ACTh{jj}>=.012; % Less than 1 voxel
        x = C.ACTh{jj}(I); x = x+rand(size(x)).*.001;
        y = C.nir_interp{jj}(I); y = y+rand(size(y)).*.001;
        subplot(1,2,jj), scatplot(x,y,'voronoi',[],100,5,1,12);
        axis square; title(cond{jj});
        AC_Mean  = mean(    C.ACTh{      jj}(I)            );
        AC_Dev   = mean(abs(C.ACTh{      jj}(I) - AC_Mean ));
        NIR_Mean = mean(    C.nir_interp{jj}               );
        NIR_Dev  = mean(abs(C.nir_interp{jj}    - NIR_Mean));
        res_row([1,3,5,7]+(jj-1)) = {AC_Mean,AC_Dev,NIR_Mean,NIR_Dev};
    end
    results(ii+1,2:end) = res_row;
    print([p3,'Previews\Graphs\',results{ii+1,1},'.tif'],'-dtiff','-r200');
    close(gcf);
    %
    % Overlay
    fMM = @(x,mn,mx) (x-mn)/(mx-mn);
    figure,imshow(nir,[]); hold on
    mn = min([C.nir_interp{1};C.nir_interp{2}]); mn = 2;
    mx = max([C.nir_interp{1};C.nir_interp{2}]); mx = 5.5;
    for jj = 1:2
       [RGB(:,1),RGB(:,2),RGB(:,3)] = ind2rgb(uint8(...
           fMM(C.nir_interp{jj},mn,mx)*255),plasma);
        patch(gca,'Vertices',C.Mkeep(jj).v2(:,1:2),'Faces',C.Mkeep(jj).f,...
        'EdgeColor','none','FaceVertexCData',RGB,'FaceColor','interp');
        clear RGB
    end
    print([p3,'Previews\Overlay\nir_',results{ii+1,1},...
        '_',num2str(mn),'_',num2str(mx),'.tif'],'-dtiff','-r200');
    close(gcf);
    
%     figure,imshow(nir,[]); hold on
%     mn = 0;
%     mx = 1.5;
%     for jj = 1:2
%        [RGB(:,1),RGB(:,2),RGB(:,3)] = ind2rgb(uint8(...
%            fMM(C.ACTh{jj},mn,mx)*255),parula);
%         patch(gca,'Vertices',C.Mkeep(jj).v2(:,1:2),'Faces',C.Mkeep(jj).f,...
%         'EdgeColor','none','FaceVertexCData',RGB,'FaceColor','interp');
%         clear RGB
%     end
%     print([p3,'Previews\Overlay\uct_',results{ii+1,1},...
%         '_',num2str(mn),'_',num2str(mx),'.tif'],'-dtiff','-r200');
%     close(gcf);
    
    looptrack(ii,L,t,results{ii+1,1});
end
% save('results.mat','results');