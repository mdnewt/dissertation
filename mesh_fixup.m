onlynew = true;
pI = [pwd,'\..\1 - Images\'];
p2 = '2 - Mesh Mapping\'    ;
p3 = '3 - Cleanups\'        ;         

d  = dir2cell(dir([p2,'*.mat']));
if onlynew
    d2 = dir2cell(dir([p3,'*.mat']));
    d  = d(~contains(d,d2));
end
L = length(d);
% CURRENTLY CONDYLES ONLY
%%
for ii = 1:L
    %%
    t = tic;
    nir = load([pI,d{ii}],'nir'); nir = nir.nir;
    M   = load([p2,d{ii}]); nir_interp = M.nir_interp; Mvis = M.Mvis;
    ACTh = M.ACTh;
    
    switch d{ii}(4)
        case 'F', fname = 'FemCond';
        case 'T', fname = 'Tibia';
    end
    figure,I = imshow(nir.(fname).post.x800,[]); hold on;
    for jj = 1:2
        Mvis(jj).b = compute_boundary(Mvis(jj).f);
        plot(Mvis(jj).v2(Mvis(jj).b,1),Mvis(jj).v2(Mvis(jj).b,2));
    end
    w = true; submask = false(size(nir.(fname).post.x800));
    while w
        w2 = true;
        while w2
            h = impoly();
            title('Mouse to keep, key to redo');
            w2 = waitforbuttonpress;
            if ~w2, sub = h.createMask(I); end
            delete(h);
        end
        submask = submask | sub;
        IMS = imshow(cat(3,single(sub),zeros(size(sub)),zeros(size(sub))));
        IMS.AlphaData = sub .* 0.3;
        title('Mouse to continue, key to end');
        w = ~waitforbuttonpress;
    end
    close(gcf);
    %
end
%%
for ii = 1:L
    %%
    t = tic;
    switch d{ii}(4)
        case 'F', fname = 'FemCond';
        case 'T', fname = 'Tibia';
    end
    nir = load([pI,d{ii}],'nir'); nir = nir.nir;
    M   = load([p2,d{ii}]); nir_interp = M.nir_interp; Mvis = M.Mvis;
    ACTh = M.ACTh;
    submask = load([p3,d{ii}],'submask'); submask = submask.submask;
    clear Mkeep vertroi vtest ftest
    for jj = 1:2
        vtest{jj} = ~interp2(single(submask),Mvis(jj).v2(:,1),Mvis(jj).v2(:,2),'Nearest');
        ftest{jj} = all(ismember(Mvis(jj).f,find(vtest{jj})),2);
       [Mkeep(jj),vertroi{jj}] = subMesh(Mvis(jj),find(ftest{jj}));
    end
    for jj = 1:2
        Mkeep(jj).b  = compute_boundary(Mkeep(jj).f);
        Mkeep(jj).v2 = Mvis(jj).v2(vertroi{jj}(:,2),:);
        nir_interp{jj} = M.nir_interp{jj}(vertroi{jj}(:,2));
        ACTh{jj} = M.ACTh{jj}(vertroi{jj}(:,2));
    end
    figure,I = imshow(nir.(fname).post.x800,[]); hold on
    for jj = 1:2
        plot(Mkeep(jj).v2(Mkeep(jj).b,1),Mkeep(jj).v2(Mkeep(jj).b,2));
    end
    print([p3,'Previews\Boundary\',d{ii}(1:end-3),'tif'],'-dtiff','-r200');
    close(gcf);
    save([p3,d{ii}],'Mkeep','ACTh','nir_interp','vtest','ftest','submask');
    looptrack(ii,L,t);
end