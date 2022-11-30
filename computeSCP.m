function [Meshes,pts] = computeSCP(Mesh,comp)

% Create Meshes structure
F      = @(x) [upper(x(1)),x(2:end)]; % Uppercase fxn
Meshes = struct('name','','v',[],'uv',[],'f',[],'n',[]);
fnames = fieldnames(comp); LL = length(fnames);
isfem  = LL==3;
for ii = 1:LL, Meshes(ii).name = F(fnames{ii}); end
Meshes(end+1).name = 'Whole';
Meshes(end).f = Mesh.f; Meshes(end).v = Mesh.v; Meshes(end).n = Mesh.n;

% Parameterize whole surface (possible for femur only)

if isfem
    uv =  embedSCP(Mesh,'generalized'); uv(:,1) = -uv(:,1); % Fix flip
    uv = [uv,zeros(length(uv),1)]     ;
    Meshes(end).uv = uv; Meshes(4).u = uv(:,1:2);
end

% Automatically break into sections, based on compartment mask
F    = @(x,M) interp3(x,M(end).v(:,1),M(end).v(:,2),M(end).v(:,3),'Nearest');
idxV = cell(1,LL);
for ii = 1:LL
    idxV{ii} = F(single(comp.(fnames{ii})),Meshes);
    idxV{ii}(isnan(idxV{ii})) = false             ;
    idxV{ii} = logical(idxV{ii})                  ;
end
F    = @(x) find(sum(ismember(Mesh.f,find(x)),2)==3);
idxF = cellfun(F,idxV,'UniformOutput',false)        ;

for ii = 1:LL
    %%
    % Perform Mesh subdivision
    M = subMesh(Mesh,idxF{ii});
    M = clipEars(M)           ;

    % Parameterize
    uv =  embedSCP(M,'generalized'); uv(:,1) = -uv(:,1); % Fix flip
    uv = [uv,zeros(length(uv),1)]  ;

    % Assign final v, f, n, uv to Meshes structure
    Meshes(ii).f  = M.f; Meshes(ii).v = M.v;
    Meshes(ii).uv = uv ; Meshes(ii).n = M.n; Meshes(ii).u = uv(:,1:2);
end

% Reorient
Meshes = cart_reorient(Meshes,idxV,idxF);

% If tibia, create composite whole tissue parameterization
if ~isfem
    %%
    if contains(Mesh.name,'L'), II = [2,1]; else, II = 1:2; end
    uv1 = Meshes(II(1)).u; uv1 = uv1 .* 1/max(uv1(:,2));
    uv2 = Meshes(II(2)).u; uv2 = uv2 .* 1/max(uv2(:,2));
    uv2(:,1) = uv2(:,1) + max(uv1(:,1)) + 0.1;
    temp = zeros(length(Meshes(end).v),2);
    temp(idxV{II(1)},:) = uv1; temp(idxV{II(2)},:) = uv2;
    Meshes(end).u = temp; Meshes(end).uv = temp; Meshes(end).uv(:,3) = 0;
end

% figure;
% subplot(1,3,1),patch('Vertices',Meshes(1).uv,'Faces',Meshes(1).f,...
%     'FaceColor','r'); axis equal;
% subplot(1,3,2),patch('Vertices',Meshes(2).uv,'Faces',Meshes(2).f,...
%     'FaceColor','g'); axis equal
% subplot(1,3,3),patch('Vertices',Meshes(3).uv,'Faces',Meshes(3).f,...
%     'FaceColor','b'); axis equal