function [uv_stack,uv_cart,bwm] = xyz_to_uvh(stack,cart,Meshes,vox,thresh,op)
% 
% 
% Input variables:
% stack  - CE-µCT image stack (original, not reoriented)
% cart   - Cartilage ROI mask (including adjacent bone and air)
% Meshes - Compartmentalized mesh structure
% vox    - Voxel size of stack
% thresh - Cartilage tissue thresholds, supplied as [low,high]
% op     - Optional binary variable; if true, displays extra information

% Check for threshold
if ~exist('thresh','var'), thresh = [200,3500]; end

% Create storage variables
LM       = nnz(~cellfun(@isempty,cat(1,{Meshes(:).uv})));
uv_stack = cell(1,LM)    ;
uv_cart  = uv_stack      ;
bwm      = uv_stack      ;

for jj = 1:LM
    %%
    % Get current mesh
    v = Meshes(jj).v; uv = Meshes(jj).uv;
    n = Meshes(jj).n; f  = Meshes(jj).f ;

    % Ensure normals are normalized
    n = n./repmat(sqrt((n(:,1).^2+n(:,2).^2+n(:,3).^2)),[1,3]);

    % Ensure normals are pointing outwards based on scaling
    scale  = mean2(abs(v  -repmat(mean(v  ),[length(v),1])));
    scale2 = mean2(abs(v+n-repmat(mean(v+n),[length(v),1])));
    if scale > scale2, n = n.*-1; end

    % Arrange for a parameterized stack +1 & -0.5 mm from the BCI at the 
    % same resolution as the original image (supplied by 'vox')
    H = ceil(0.5/vox);
    Q = -H:1:H   ;

    % Interpolate intensity and mask values at each vertex
    I_data = zeros(length(v),length(Q));
    M_data = zeros(length(v),length(Q));

    for q = 1:length(Q)
        %%
        v_norm = v+n.*Q(q);
        I_data(:,q) = interp3(stack,v_norm(:,1),v_norm(:,2),v_norm(:,3));
%             M_data(:,q) = interp3(cart ,v_norm(:,1),v_norm(:,2),v_norm(:,3),...
%                                   'nearest');
        M_data(:,q) = round(interp3(cart ,v_norm(:,1),v_norm(:,2),v_norm(:,3)));
    end
    I_data = reshape(I_data,[],1);
    M_data = reshape(M_data,[],1);

    % Move min to 0
    mn      = min(uv,[],1) ;
    uv(:,1) = uv(:,1)-mn(1);
    uv(:,2) = uv(:,2)-mn(2);
    mx      = max(uv,[],1) ;

    % Compute mesh boundary
    b = compute_boundary(f); % Can result in error if there are holes in the mesh

    % Determine relative change in edge length v --> uv
    Mesh = makeMesh(v,f,n); e = Mesh.e;
    eXYZ = edgeStats(v,e); eUV = edgeStats(uv,e);

    % Create a matrix sized so that the effective 'voxel size' of the
    % parameterization is roughly equal to that of the original stack
    matSize = ceil([(eXYZ/eUV).*mx(1:2),length(Q)]);

    % Resize uv
    uv2       = zeros(size(uv))             ;
    uv2(:,1)  = uv(:,1).*(matSize(1)/mx(1)) ;
    uv2(:,2)  = uv(:,2).*(matSize(2)/mx(2)) ;
    uv2(:,1:2)= uv2(:,1:2) + 2              ;

    % Add height information
    uv_h    = repmat(uv2,[length(Q),1])  ;
    for kk = 1:length(Q)
        uv_h(1+(kk-1)*length(uv2):kk*length(uv2),3) = kk;
    end

    % Discretize
   [x,y,z]   = meshgrid(1:matSize(1)+3,1:matSize(2)+3,1:matSize(3));
    sInt     = scatteredInterpolant(uv_h(:,1),uv_h(:,2),uv_h(:,3),...
        I_data,'linear','none');
    temp_stack = sInt(x,y,z);
%         sInt     = scatteredInterpolant(uv_h(:,1),uv_h(:,2),uv_h(:,3),...
%             M_data,'nearest','none');
    sInt     = scatteredInterpolant(uv_h(:,1),uv_h(:,2),uv_h(:,3),...
        M_data,'linear','none');
    temp_cart = round(sInt(x,y,z));
    temp_cart(isnan(temp_cart))=0 ;
    temp_cart = logical(temp_cart);

    % Remove values outside of the boundary
    sz = size(temp_stack);
    bwMask = poly2mask(uv2(b,1),uv2(b,2),sz(1),sz(2));
    perim  = bwperim(bwMask);
    temp_stack(~repmat(bwMask,[1,1,sz(3)])) = NaN  ;
    
    % Minor clean-up on uv cartilage mask
    temp_cart( ~repmat(bwMask,[1,1,sz(3)])) = false;
    temp_cart(isnan(temp_cart ))=false;
    temp_cart = imclose(temp_cart,strel('disk',3));
    temp_cart(isnan(temp_stack))=false;

    % Directional dilation to fill in missing cartilage pixels
    se = false(3,3,3); se(2,2,2:3) = true;
    w = Inf; count = 0;
    thresh1 = temp_stack > thresh(1); thresh2 = temp_stack < thresh(2);
    while (w > 10 || w < 0) && count < 10
        temp = imdilate(temp_cart,se) & thresh1 & thresh2;
        w = nnz(temp) - nnz(temp_cart);
        if exist('op','var') && op, fprintf([num2str(w),' ']); end
        temp_cart = temp;
        count = count+1;
    end
    if exist('op','var') && op, fprintf('\n'); end

    % Perform some extra cartilage filling
    temp_cart = imclose(temp_cart,strel('disk',3));
    temp_cart = imfill(temp_cart,'holes');
    temp_cart(isnan(temp_stack))=false;
    
    % Remove mask voxels in subchondral bone
    temp_cart(:,:,1:floor(sz(3)/2-1)) = false;

    % Purify and wipe out minor features
    threshV   = round(max(regionprops3(temp_cart,'Volume').Volume) / 20);
    try
        temp_cart = bwareaopen(temp_cart,threshV);
    catch
    end
    temp_cart = imopen(temp_cart,strel('disk',5));
    temp_cart = bwareaopen(temp_cart,threshV);

    % Flip for easier viewing
    temp_stack = flip(temp_stack,3);
    temp_cart  = flip(temp_cart ,3);

    % Write to cell
    uv_stack{jj} = temp_stack;
    uv_cart{ jj} = temp_cart ;
    bwm{     jj} = bwMask    ;

end