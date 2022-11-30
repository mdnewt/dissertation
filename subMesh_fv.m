function [f,v,n] = subMesh_fv(f,v,v_roi,n)
%[ submesh, vertmap, facemap ] = subMesh( mesh, face_roi )
%   This is an altered version of the original subMesh function from the
%   matlabmesh toolbox. In this version, it is not necessary to create a
%   'mesh' object, but rather faces and vertices from isosurface can be
%   directly supplied.

% Which vertices are included?
vidx = find(v_roi);
newv = v(v_roi,:) ;

% Which full faces are defined by included vertices?
f_roi = all(ismember(f,vidx),2);
newf  = f(f_roi,:)             ;

% Are there any extraneous vertices that are not part of a face?
vtest = ismember(vidx,newf);
newv  = newv(vtest,:);
vidx  = vidx(vtest)  ;

if exist('n','var'), n = n(v_roi,:); end

% % Rewrite face indices - Old slow code
% newf2 = zeros(size(newf));
% for ii = 1:length(vidx)
%     fmap = newf==vidx(ii);
%     newf2(fmap) = ii     ;
% end

% Rewrite face indices - Fast new code
sz = size(newf);
I  = reshape(newf',[numel(newf),1]);
[~,~,I2] = unique(I);
newf2 = reshape(I2,sz(2:-1:1))';

f = newf2;
v = newv ;