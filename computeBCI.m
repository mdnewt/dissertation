function [f,v] = computeBCI(stack,cartroi,FT)
 % FT - 0 = Tibia, 1 = Femur
 
 % Get bone mask
 bone = bwareaopen(stack > 3500,500); % Threshold and remove isolated flecks
 bone = imfill(imclose(bone,ones(5,5,5)),'holes'); % Fill pores

 % Get the part of the ROI mask which does not contain bone
 other = imfill(imclose(bone | cartroi,ones(5,5,5)),'holes');
 
 % Purify to largest single blob for femur; largest two blobs for tibia
 if FT, other = purify(other & ~bone);
 else , temp  = purify(other & ~bone); other = temp | purify(other & ~bone & ~temp);
 end

 %% 3) Isolate BCI

 % Compute meshes of bone and its 1 pixel dilation
 sz      = size(bone);
[x,y,z]  = meshgrid(1:sz(2),1:sz(1),1:sz(3));
 if FT, border1 = purify(imdilate(other,ones(3,3,3)) & bone  );
        border2 = purify(imdilate(bone ,ones(3,3,3)) & other );
 else , border1 = purify(imdilate(other,ones(3,3,3)) & bone  );
        border1 = border1 | purify(imdilate(other,ones(3,3,3))...
                & bone  & ~border1);
        border2 = purify(imdilate(bone ,ones(3,3,3)) & other );
        border2 = border2 | purify(imdilate(bone ,ones(3,3,3))...
                & other & ~border2);
 end
[F ,V ]  = isosurface(x,y,z,border1); Vr  = round(V .*2)./2;
[FF,VV]  = isosurface(x,y,z,border2); VVr = round(VV.*2)./2;

 % Isolate shared vertices between bone and its 1 pixel dilation
 Vi = ismember(Vr,VVr,'rows');

 % Create a submesh of the shared vertices
[f,v] = subMesh_fv(F,V,Vi);

%  % Visualize final BCI with bone
%  figure,  patch('Faces',FF,'Vertices',VV,'FaceColor','r','EdgeColor','none');
%  hold on; patch('Faces',f ,'Vertices',v ,'FaceColor','b','EdgeColor','none');
%  legend('Bone','BCI');
%  axis equal
%  title([file,' ',loc]);
