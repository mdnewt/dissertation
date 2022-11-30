function [mask,xyz,R] = bwtrim(mask,buffer,R)

 if ~exist('buffer','var'), buffer = 0; end
 
 sz     = size(mask)            ;
[x,y,z] = ind2sub(sz,find(mask));

 x = [min(x)-buffer,max(x)+buffer]; x(x > sz(1)) = sz(1); x(x < 1) = 1;
 y = [min(y)-buffer,max(y)+buffer]; y(y > sz(2)) = sz(2); y(y < 1) = 1;
%  z = [min(z)-buffer,max(z)+buffer]; z(z > sz(3)) = sz(3); z(z < 1) = 1;

 mask = mask(x(1):x(2),y(1):y(2),z(1):z(2));

 % Optional generation of updated reference frame
 if nargin==3 && nargout==3
	sz_box = repmat(sz,[2,1])+0.5;
    XYZ = [x',y',z'] + repmat([-0.5;0.5],[1,3]);
    XYZ(XYZ < 0.5   ) = 0.5                 ;
    XYZ(XYZ > sz_box) = sz_box(XYZ > sz_box);
   [X,Y,Z] = intrinsicToWorld(R,[XYZ(3),XYZ(4)],[XYZ(1),XYZ(2)],[XYZ(5),XYZ(6)]);
    temp = R;
    temp.XWorldLimits = X; temp.YWorldLimits = Y; temp.ZWorldLimits = Z;
    temp.ImageSize = diff(XYZ);
    R = temp;
 end
 
 % Optional generation of xyz structure
 if nargout>=2, xyz = struct('x',x,'y',y,'z',z); end
